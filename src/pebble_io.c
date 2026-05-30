#include "pebble/pebble_io.h"

#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char *chrom;
    int min_start;
    int max_end;
} chrom_extent_t;

typedef struct {
    chrom_extent_t *items;
    size_t count;
    size_t cap;
} chrom_extent_list_t;

static pebble_io_status_t append_extent(chrom_extent_list_t *list, const char *chrom, int start, int end)
{
    for (size_t i = 0; i < list->count; i++) {
        if (strcmp(list->items[i].chrom, chrom) == 0) {
            if (start < list->items[i].min_start) {
                list->items[i].min_start = start;
            }
            if (end > list->items[i].max_end) {
                list->items[i].max_end = end;
            }
            return PEBBLE_IO_OK;
        }
    }

    if (list->count == list->cap) {
        size_t new_cap = list->cap == 0U ? 4U : list->cap * 2U;
        chrom_extent_t *grown = (chrom_extent_t *)realloc(list->items, new_cap * sizeof(chrom_extent_t));
        if (grown == NULL) {
            return PEBBLE_IO_ERR_OOM;
        }
        list->items = grown;
        list->cap = new_cap;
    }

    list->items[list->count].chrom = strdup(chrom);
    if (list->items[list->count].chrom == NULL) {
        return PEBBLE_IO_ERR_OOM;
    }
    list->items[list->count].min_start = start;
    list->items[list->count].max_end = end;
    list->count++;
    return PEBBLE_IO_OK;
}

static void free_extent_list(chrom_extent_list_t *list)
{
    if (list == NULL) {
        return;
    }
    for (size_t i = 0; i < list->count; i++) {
        free(list->items[i].chrom);
    }
    free(list->items);
    list->items = NULL;
    list->count = 0;
    list->cap = 0;
}

static int is_header_line(const char *line)
{
    if (line[0] == '#' || line[0] == '\0') {
        return 1;
    }
    if (strncmp(line, "track", 5) == 0) {
        return 1;
    }
    if (strncmp(line, "browser", 7) == 0) {
        return 1;
    }
    return 0;
}

static pebble_io_status_t parse_interval_line(
    const char *line,
    char *chrom,
    int *start,
    int *end,
    double *value,
    int use_score_field)
{
    char name[256];
    char score_text[64];
    int fields = 0;

    if (use_score_field) {
        fields = sscanf(
            line,
            "%255s\t%d\t%d\t%255[^\t]\t%63[^\t\n]",
            chrom,
            start,
            end,
            name,
            score_text
        );
        if (fields >= 5) {
            *value = atof(score_text);
            return PEBBLE_IO_OK;
        }
        if (fields >= 3) {
            *value = 1.0;
            return PEBBLE_IO_OK;
        }
    } else {
        fields = sscanf(line, "%255s\t%d\t%d\t%lf", chrom, start, end, value);
        if (fields == 4) {
            return PEBBLE_IO_OK;
        }
    }

    (void)name;
    return PEBBLE_IO_ERR_PARSE;
}

static int16_t clamp_to_int16(double value)
{
    if (value >= (double)INT16_MAX) {
        return INT16_MAX;
    }
    if (value <= (double)INT16_MIN) {
        return INT16_MIN;
    }
    return (int16_t)value;
}

static pebble_io_status_t build_batch_from_extents(
    const chrom_extent_list_t *extents,
    pebble_coverage_batch_t *out)
{
    out->items = (pebble_coverage_t *)calloc(extents->count, sizeof(pebble_coverage_t));
    if (out->items == NULL) {
        return PEBBLE_IO_ERR_OOM;
    }
    out->count = extents->count;

    for (size_t i = 0; i < extents->count; i++) {
        size_t length = (size_t)(extents->items[i].max_end - extents->items[i].min_start);
        out->items[i].chrom = strdup(extents->items[i].chrom);
        out->items[i].start_offset = extents->items[i].min_start;
        out->items[i].length = length;
        out->items[i].coverage = (int16_t *)calloc(length, sizeof(int16_t));
        if (out->items[i].chrom == NULL || out->items[i].coverage == NULL) {
            pebble_coverage_batch_free(out);
            return PEBBLE_IO_ERR_OOM;
        }
    }

    return PEBBLE_IO_OK;
}

static int find_chrom_index(const pebble_coverage_batch_t *batch, const char *chrom)
{
    for (size_t i = 0; i < batch->count; i++) {
        if (strcmp(batch->items[i].chrom, chrom) == 0) {
            return (int)i;
        }
    }
    return -1;
}

static pebble_io_status_t fill_batch_from_file(
    const char *path,
    int use_score_field,
    const char *chrom_filter,
    pebble_coverage_batch_t *batch)
{
    FILE *input = fopen(path, "r");
    char line[4096];

    if (input == NULL) {
        return PEBBLE_IO_ERR_IO;
    }

    while (fgets(line, sizeof(line), input) != NULL) {
        char chrom[256];
        int start = 0;
        int end = 0;
        double value = 0.0;
        pebble_io_status_t parse_status;
        int chrom_idx;

        if (is_header_line(line)) {
            continue;
        }

        parse_status = parse_interval_line(
            line,
            chrom,
            &start,
            &end,
            &value,
            use_score_field
        );
        if (parse_status != PEBBLE_IO_OK) {
            fclose(input);
            return parse_status;
        }
        if (end <= start) {
            continue;
        }
        if (chrom_filter != NULL && strcmp(chrom, chrom_filter) != 0) {
            continue;
        }

        chrom_idx = find_chrom_index(batch, chrom);
        if (chrom_idx < 0) {
            continue;
        }

        {
            pebble_coverage_t *entry = &batch->items[(size_t)chrom_idx];
            int16_t clipped = clamp_to_int16(value);
            int local_start = start - entry->start_offset;
            int local_end = end - entry->start_offset;

            if (local_start < 0) {
                local_start = 0;
            }
            if (local_end > (int)entry->length) {
                local_end = (int)entry->length;
            }
            for (int i = local_start; i < local_end; i++) {
                entry->coverage[(size_t)i] = clipped;
            }
        }
    }

    if (ferror(input)) {
        fclose(input);
        return PEBBLE_IO_ERR_IO;
    }

    fclose(input);
    return PEBBLE_IO_OK;
}

static pebble_io_status_t read_interval_file(
    const char *path,
    const char *chrom_filter,
    int use_score_field,
    pebble_coverage_batch_t *out)
{
    FILE *input = fopen(path, "r");
    char line[4096];
    chrom_extent_list_t extents = {0};
    pebble_io_status_t status;

    if (path == NULL || out == NULL) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }

    out->items = NULL;
    out->count = 0;

    if (input == NULL) {
        return PEBBLE_IO_ERR_IO;
    }

    while (fgets(line, sizeof(line), input) != NULL) {
        char chrom[256];
        int start = 0;
        int end = 0;
        double value = 0.0;

        if (is_header_line(line)) {
            continue;
        }

        status = parse_interval_line(
            line,
            chrom,
            &start,
            &end,
            &value,
            use_score_field
        );
        if (status != PEBBLE_IO_OK) {
            fclose(input);
            free_extent_list(&extents);
            return status;
        }
        if (end <= start) {
            continue;
        }
        if (chrom_filter != NULL && strcmp(chrom, chrom_filter) != 0) {
            continue;
        }

        status = append_extent(&extents, chrom, start, end);
        if (status != PEBBLE_IO_OK) {
            fclose(input);
            free_extent_list(&extents);
            return status;
        }
    }

    if (ferror(input)) {
        fclose(input);
        free_extent_list(&extents);
        return PEBBLE_IO_ERR_IO;
    }
    fclose(input);

    if (extents.count == 0U) {
        free_extent_list(&extents);
        return PEBBLE_IO_ERR_PARSE;
    }

    status = build_batch_from_extents(&extents, out);
    free_extent_list(&extents);
    if (status != PEBBLE_IO_OK) {
        return status;
    }

    return fill_batch_from_file(path, use_score_field, chrom_filter, out);
}

void pebble_coverage_batch_free(pebble_coverage_batch_t *batch)
{
    if (batch == NULL) {
        return;
    }
    for (size_t i = 0; i < batch->count; i++) {
        free(batch->items[i].chrom);
        free(batch->items[i].coverage);
    }
    free(batch->items);
    batch->items = NULL;
    batch->count = 0;
}

pebble_io_status_t pebble_read_bedgraph(
    const char *path,
    const char *chrom_filter,
    pebble_coverage_batch_t *out)
{
    return read_interval_file(path, chrom_filter, 0, out);
}

pebble_io_status_t pebble_read_bed(
    const char *path,
    const char *chrom_filter,
    pebble_coverage_batch_t *out)
{
    return read_interval_file(path, chrom_filter, 1, out);
}

pebble_io_status_t pebble_read_coverage(
    const char *path,
    pebble_input_format_t format,
    const char *chrom_filter,
    pebble_coverage_batch_t *out)
{
    switch (format) {
    case PEBBLE_FORMAT_BEDGRAPH:
        return pebble_read_bedgraph(path, chrom_filter, out);
    case PEBBLE_FORMAT_BED:
        return pebble_read_bed(path, chrom_filter, out);
    default:
        return PEBBLE_IO_ERR_UNSUPPORTED;
    }
}

pebble_io_status_t pebble_write_bedgraph(
    FILE *out,
    const char *chrom,
    int start_offset,
    const pebble_config_t *config,
    const double *values,
    size_t value_count)
{
    if (out == NULL || chrom == NULL || config == NULL || values == NULL) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }

    for (size_t idx = 0; idx < value_count; idx++) {
        int center = (int)(idx * (size_t)config->step_size) + (config->window_size / 2) + start_offset;
        int start = center - (config->step_size / 2);
        int end = center + (config->step_size / 2);

        if (fprintf(out, "%s\t%d\t%d\t%.2f\n", chrom, start, end, values[idx]) < 0) {
            return PEBBLE_IO_ERR_IO;
        }
    }

    return PEBBLE_IO_OK;
}

pebble_io_status_t pebble_write_bedgraph_file(
    const char *path,
    const char *chrom,
    int start_offset,
    const pebble_config_t *config,
    const double *values,
    size_t value_count)
{
    FILE *out = stdout;

    if (path != NULL) {
        out = fopen(path, "w");
        if (out == NULL) {
            return PEBBLE_IO_ERR_IO;
        }
    }

    {
        pebble_io_status_t status = pebble_write_bedgraph(
            out,
            chrom,
            start_offset,
            config,
            values,
            value_count
        );
        if (path != NULL) {
            if (fclose(out) != 0 && status == PEBBLE_IO_OK) {
                status = PEBBLE_IO_ERR_IO;
            }
        }
        return status;
    }
}

const char *pebble_io_status_string(pebble_io_status_t status)
{
    switch (status) {
    case PEBBLE_IO_OK:
        return "ok";
    case PEBBLE_IO_ERR_INVALID_ARG:
        return "invalid argument";
    case PEBBLE_IO_ERR_IO:
        return "io error";
    case PEBBLE_IO_ERR_PARSE:
        return "parse error";
    case PEBBLE_IO_ERR_OOM:
        return "out of memory";
    case PEBBLE_IO_ERR_UNSUPPORTED:
        return "unsupported format";
    default:
        return "unknown error";
    }
}
