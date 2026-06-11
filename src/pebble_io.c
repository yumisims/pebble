#include "pebble/pebble_io.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
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

static int copy_field(const char *start, const char *end, char *out, size_t out_cap)
{
    size_t len;

    if (start == NULL || end == NULL || out == NULL || out_cap == 0U) {
        return -1;
    }
    if (end < start) {
        return -1;
    }

    len = (size_t)(end - start);
    if (len >= out_cap) {
        return -1;
    }

    memcpy(out, start, len);
    out[len] = '\0';
    return 0;
}

static int parse_int_field(const char *text, int *out)
{
    char *end = NULL;
    long value;

    if (text == NULL || out == NULL || text[0] == '\0') {
        return -1;
    }

    errno = 0;
    value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        return -1;
    }
    if (value < (long)INT_MIN || value > (long)INT_MAX) {
        return -1;
    }

    *out = (int)value;
    return 0;
}

static int parse_double_field(const char *text, double *out)
{
    char *end = NULL;

    if (text == NULL || out == NULL || text[0] == '\0') {
        return -1;
    }

    errno = 0;
    *out = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0') {
        return -1;
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
    char fields[6][256];
    size_t field_count = 0U;
    const char *cursor = line;
    const char *field_end;

    while (field_count < 6U) {
        field_end = strchr(cursor, '\t');
        if (field_end == NULL) {
            size_t tail = strcspn(cursor, "\r\n");
            if (copy_field(cursor, cursor + (ptrdiff_t)tail, fields[field_count], sizeof(fields[0])) != 0) {
                return PEBBLE_IO_ERR_PARSE;
            }
            if (fields[field_count][0] != '\0') {
                field_count++;
            }
            break;
        }

        if (copy_field(cursor, field_end, fields[field_count], sizeof(fields[0])) != 0) {
            return PEBBLE_IO_ERR_PARSE;
        }
        field_count++;
        cursor = field_end + 1;
    }

    if (field_count < 3U) {
        return PEBBLE_IO_ERR_PARSE;
    }
    if (parse_int_field(fields[1], start) != 0 || parse_int_field(fields[2], end) != 0) {
        return PEBBLE_IO_ERR_PARSE;
    }

    if (use_score_field) {
        if (field_count >= 5U) {
            if (parse_double_field(fields[4], value) != 0) {
                return PEBBLE_IO_ERR_PARSE;
            }
            memcpy(chrom, fields[0], 256U);
            return PEBBLE_IO_OK;
        }
        if (field_count == 4U) {
            if (parse_double_field(fields[3], value) != 0) {
                return PEBBLE_IO_ERR_PARSE;
            }
            memcpy(chrom, fields[0], 256U);
            return PEBBLE_IO_OK;
        }
        *value = 1.0;
        memcpy(chrom, fields[0], 256U);
        return PEBBLE_IO_OK;
    }

    if (field_count >= 4U) {
        if (parse_double_field(fields[3], value) != 0) {
            return PEBBLE_IO_ERR_PARSE;
        }
        memcpy(chrom, fields[0], 256U);
        return PEBBLE_IO_OK;
    }

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
        size_t length;

        if (extents->items[i].max_end <= 0) {
            pebble_coverage_batch_free(out);
            return PEBBLE_IO_ERR_PARSE;
        }
        if ((size_t)extents->items[i].max_end > SIZE_MAX / sizeof(int16_t)) {
            pebble_coverage_batch_free(out);
            return PEBBLE_IO_ERR_OOM;
        }
        length = (size_t)extents->items[i].max_end;
        out->items[i].chrom = strdup(extents->items[i].chrom);
        out->items[i].start_offset = 0;
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
    int *prev_end = NULL;
    int16_t *prev_value = NULL;

    if (input == NULL) {
        return PEBBLE_IO_ERR_IO;
    }

    if (batch->count > 0U) {
        prev_end = (int *)calloc(batch->count, sizeof(int));
        prev_value = (int16_t *)calloc(batch->count, sizeof(int16_t));
        if (prev_end == NULL || prev_value == NULL) {
            free(prev_end);
            free(prev_value);
            fclose(input);
            return PEBBLE_IO_ERR_OOM;
        }
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
            free(prev_end);
            free(prev_value);
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
            int local_start = start;
            int local_end = end;
            int16_t gap_fill;
            int i;

            if (local_end > (int)entry->length) {
                local_end = (int)entry->length;
            }
            if (local_start < 0) {
                local_start = 0;
            }

            if (start > prev_end[chrom_idx]) {
                gap_fill = (prev_end[chrom_idx] == 0) ? clipped : prev_value[chrom_idx];
                for (i = prev_end[chrom_idx]; i < start && i < (int)entry->length; i++) {
                    entry->coverage[(size_t)i] = gap_fill;
                }
            }

            for (i = local_start; i < local_end; i++) {
                entry->coverage[(size_t)i] = clipped;
            }

            if (end > prev_end[chrom_idx]) {
                prev_end[chrom_idx] = end;
            }
            prev_value[chrom_idx] = clipped;
        }
    }

    free(prev_end);
    free(prev_value);

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

    status = fill_batch_from_file(path, use_score_field, chrom_filter, out);
    if (status != PEBBLE_IO_OK) {
        pebble_coverage_batch_free(out);
    }
    return status;
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

static int compare_coverage_chrom(const void *left, const void *right)
{
    const pebble_coverage_t *a = (const pebble_coverage_t *)left;
    const pebble_coverage_t *b = (const pebble_coverage_t *)right;
    int chrom_cmp = strcmp(a->chrom, b->chrom);

    if (chrom_cmp != 0) {
        return chrom_cmp;
    }
    if (a->start_offset < b->start_offset) {
        return -1;
    }
    if (a->start_offset > b->start_offset) {
        return 1;
    }
    return 0;
}

void pebble_coverage_batch_sort(pebble_coverage_batch_t *batch)
{
    if (batch == NULL || batch->items == NULL || batch->count <= 1U) {
        return;
    }

    qsort(batch->items, batch->count, sizeof(pebble_coverage_t), compare_coverage_chrom);
}

static int compare_bedgraph_record(const void *left, const void *right)
{
    const pebble_bedgraph_record_t *a = (const pebble_bedgraph_record_t *)left;
    const pebble_bedgraph_record_t *b = (const pebble_bedgraph_record_t *)right;
    int chrom_cmp = strcmp(a->chrom, b->chrom);

    if (chrom_cmp != 0) {
        return chrom_cmp;
    }
    if (a->start < b->start) {
        return -1;
    }
    if (a->start > b->start) {
        return 1;
    }
    return 0;
}

void pebble_bedgraph_batch_free(pebble_bedgraph_batch_t *batch)
{
    size_t i;

    if (batch == NULL) {
        return;
    }
    if (batch->items != NULL) {
        for (i = 0; i < batch->count; i++) {
            free(batch->items[i].chrom);
        }
        free(batch->items);
    }
    batch->items = NULL;
    batch->count = 0;
}

void pebble_bedgraph_batch_sort(pebble_bedgraph_batch_t *batch)
{
    if (batch == NULL || batch->items == NULL || batch->count <= 1U) {
        return;
    }

    qsort(batch->items, batch->count, sizeof(pebble_bedgraph_record_t), compare_bedgraph_record);
}

pebble_io_status_t pebble_read_interval_records(
    const char *path,
    const char *chrom_filter,
    int use_score_field,
    pebble_bedgraph_batch_t *out)
{
    FILE *input;
    char line[4096];
    size_t cap = 0;
    pebble_io_status_t status = PEBBLE_IO_OK;

    if (path == NULL || out == NULL) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }

    out->items = NULL;
    out->count = 0;

    input = fopen(path, "r");
    if (input == NULL) {
        return PEBBLE_IO_ERR_IO;
    }

    while (fgets(line, sizeof(line), input) != NULL) {
        char chrom[256];
        int start = 0;
        int end = 0;
        double value = 0.0;
        pebble_bedgraph_record_t *items;
        char *chrom_copy;

        if (is_header_line(line)) {
            continue;
        }

        status = parse_interval_line(line, chrom, &start, &end, &value, use_score_field);
        if (status != PEBBLE_IO_OK) {
            break;
        }
        if (end <= start) {
            continue;
        }
        if (chrom_filter != NULL && strcmp(chrom, chrom_filter) != 0) {
            continue;
        }

        if (out->count == cap) {
            size_t new_cap = cap == 0U ? 64U : cap * 2U;
            items = (pebble_bedgraph_record_t *)realloc(out->items, new_cap * sizeof(pebble_bedgraph_record_t));
            if (items == NULL) {
                status = PEBBLE_IO_ERR_OOM;
                break;
            }
            out->items = items;
            cap = new_cap;
        }

        chrom_copy = strdup(chrom);
        if (chrom_copy == NULL) {
            status = PEBBLE_IO_ERR_OOM;
            break;
        }

        out->items[out->count].chrom = chrom_copy;
        out->items[out->count].start = start;
        out->items[out->count].end = end;
        out->items[out->count].value = value;
        out->count++;
    }

    if (ferror(input)) {
        status = PEBBLE_IO_ERR_IO;
    }
    if (status == PEBBLE_IO_OK && out->count == 0U) {
        status = PEBBLE_IO_ERR_PARSE;
    }

    fclose(input);
    if (status != PEBBLE_IO_OK) {
        pebble_bedgraph_batch_free(out);
    }
    return status;
}

void pebble_bedgraph_batch_extend_to_zero(pebble_bedgraph_batch_t *batch)
{
    size_t i = 0;

    if (batch == NULL || batch->items == NULL) {
        return;
    }

    while (i < batch->count) {
        size_t j = i + 1U;

        while (j < batch->count && strcmp(batch->items[j].chrom, batch->items[i].chrom) == 0) {
            j++;
        }
        if (batch->items[i].start > 0) {
            batch->items[i].start = 0;
        }
        i = j;
    }
}

pebble_io_status_t pebble_read_bedgraph_records(
    const char *path,
    const char *chrom_filter,
    pebble_bedgraph_batch_t *out)
{
    return pebble_read_interval_records(path, chrom_filter, 0, out);
}

pebble_io_status_t pebble_write_bedgraph_records(
    FILE *out,
    const pebble_bedgraph_record_t *records,
    size_t count)
{
    size_t i;

    if (out == NULL || records == NULL) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }

    for (i = 0; i < count; i++) {
        if (
            fprintf(
                out,
                "%s\t%d\t%d\t%d\n",
                records[i].chrom,
                records[i].start,
                records[i].end,
                pebble_round_coverage(records[i].value)
            ) < 0
        ) {
            return PEBBLE_IO_ERR_IO;
        }
    }

    return PEBBLE_IO_OK;
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

void pebble_output_interval(
    size_t step_index,
    int start_offset,
    const pebble_config_t *config,
    int *interval_start,
    int *interval_end)
{
    int window_start;

    if (config == NULL || interval_start == NULL || interval_end == NULL) {
        return;
    }

    window_start = (int)(step_index * (size_t)config->step_size) + start_offset;
    *interval_start = window_start;
    *interval_end = window_start + config->step_size;
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

static int round_coverage(double value)
{
    if (value >= 0.0) {
        return (int)(value + 0.5);
    }
    return (int)(value - 0.5);
}

int pebble_round_coverage(double value)
{
    return round_coverage(value);
}

void pebble_smoothed_genome_average_add(
    int start_offset,
    const pebble_config_t *config,
    const double *values,
    size_t value_count,
    double *weighted_sum,
    size_t *total_bases)
{
    size_t idx;

    if (config == NULL || values == NULL || weighted_sum == NULL || total_bases == NULL) {
        return;
    }

    for (idx = 0; idx < value_count; idx++) {
        int start = 0;
        int end = 0;
        size_t span;

        pebble_output_interval(idx, start_offset, config, &start, &end);
        span = (size_t)(end - start);
        if (span == 0U) {
            continue;
        }

        *weighted_sum += values[idx] * (double)span;
        *total_bases += span;
    }
}

double pebble_coverage_normalise_factor(double genome_average)
{
    if (genome_average <= 0.0) {
        return 1.0;
    }
    return genome_average / 2.0;
}

void pebble_normalise_smoothed_values(
    double *values,
    size_t value_count,
    double normalise_factor)
{
    size_t idx;

    if (values == NULL || normalise_factor <= 0.0) {
        return;
    }

    for (idx = 0; idx < value_count; idx++) {
        values[idx] /= normalise_factor;
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
        int start = 0;
        int end = 0;

        pebble_output_interval(idx, start_offset, config, &start, &end);

        if (fprintf(out, "%s\t%d\t%d\t%d\n", chrom, start, end, round_coverage(values[idx])) < 0) {
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
