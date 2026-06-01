#include "pebble/pebble_bigwig.h"

#include "bigWig.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PEBBLE_BIGWIG_MAX_ZOOMS 10

struct pebble_bigwig_writer {
    bigWigFile_t *fp;
    char *current_chrom;
    int has_block;
};

static char *trim_inplace(char *text)
{
    char *end;

    if (text == NULL) {
        return NULL;
    }

    while (*text != '\0' && isspace((unsigned char)*text)) {
        text++;
    }
    if (*text == '\0') {
        return text;
    }

    end = text + strlen(text) - 1;
    while (end > text && isspace((unsigned char)*end)) {
        *end = '\0';
        end--;
    }

    return text;
}

static pebble_io_status_t parse_chrom_sizes_line(
    char *line,
    char **name_out,
    uint32_t *length_out)
{
    char *tab;
    char *name;
    char *length_text;
    unsigned long length;
    char *end;

    *name_out = NULL;
    line = trim_inplace(line);
    if (line[0] == '\0' || line[0] == '#') {
        return PEBBLE_IO_OK;
    }

    tab = strchr(line, '\t');
    if (tab == NULL) {
        tab = strchr(line, ' ');
    }
    if (tab == NULL) {
        return PEBBLE_IO_ERR_PARSE;
    }

    *tab = '\0';
    name = trim_inplace(line);
    length_text = trim_inplace(tab + 1);
    if (name[0] == '\0' || length_text[0] == '\0') {
        return PEBBLE_IO_ERR_PARSE;
    }

    length = strtoul(length_text, &end, 10);
    if (end == length_text || *end != '\0' || length == 0UL || length > UINT32_MAX) {
        return PEBBLE_IO_ERR_PARSE;
    }

    *name_out = strdup(name);
    if (*name_out == NULL) {
        return PEBBLE_IO_ERR_OOM;
    }
    *length_out = (uint32_t)length;
    return PEBBLE_IO_OK;
}

void pebble_chrom_sizes_free(pebble_chrom_sizes_t *sizes)
{
    size_t i;

    if (sizes == NULL) {
        return;
    }

    if (sizes->names != NULL) {
        for (i = 0; i < sizes->count; i++) {
            free(sizes->names[i]);
        }
        free(sizes->names);
    }
    free(sizes->lengths);
    sizes->names = NULL;
    sizes->lengths = NULL;
    sizes->count = 0;
}

pebble_io_status_t pebble_read_chrom_sizes(const char *path, pebble_chrom_sizes_t *out)
{
    FILE *input;
    char line[4096];
    size_t capacity = 0;
    pebble_io_status_t status = PEBBLE_IO_OK;

    if (path == NULL || out == NULL) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }

    out->names = NULL;
    out->lengths = NULL;
    out->count = 0;

    input = fopen(path, "r");
    if (input == NULL) {
        return PEBBLE_IO_ERR_IO;
    }

    while (fgets(line, sizeof(line), input) != NULL) {
        char *name = NULL;
        uint32_t length = 0;
        char **names;
        uint32_t *lengths;

        status = parse_chrom_sizes_line(line, &name, &length);
        if (status != PEBBLE_IO_OK) {
            if (name != NULL) {
                free(name);
            }
            if (status == PEBBLE_IO_OK) {
                continue;
            }
            break;
        }
        if (name == NULL) {
            continue;
        }

        if (out->count == capacity) {
            size_t new_capacity = capacity == 0U ? 16U : capacity * 2U;
            names = (char **)realloc(out->names, new_capacity * sizeof(char *));
            lengths = (uint32_t *)realloc(out->lengths, new_capacity * sizeof(uint32_t));
            if (names == NULL || lengths == NULL) {
                free(name);
                free(names);
                free(lengths);
                status = PEBBLE_IO_ERR_OOM;
                break;
            }
            out->names = names;
            out->lengths = lengths;
            capacity = new_capacity;
        }

        out->names[out->count] = name;
        out->lengths[out->count] = length;
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
        pebble_chrom_sizes_free(out);
    }
    return status;
}

int pebble_chrom_sizes_contains(const pebble_chrom_sizes_t *sizes, const char *chrom)
{
    size_t i;

    if (sizes == NULL || chrom == NULL) {
        return 0;
    }

    for (i = 0; i < sizes->count; i++) {
        if (strcmp(sizes->names[i], chrom) == 0) {
            return 1;
        }
    }
    return 0;
}

pebble_output_format_t pebble_infer_output_format(const char *path)
{
    size_t len;

    if (path == NULL) {
        return PEBBLE_OUTPUT_BEDGRAPH;
    }

    len = strlen(path);
    if (len >= 3 && strcmp(path + len - 3, ".bw") == 0) {
        return PEBBLE_OUTPUT_BIGWIG;
    }
    if (len >= 8 && strcmp(path + len - 8, ".bigwig") == 0) {
        return PEBBLE_OUTPUT_BIGWIG;
    }
    return PEBBLE_OUTPUT_BEDGRAPH;
}

pebble_bigwig_writer_t *pebble_bigwig_writer_create(
    const char *path,
    const pebble_chrom_sizes_t *sizes)
{
    pebble_bigwig_writer_t *writer;

    if (path == NULL || sizes == NULL || sizes->count == 0U) {
        return NULL;
    }

    writer = (pebble_bigwig_writer_t *)calloc(1, sizeof(*writer));
    if (writer == NULL) {
        return NULL;
    }

    if (bwInit(1U << 17) != 0) {
        free(writer);
        return NULL;
    }

    writer->fp = bwOpen(path, NULL, "w");
    if (writer->fp == NULL) {
        bwCleanup();
        free(writer);
        return NULL;
    }

    if (bwCreateHdr(writer->fp, PEBBLE_BIGWIG_MAX_ZOOMS) != 0) {
        bwClose(writer->fp);
        bwCleanup();
        free(writer);
        return NULL;
    }

    writer->fp->cl = bwCreateChromList(
        (const char *const *)sizes->names,
        sizes->lengths,
        (int64_t)sizes->count
    );
    if (writer->fp->cl == NULL) {
        bwClose(writer->fp);
        bwCleanup();
        free(writer);
        return NULL;
    }

    if (bwWriteHdr(writer->fp) != 0) {
        bwClose(writer->fp);
        bwCleanup();
        free(writer);
        return NULL;
    }

    return writer;
}

pebble_io_status_t pebble_bigwig_write_bedgraph(
    pebble_bigwig_writer_t *writer,
    const char *chrom,
    int start_offset,
    const pebble_config_t *config,
    const double *values,
    size_t value_count)
{
    uint32_t *starts = NULL;
    uint32_t *ends = NULL;
    float *float_values = NULL;
    const char **chrom_ptrs = NULL;
    size_t i;
    int bw_status;
    pebble_io_status_t status = PEBBLE_IO_OK;

    if (writer == NULL || writer->fp == NULL || chrom == NULL || config == NULL || values == NULL) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }
    if (value_count == 0U) {
        return PEBBLE_IO_OK;
    }
    if (value_count > UINT32_MAX) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }

    starts = (uint32_t *)malloc(value_count * sizeof(uint32_t));
    ends = (uint32_t *)malloc(value_count * sizeof(uint32_t));
    float_values = (float *)malloc(value_count * sizeof(float));
    chrom_ptrs = (const char **)malloc(value_count * sizeof(char *));
    if (starts == NULL || ends == NULL || float_values == NULL || chrom_ptrs == NULL) {
        status = PEBBLE_IO_ERR_OOM;
        goto cleanup;
    }

    for (i = 0; i < value_count; i++) {
        int interval_start = 0;
        int interval_end = 0;

        pebble_output_interval(i, start_offset, config, &interval_start, &interval_end);
        if (interval_start < 0 || interval_end < 0) {
            status = PEBBLE_IO_ERR_INVALID_ARG;
            goto cleanup;
        }

        starts[i] = (uint32_t)interval_start;
        ends[i] = (uint32_t)interval_end;
        float_values[i] = (float)pebble_round_coverage(values[i]);
        chrom_ptrs[i] = chrom;
    }

    if (!writer->has_block || writer->current_chrom == NULL || strcmp(writer->current_chrom, chrom) != 0) {
        bw_status = bwAddIntervals(
            writer->fp,
            chrom_ptrs,
            starts,
            ends,
            float_values,
            (uint32_t)value_count
        );
        if (writer->current_chrom != NULL) {
            free(writer->current_chrom);
        }
        writer->current_chrom = strdup(chrom);
        if (writer->current_chrom == NULL) {
            status = PEBBLE_IO_ERR_OOM;
            goto cleanup;
        }
        writer->has_block = 1;
    } else {
        bw_status = bwAppendIntervals(
            writer->fp,
            starts,
            ends,
            float_values,
            (uint32_t)value_count
        );
    }

    if (bw_status != 0) {
        status = PEBBLE_IO_ERR_IO;
    }

cleanup:
    free(starts);
    free(ends);
    free(float_values);
    free(chrom_ptrs);
    return status;
}

static pebble_io_status_t bigwig_write_interval_block(
    pebble_bigwig_writer_t *writer,
    const char *chrom,
    const uint32_t *starts,
    const uint32_t *ends,
    const float *values,
    size_t count)
{
    const char **chrom_ptrs;
    size_t i;
    int bw_status;

    if (count == 0U) {
        return PEBBLE_IO_OK;
    }
    if (count > UINT32_MAX) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }

    chrom_ptrs = (const char **)malloc(count * sizeof(char *));
    if (chrom_ptrs == NULL) {
        return PEBBLE_IO_ERR_OOM;
    }
    for (i = 0; i < count; i++) {
        chrom_ptrs[i] = chrom;
    }

    if (!writer->has_block || writer->current_chrom == NULL || strcmp(writer->current_chrom, chrom) != 0) {
        bw_status = bwAddIntervals(
            writer->fp,
            chrom_ptrs,
            starts,
            ends,
            values,
            (uint32_t)count
        );
        if (writer->current_chrom != NULL) {
            free(writer->current_chrom);
        }
        writer->current_chrom = strdup(chrom);
        if (writer->current_chrom == NULL) {
            free(chrom_ptrs);
            return PEBBLE_IO_ERR_OOM;
        }
        writer->has_block = 1;
    } else {
        bw_status = bwAppendIntervals(
            writer->fp,
            starts,
            ends,
            values,
            (uint32_t)count
        );
    }

    free(chrom_ptrs);
    return bw_status == 0 ? PEBBLE_IO_OK : PEBBLE_IO_ERR_IO;
}

pebble_io_status_t pebble_bigwig_write_records(
    pebble_bigwig_writer_t *writer,
    const pebble_bedgraph_record_t *records,
    size_t count)
{
    size_t i = 0;
    pebble_io_status_t status = PEBBLE_IO_OK;

    if (writer == NULL || writer->fp == NULL || records == NULL) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }

    while (i < count && status == PEBBLE_IO_OK) {
        const char *chrom = records[i].chrom;
        size_t j = i + 1U;

        while (j < count && strcmp(records[j].chrom, chrom) == 0) {
            j++;
        }

        {
            size_t block_count = j - i;
            uint32_t *starts = (uint32_t *)malloc(block_count * sizeof(uint32_t));
            uint32_t *ends = (uint32_t *)malloc(block_count * sizeof(uint32_t));
            float *values = (float *)malloc(block_count * sizeof(float));
            size_t k;

            if (starts == NULL || ends == NULL || values == NULL) {
                free(starts);
                free(ends);
                free(values);
                return PEBBLE_IO_ERR_OOM;
            }

            for (k = 0; k < block_count; k++) {
                const pebble_bedgraph_record_t *rec = &records[i + k];

                if (rec->start < 0 || rec->end < 0 || rec->end <= rec->start) {
                    free(starts);
                    free(ends);
                    free(values);
                    return PEBBLE_IO_ERR_INVALID_ARG;
                }
                starts[k] = (uint32_t)rec->start;
                ends[k] = (uint32_t)rec->end;
                values[k] = (float)rec->value;
            }

            status = bigwig_write_interval_block(writer, chrom, starts, ends, values, block_count);
            free(starts);
            free(ends);
            free(values);
        }

        i = j;
    }

    return status;
}

pebble_io_status_t pebble_bigwig_writer_close(pebble_bigwig_writer_t *writer)
{
    pebble_io_status_t status = PEBBLE_IO_OK;

    if (writer == NULL) {
        return PEBBLE_IO_ERR_INVALID_ARG;
    }

    if (writer->fp != NULL) {
        bwClose(writer->fp);
        writer->fp = NULL;
    }

    bwCleanup();
    free(writer->current_chrom);
    free(writer);
    return status;
}
