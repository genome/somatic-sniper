#include "output_format.h"
#include "output_bed.h"
#include "output_classic.h"
#include "output_vcf.h"

#include <stdio.h>
#include <stdlib.h>

const static struct {
    const char *name;
    output_header_fn header_fn;
    output_fn output_fn;
} _formats[] = {
    { "classic", &output_classic_header, &output_classic },
    { "vcf", &output_vcf_header, &output_vcf },
    { "bed", &output_bed_header, &output_bed }
};
const static uint32_t _n_formats = sizeof(_formats)/sizeof(_formats[0]);

output_formatter_t output_formatter_create(const char* name, FILE* fh) {
    uint32_t i;
    for (i = 0; i < _n_formats; ++i) {
        if (strcmp(name, _formats[i].name) == 0) {
            output_formatter_t rv;
            rv.fh = fh;
            rv.header_fn = _formats[i].header_fn;
            rv.output_fn = _formats[i].output_fn;
            return rv;
        }
    }
    fprintf(stderr, "unknown output format: '%s'. Abort!\n", name);
    exit(1);
}

void output_formatter_write(const output_formatter_t *formatter, const sniper_output_t *p) {
    formatter->output_fn(formatter->fh, p);
    fflush(formatter->fh);
}

uint32_t n_output_formatters() {
    return _n_formats;
}

const char* output_formatter_name(uint32_t idx) {
    if (idx >= _n_formats) {
        fprintf(stderr, "Request for output formatter '%u' is out of range\n", idx);
        exit(1);
    }
    return _formats[idx].name;
}
