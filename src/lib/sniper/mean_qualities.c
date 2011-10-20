#include "mean_qualities.h"

#include <bam.h>
#include <assert.h>

void get_dqstats(
        const bam_pileup1_t* buf,
        int n_reads,
        int ref_base,
        uint32_t wanted_bases,
        dqstats_t *dqs
        )
{
    int i, j;
    int base;

    memset(dqs, 0, sizeof(dqstats_t));

    for (i = 0; i < n_reads; ++i) {
        if (buf[i].is_del || buf[i].b->core.flag&BAM_FUNMAP)
            continue;

        ++dqs->total_depth;
        dqs->total_mean_mapQ += buf[i].b->core.qual;

        base = bam1_seqi(bam1_seq(buf[i].b), buf[i].qpos);
        if (base & ref_base)
            ++dqs->dp4[bam1_strand(buf[i].b)];
        if (base & ~ref_base)
            ++dqs->dp4[2+bam1_strand(buf[i].b)];

        for (j = 0; j < 4; ++j) {
            int value = 1 << j;
            if (base & value & wanted_bases) {
                dqs->mean_baseQ[j] += bam1_qual(buf[i].b)[buf[i].qpos];
                dqs->mean_mapQ[j] += buf[i].b->core.qual;
                ++dqs->base_occ[j];
            }
        }
    }

    for (i = 0; i < 4; ++i) {
        if (dqs->base_occ[i] > 0) {
            dqs->mean_baseQ[i] = dqs->mean_baseQ[i]/(double)dqs->base_occ[i] + .499;
            dqs->mean_mapQ[i] = dqs->mean_mapQ[i]/(double)dqs->base_occ[i] + .499;
        }
    }

    if (dqs->total_depth > 0)
        dqs->total_mean_mapQ = dqs->total_mean_mapQ / (double)dqs->total_depth + .499;
}

void print_mean_quality_values(FILE* fh, int bases, const uint32_t values[4]) {
    int i;
    int need_comma = 0;
    for (i = 0; i < 4; ++i) {
        int value = 1 << i;
        if (bases & value) {
            if (need_comma) fputc(',', fh);
            /* fprintf(fh, "QUAL_%c=%d", bam_nt16_rev_table[value], values[i]); */
            fprintf(fh, "%d", values[i]);
            need_comma = 1;
        }
    }
    if (!need_comma) fputc('0', fh);
}

void print_base_count(FILE* fh, int bases, const uint32_t counts[4]) {
    int i;
    int need_comma = 0;
    for (i = 0; i < 4; ++i) {
        int value = 1 << i;
        if (bases & value) {
            if (need_comma) fputc(',', fh);
            /* fprintf(fh, "COUNT_%c=%d:",bam_nt16_rev_table[value], counts[i]); */
            fprintf(fh, "%d", counts[i]);
            need_comma = 1;
        }
    }
    if (!need_comma) fputc('0', fh);
}

void print_dp4(FILE* fh, const uint32_t dp4[4]) {
    fprintf(fh, "%d,%d,%d,%d", dp4[0], dp4[1], dp4[2], dp4[3]);
}
