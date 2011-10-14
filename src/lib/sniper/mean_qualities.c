#include "mean_qualities.h"

#include <bam.h>
#include <assert.h>

void mean_quality_values(
        const bam_pileup1_t* buf,
        int n_reads,
        int ref_base,
        uint32_t wanted_bases,
        uint32_t mean_baseQ[4],
        uint32_t mean_mapQ[4],
        uint32_t base_occ[4],
        dp4_t *dp4
        )
{
    int i, j;
    int base;

    memset(&mean_baseQ[0], 0, 4*sizeof(mean_baseQ[0]));
    memset(&mean_mapQ[0], 0, 4*sizeof(mean_mapQ[0]));
    memset(&base_occ[0], 0, 4*sizeof(base_occ[0]));
    memset(&base_occ[0], 0, 4*sizeof(base_occ[0]));
    memset(dp4, 0, sizeof(*dp4));

    for (i = 0; i < n_reads; ++i) {
        if (buf[i].is_del || buf[i].b->core.flag&BAM_FUNMAP)
            continue;

        base = bam1_seqi(bam1_seq(buf[i].b), buf[i].qpos);
        if (base == ref_base)
            ++dp4->supporting_ref[ bam1_strand(buf[i].b) ];
        else
            ++dp4->variants[ bam1_strand(buf[i].b) ];

        for (j = 0; j < 4; ++j) {
            int value = 1 << j;
            if (base & value & wanted_bases) {
                mean_baseQ[j] += bam1_qual(buf[i].b)[buf[i].qpos];
                mean_mapQ[j] += buf[i].b->core.qual;
                ++base_occ[j];
            }
        }
    }

    for (i = 0; i < 4; ++i) {
        if (base_occ[i] > 0) {
            mean_baseQ[i] = mean_baseQ[i]/(double)base_occ[i] + .499;
            mean_mapQ[i] = mean_mapQ[i]/(double)base_occ[i] + .499;
        }
    }
}

void print_mean_quality_values(FILE* fh, int bases, uint32_t values[4]) {
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

void print_base_count(FILE* fh, int bases, uint32_t counts[4]) {
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

void print_dp4(FILE* fh, const dp4_t* dp4) {
    fprintf(fh, "%d,%d,%d,%d",
        dp4->supporting_ref[0], dp4->supporting_ref[1],
        dp4->variants[0], dp4->variants[1]
    );
}
