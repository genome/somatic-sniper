#include <bam.h>
#include <assert.h>

/* note: mean_baseQ and mean_mapQ MUST be initialized to 0s or you are going to have a bad day */
void mean_quality_values(
        const bam_pileup1_t* buf,
        int n_reads,
        uint32_t wanted_bases,
        uint32_t mean_baseQ[4],
        uint32_t count_baseQ[4],
        uint32_t mean_mapQ[4],
        uint32_t count_mapQ[4]
        )
{
    int i, j;
    int total = 0;

    memset(&mean_baseQ[0], 0, 4*sizeof(mean_baseQ[0]));
    memset(&mean_mapQ[0], 0, 4*sizeof(mean_mapQ[0])); 
    memset(&count_baseQ[0], 0, 4*sizeof(count_baseQ[0]));  
    memset(&count_mapQ[0], 0, 4*sizeof(count_mapQ[0])); 

    for (i = 0; i < n_reads; ++i) {
        if (buf[i].is_del)
            continue;

        int base = bam1_seqi(bam1_seq(buf[i].b), buf[i].qpos);
        for (j = 0; j < 4; ++j) {
            int value = 1 << j;
            if (base & value & wanted_bases) {
                mean_baseQ[j] += bam1_qual(buf[i].b)[buf[i].qpos];
                ++count_baseQ[j];

                mean_mapQ[j] += buf[i].b->core.qual;
                ++count_mapQ[j];
            }
        }
    }

    for (i = 0; i < 4; ++i) {
        if (count_baseQ[i] > 0)
            mean_baseQ[i] = mean_baseQ[i]/(double)count_baseQ[i] + .499;
        if (count_mapQ[i] > 0)
            mean_mapQ[i] = mean_mapQ[i]/(double)count_mapQ[i] + .499;
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
    if (!need_comma) fputc('-', fh);
}

void print_base_count(FILE* fh, int bases, uint32_t counts[4]) {
    int i;
    int need_comma = 0;
    for (i = 0; i < 4; ++i) {
        int value = 1 << i;
        if (bases & value) {
            if (need_comma) fputc(',', fh);
            /* fprintf(fh, "COUNT_%c=%d",bam_nt16_rev_table[value], counts[i]); */
            fprintf(fh, "%d", counts[i]);
            need_comma = 1;
        }
    }
    if (!need_comma) fputc('-', fh);
}
