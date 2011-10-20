#include "output_classic.h"

#include "dqstats.h"
#include <stdio.h>

void output_classic_header(FILE* fh, const header_data_t* h) {
}

void output_classic(FILE *fh, const sniper_output_t *p) {
    fprintf(fh, "%s\t%d\t%c\t%c\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",
        p->seq_name,
        p->pos + 1,
        p->ref_base,
        bam_nt16_rev_table[p->tumor.genotype],
        bam_nt16_rev_table[p->normal.genotype],
        p->tumor.somatic_score,
        p->tumor.consensus_quality,
        p->tumor.variant_allele_quality,
        p->tumor.dqstats.total_mean_mapQ,
        p->normal.consensus_quality,
        p->normal.variant_allele_quality,
        p->normal.dqstats.total_mean_mapQ,
        p->tumor.dqstats.total_depth,
        p->normal.dqstats.total_depth
    );


    /* mean {map,base} quality for tumor */
    print_mean_quality_values(fh, p->ref_base4, p->tumor.dqstats.mean_baseQ);
    fputc('\t', fh);
    print_mean_quality_values(fh, p->ref_base4, p->tumor.dqstats.mean_mapQ);
    fputc('\t', fh);
    print_base_count(fh, p->ref_base4, p->tumor.dqstats.base_occ);
    fputc('\t', fh);
    print_mean_quality_values(fh, ~p->ref_base4 & p->tumor.genotype, p->tumor.dqstats.mean_baseQ);
    fputc('\t', fh);
    print_mean_quality_values(fh, ~p->ref_base4 & p->tumor.genotype, p->tumor.dqstats.mean_mapQ);
    fputc('\t', fh);
    print_base_count(fh, ~p->ref_base4 & p->tumor.genotype, p->tumor.dqstats.base_occ);
    fputc('\t', fh);
    print_dp4(fh, p->tumor.dqstats.dp4);
    fputc('\t', fh);

    /* mean {map,base} quality for normal */
    print_mean_quality_values(fh, p->ref_base4, p->normal.dqstats.mean_baseQ);
    fputc('\t', fh);
    print_mean_quality_values(fh, p->ref_base4, p->normal.dqstats.mean_mapQ);
    fputc('\t', fh);
    print_base_count(fh, p->ref_base4, p->normal.dqstats.base_occ);
    fputc('\t', fh);
    print_mean_quality_values(fh, ~p->ref_base4&p->normal.genotype, p->normal.dqstats.mean_baseQ);
    fputc('\t', fh);
    print_mean_quality_values(fh, ~p->ref_base4&p->normal.genotype, p->normal.dqstats.mean_mapQ);
    fputc('\t', fh);
    print_base_count(fh, ~p->ref_base4&p->normal.genotype, p->normal.dqstats.base_occ);
    fputc('\t', fh);
    print_dp4(fh, p->normal.dqstats.dp4);
    fputc('\n', fh);
}
