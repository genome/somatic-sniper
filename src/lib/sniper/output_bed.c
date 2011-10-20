#include "output_bed.h"

#include <stdio.h>

void output_bed_header(FILE* fh, const header_data_t* h) {
}

void output_bed(FILE *fh, const sniper_output_t *p) {
    fprintf(fh, "%s\t%d\t%d\t%c/%c\t%d\t%d\n",
        p->seq_name,
        p->pos,
        p->pos + 1,
        p->ref_base,
        bam_nt16_rev_table[p->tumor.genotype],
        p->tumor.somatic_score,
        p->tumor.dqstats.total_depth
    );
}
