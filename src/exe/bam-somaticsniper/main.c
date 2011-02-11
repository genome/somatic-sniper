#include "faidx.h"
#include "khash.h"
#include "kstring.h"
#include "sam.h"
#include "sniper/mean_qualities.h"
#include "sniper/sniper_maqcns.h"
#include "sniper/somatic_sniper.h"

#include <bam.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

void usage(const char* progname, pu_data2_t* d) {
    /* we dont like basename(3) */
    const char* pn = strrchr(progname, '/');
    if (pn == NULL)
        pn = progname;
    else
        ++pn;

    fprintf(stderr, "\n");
    fprintf(stderr, "%s [options] -f <ref.fasta> <tumor.bam> <normal.bam> <snp_output_file>\n\n", pn);
    fprintf(stderr, "Required Option: \n");
    fprintf(stderr, "        -f FILE   REQUIRED reference sequence in the FASTA format\n\n");
    fprintf(stderr, "Options: \n");
    fprintf(stderr, "        -q INT    filtering reads with mapping quality less than INT [%d]\n", d->mapQ);
    fprintf(stderr, "        -Q INT    filtering somatic snp output with somatic quality less than  INT [15]\n");
    fprintf(stderr, "        -p FLAG   disable priors in the somatic calculation. Increases sensitivity for solid tumors\n");
    fprintf(stderr, "        -T FLOAT  theta in maq consensus calling model (for -c/-g) [%f]\n", d->c->theta);
    fprintf(stderr, "        -N INT    number of haplotypes in the sample (for -c/-g) [%d]\n", d->c->n_hap);

    fprintf(stderr, "        -r FLOAT  prior of a difference between two haplotypes (for -c/-g) [%f]\n", d->c->het_rate);
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
    int c;
    char *fn_fa = 0;
    pu_data2_t *d = (pu_data2_t*)calloc(1, sizeof(pu_data2_t));
    d->min_somatic_qual=0;
    d->tid = -1; d->mask = BAM_DEF_MASK; d->mapQ = 0;
    d->c = sniper_maqcns_init();
    int use_priors = 1;
    while ((c = getopt(argc, argv, "f:T:N:r:I:G:q:Q:p")) >= 0) {
        switch (c) {
            case 'f': fn_fa = strdup(optarg); break;
            case 'T': d->c->theta = atof(optarg); break;
            case 'N': d->c->n_hap = atoi(optarg); break;
            case 'r': d->c->het_rate = atoi(optarg); break;
            case 'q': d->mapQ = atoi(optarg); break;
            case 'Q': d->min_somatic_qual = atoi(optarg); break;
            case 'p': use_priors = 0; break;
            default: fprintf(stderr, "Unrecognizd option '-%c'.\n", c); return 1;
        }
    }
    if (optind == argc) {
        usage(argv[0], d);
        free(fn_fa); sniper_maqcns_destroy(d->c); free(d->ido); free(d);
        return 1;
    }
    if (fn_fa) {
        d->fai = fai_load(fn_fa);
    }
    else {
        fprintf(stderr, "You MUST specify a reference sequence. It isn't optional.\n");
        sniper_maqcns_destroy(d->c);
        free(d->ido);
        free(d);
        exit(1);
    }
    free(fn_fa);
    sniper_maqcns_prepare(d->c);
    fprintf(stderr,"Preparing to snipe some somatics\n");
    if(use_priors) {
        fprintf(stderr,"Using prior probabilities\n");
        makeSoloPrior();
    }
    bamFile fp1, fp2;
    qAddTableInit();
    fp1 = (strcmp(argv[optind], "-") == 0)? bam_dopen(fileno(stdin), "r") : bam_open(argv[optind], "r");
    fprintf(stderr, "Normal bam is %s\n", argv[optind+1]);
    fprintf(stderr, "Tumor bam is %s\n", argv[optind]);
    d->h1 = bam_header_read(fp1);
    sam_header_parse_rg(d->h1);
    fp2 = bam_open(argv[optind+1], "r");
    d->h2 = bam_header_read(fp2);
    sam_header_parse_rg(d->h2);
    FILE* snp_fh = fopen(argv[optind+2], "w");
    FILE* indel_fh = NULL; 

    if(snp_fh) {
        //NEED TO ADD IN AN ACTUAL FUNCTION NAME HERE
        bam_sspileup_file(fp1, fp2, d->mask, d->mapQ, glf_somatic, d, snp_fh, indel_fh);
    }
    else {
        fprintf(stderr, "Unable to open snp file!!!!!!!!!\n");
        exit(1);
    }

    bam_close(fp1);
    bam_close(fp2);
    bam_header_destroy(d->h1);
    bam_header_destroy(d->h2);
    if (d->fai) fai_destroy(d->fai);
    sniper_maqcns_destroy(d->c);
    free(d->ido); free(d->ref); free(d);
    fclose(snp_fh);
    if (indel_fh)
        fclose(indel_fh);
    return 0;
}
