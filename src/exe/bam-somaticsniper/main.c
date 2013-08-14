#include "faidx.h"
#include "khash.h"
#include "kstring.h"
#include "sam.h"
#include "sniper/sniper_maqcns.h"
#include "sniper/somatic_sniper.h"
#include "sniper/output_format.h"
#include "version.h"

#include <bam.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

static const char *_default_normal_sample_id = "NORMAL";
static const char *_default_tumor_sample_id = "TUMOR";
static const char *_default_output_format = "classic";

void version_info() {
    printf("Somatic Sniper version (%s) (commit %s)", __g_prog_version, __g_commit_hash);
    if (strlen(__g_build_type) > 0)
        printf(" (%s)", __g_build_type);
    printf("\n");
}

void usage(const char* progname, pu_data2_t* d) {
    /* we dont like basename(3) */
    const char* pn = strrchr(progname, '/');
    int i;
    int n_formats = n_output_formatters();
    if (pn == NULL)
        pn = progname;
    else
        ++pn;

    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s [options] -f <ref.fasta> <tumor.bam> <normal.bam> <snp_output_file>\n\n", pn);
    fprintf(stderr, "Required Option: \n");
    fprintf(stderr, "        -f FILE   REQUIRED reference sequence in the FASTA format\n\n");
    fprintf(stderr, "Options: \n");
    fprintf(stderr, "        -v        Display version information\n\n");
    fprintf(stderr, "        -q INT    filtering reads with mapping quality less than INT [%d]\n", d->mapQ);
    fprintf(stderr, "        -Q INT    filtering somatic snv output with somatic quality less than  INT [%d]\n", d->min_somatic_qual);
    fprintf(stderr, "        -p FLAG   disable priors in the somatic calculation. Increases sensitivity for solid tumors\n");
    fprintf(stderr, "        -J FLAG   Use prior probabilities accounting for the somatic mutation rate\n");
    fprintf(stderr, "        -s FLOAT  prior probability of a somatic mutation (implies -J) [%f]\n",d->somatic_mutation_rate);
    fprintf(stderr, "        -T FLOAT  theta in maq consensus calling model (for -c/-g) [%f]\n", d->c->theta);
    fprintf(stderr, "        -N INT    number of haplotypes in the sample (for -c/-g) [%d]\n", d->c->n_hap);
    fprintf(stderr, "        -r FLOAT  prior of a difference between two haplotypes (for -c/-g) [%f]\n", d->c->het_rate);
    fprintf(stderr, "        -n STRING normal sample id (for VCF header) [%s]\n", _default_normal_sample_id);
    fprintf(stderr, "        -t STRING tumor sample id (for VCF header) [%s]\n", _default_tumor_sample_id);
    fprintf(stderr, "        -F STRING select output format [%s]\n", _default_output_format);
    fprintf(stderr, "           Available formats:\n");
    for (i = 0; i < n_formats; ++i) {
        fprintf(stderr, "             %s\n", output_formatter_name(i));
    }
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
    int c;
    const char* normal_sample_id = _default_normal_sample_id;
    const char* tumor_sample_id = _default_tumor_sample_id;
    const char *fn_fa = 0;
    pu_data2_t *d = (pu_data2_t*)calloc(1, sizeof(pu_data2_t));
    d->min_somatic_qual=15;
    d->tid = -1; d->mask = BAM_DEF_MASK; d->mapQ = 0;
    d->c = sniper_maqcns_init();
    int use_priors = 1;
    d->use_joint_priors = 0;
    d->somatic_mutation_rate = 0.01;
    const char *output_format = "classic";

    while ((c = getopt(argc, argv, "n:t:vf:T:N:r:I:G:q:Q:pJs:F:")) >= 0) {
        switch (c) {
            case 'f': fn_fa = optarg; break;
            case 'T': d->c->theta = atof(optarg); break;
            case 'N': d->c->n_hap = atoi(optarg); break;
            case 'r': d->c->het_rate = atof(optarg); break;
            case 'q': d->mapQ = atoi(optarg); break;
            case 'Q': d->min_somatic_qual = atoi(optarg); break;
            case 'F': output_format = optarg; break;
            case 'p': use_priors = 0; break;
            case 'J': d->use_joint_priors = 1; break;
            case 's': d->somatic_mutation_rate = atof(optarg); d->use_joint_priors = 1; break;
            case 'v': version_info(); exit(0); break;
            case 't': tumor_sample_id = optarg; break;
            case 'n': normal_sample_id = optarg; break;
            default: fprintf(stderr, "Unrecognizd option '-%c'.\n", c); return 1;
        }
    }

    if (optind == argc) {
        usage(argv[0], d);
        sniper_maqcns_destroy(d->c); free(d);
        return 1;
    }
    if (fn_fa) {
        d->fai = fai_load(fn_fa);
    }
    else {
        fprintf(stderr, "You MUST specify a reference sequence. It isn't optional.\n");
        sniper_maqcns_destroy(d->c);
        free(d);
        exit(1);
    }
    if(d->use_joint_priors) {
        fprintf(stderr,"Using priors accounting for somatic mutation rate. Prior probability of a somatic mutation is %f\n",d->somatic_mutation_rate);
        make_joint_prior(d->somatic_mutation_rate);
    }

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
    /* this will exit if the format name is invalid */
    output_formatter_t fmt = output_formatter_create(output_format, snp_fh);
    d->output_formatter = &fmt;
    if(snp_fh) {
        header_data_t hdr;
        hdr.refseq = fn_fa;
        hdr.normal_sample_id = normal_sample_id;
        hdr.tumor_sample_id = tumor_sample_id;
        d->output_formatter->header_fn(snp_fh, &hdr);
        bam_sspileup_file(fp1, fp2, d->mask, d->mapQ, glf_somatic, d, snp_fh);
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
    free(d->ref); free(d);
    fclose(snp_fh);
    return 0;
}
