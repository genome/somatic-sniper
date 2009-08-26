#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include "sam.h"
#include "faidx.h"
#include "sniper_maqcns.h"
#include "khash.h"
#include "kstring.h"
#include "somatic_sniper.h"

typedef int *indel_list_t;
KHASH_MAP_INIT_INT64(64, indel_list_t)
int get_next_pos(bam_plbuf_t *buf,bamFile fp); 
#define BAM_PLF_SIMPLE     0x01
#define BAM_PLF_CNS        0x02
#define BAM_PLF_INDEL_ONLY 0x04
#define BAM_PLF_GLF        0x08
#define BAM_PLF_VAR_ONLY   0x10
#define HOMO_INDEL1 0
#define HOMO_INDEL2 1
#define HET_INDEL   2

static int qAddTable[1024];
double THETA = 0.001 ;      /* population scaled mutation rate */
int minimum_somatic_qual = 4; //minimum somatic phred score in order to report the site
static int prior[16][10] ;  /* index over reference base, genotype */

#define qAdd(x,y)  (x - qAddTable[512+y-x])

typedef struct {
    bam_header_t *h1;
    bam_header_t *h2;
    sniper_maqcns_t *c;
    sniper_maqindel_opt_t *ido;
    faidx_t *fai;
    khash_t(64) *hash;
    uint32_t format;
    int tid, len, last_pos;
    int mask;
    int mapQ;
    int min_somatic_qual;//for limiting snp calls in somatic sniper
    char *ref;
    glfFile fp; // for glf output only
} pu_data2_t;

int prob_indel(glf3_t* indel, int prior_prob, int likelihood_flag);
sniper_maqindel_ret_t *sniper_somaticindelscore(int n, int pos, const sniper_maqindel_opt_t *mi, const bam_pileup1_t *pl, const char *ref, sniper_maqindel_ret_t *tumor_indel);
char **__bam_get_lines(const char *fn, int *_n);
void bam_init_header_hash(bam_header_t *header);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);


void calculatePosteriors(glf1_t *g, int lkResult[]) {
    unsigned char refBase = g->ref_base;
    int qSum = 255;
    int qMin = 1000;
    int j;

    //Calculate Posteriors
    for (j = 0 ; j < 10 ; ++j) { 
        int x = g->lk[j] + prior[refBase][j];
        qSum = qAdd (x, qSum) ;
        if (x < qMin) qMin = x ;
        lkResult[j] = x ;
    }
    //fprintf(stderr,"qSum = %d\nqMin = %d\n",qSum,qMin);
    //iprintLikelihoods(lkResult);
    for (j = 0 ; j < 10 ; ++j) {
        lkResult[j] -= qSum;
        //lkResult[j] -= qMin;
        if(lkResult[j] > 255) {
            lkResult[j] = 255;
        }
    }

}
static void qAddTableInit (void)
{
    int i ;
    for (i = 0 ; i < 1000 ; ++i)
    { double e = 1 + expPhred(i-512) ;
        qAddTable[i] = -logPhred(e) ;
    }
}

glf3_t *sniper_maqindel2glf(sniper_maqindel_ret_t *r, int n) {
    glf3_t *g3;
    g3 = glf3_init1();

    int het = 3 * n, min;
    min = het;
    if (min > r->gl[0]) min = r->gl[0];
    if (min > r->gl[1]) min = r->gl[1];
    g3->ref_base = 0;
    g3->rtype = GLF3_RTYPE_INDEL;
    memset(g3->lk, 0, 10);
    g3->lk[0] = r->gl[0] - min < 255? r->gl[0] - min : 255;
    g3->lk[1] = r->gl[1] - min < 255? r->gl[1] - min : 255;
    g3->lk[2] = het - min < 255? het - min : 255;
    g3->offset = 0;
    g3->indel_len[0] = r->indel1;
    g3->indel_len[1] = r->indel2;
    g3->min_lk = min < 255? min : 255;
    g3->max_len = (abs(r->indel1) > abs(r->indel2)? abs(r->indel1) : abs(r->indel2)) + 1;
    g3->indel_seq[0] = strdup(r->s[0]+1);
    g3->indel_seq[1] = strdup(r->s[1]+1);
    return g3;
}



static int glf_somatic(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE *snp_fh, FILE *indel_fh) {
    //hacked copy from function gl3_func behavior to get a g with 10 probabilities to do somatic probability calculation    
    pu_data2_t *d = (pu_data2_t*)data;
    sniper_maqindel_ret_t *r = 0;
    if (d->fai && (int)tid != d->tid) {
        free(d->ref);
        d->ref = fai_fetch(d->fai, d->h1->target_name[tid], &d->len);
        d->tid = tid;
    }
    int rb = (d->ref && (int)pos < d->len)? d->ref[pos] : 'N';

    int lkTumor[10], lkNormal[10];

    int qPosteriorSum = 255;

    glf1_t *gTumor =sniper_maqcns_glfgen(n1, pl1, bam_nt16_table[rb], d->c);
    glf1_t *gNormal =sniper_maqcns_glfgen(n2, pl2, bam_nt16_table[rb], d->c);
    //now we have the filled g1,g2 to compare with code from larsons's glfSomatic
    if (rb != 'N' && gTumor->depth > 0 && gNormal->depth > 0) {
        //calculate tumor posteriors
        //fprintf(stderr,"---Tumor---\n");
        //printGLF(&gTumor);

        uint32_t x;
        x = sniper_maqcns_call(n1, pl1, d->c);
        int  ref_q, rb4 = bam_nt16_table[rb];
        if (rb4 != 15 && x>>28 != 15 && x>>28 != rb4) { // a SNP
            ref_q = 0;
            if (rb4 != 15 && x>>28 != 15 && x>>28 != rb4) { // a SNP
                ref_q = ((x>>24&0xf) == rb4)? x>>8&0xff : (x>>8&0xff) + (x&0xff);
                if (ref_q > 255) ref_q = 255;
            }

            calculatePosteriors(gTumor, lkTumor);
            //iprintLikelihoods(lkTumor);
            //fprintf(stderr,"---Normal---\n");
            //printGLF(&gNormal);
            calculatePosteriors(gNormal, lkNormal);
            //iprintLikelihoods(lkNormal);
            int j;
            for(j = 0; j < 10; j++) {
                qPosteriorSum = qAdd(qPosteriorSum,(lkTumor[j] + lkNormal[j]));
            }

            // int result = qAdd(0,-qPosteriorSum);
            if(d->min_somatic_qual <= qPosteriorSum) {  
                fprintf(snp_fh, "%s\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n",d->h1->target_name[tid], pos + 1 , rb, bam_nt16_rev_table[x>>28], qPosteriorSum, x>>8&0xff, ref_q, x>>16&0xff, n1, n2);
            }
        }    
        r = sniper_maqindel(n1, pos, d->ido, pl1, d->ref, 0,0, d->h1);
        if (r) {
            sniper_maqindel_ret_t *q=0;
            glf3_t *tumor_g3=0, *normal_g3=0;
            tumor_g3=sniper_maqindel2glf(r, n1);
            q=sniper_somaticindelscore(n2, pos, d->ido, pl2, d->ref,r);
            q->libs1=0;
            q->libs2=0;
            normal_g3=sniper_maqindel2glf(q, n2);
            int prior1 = tumor_g3->indel_len[0] ? 43 : 0;
            int prior2 = tumor_g3->indel_len[1] ? 43 : 0;
            int prior3 = prior1 && prior2 ? 80 : 40;
            int hetindeltumor = prob_indel(tumor_g3, prior3, HET_INDEL);
            int hom2indeltumor = prob_indel(tumor_g3, prior2, HOMO_INDEL2);
            int hom1indeltumor = prob_indel(tumor_g3, prior1, HOMO_INDEL1);
            int hetindelnormal = prob_indel(normal_g3, prior3, HET_INDEL);
            int hom2indelnormal = prob_indel(normal_g3, prior2, HOMO_INDEL2);
            int hom1indelnormal = prob_indel(normal_g3, prior1, HOMO_INDEL1);
            qPosteriorSum = 255;
            int het_sum = (hetindeltumor+hetindelnormal) > 255 ? 255 : hetindeltumor+hetindelnormal; 
            int hom1_sum = (hom1indeltumor+hom1indelnormal) > 255 ? 255 : hom1indeltumor+hom1indelnormal;  
            int hom2_sum = (hom2indeltumor+hom2indelnormal) > 255 ? 255 : hom2indeltumor+hom2indelnormal;  
            //fprintf(stdout,"het_sum:%d\thom1_sum%d\t:hom2_sum%d\n", het_sum, hom1_sum, hom2_sum);
            qPosteriorSum = qAdd(qPosteriorSum, het_sum);
            qPosteriorSum = qAdd(qPosteriorSum, hom1_sum);
            qPosteriorSum = qAdd(qPosteriorSum, hom2_sum);
            //if(minimum_somatic_qual >= qPosteriorSum) {  
            //print general info on the site in common among both tumor and normal
            fprintf(indel_fh, "%s\t%d\t*\t%d\t%s\t%s\t%d\t%d\t",d->h1->target_name[tid], pos + 1,qPosteriorSum, tumor_g3->indel_len[0] ? tumor_g3->indel_seq[0] : "*", tumor_g3->indel_len[1] ? tumor_g3->indel_seq[1] : "*", tumor_g3->indel_len[0], tumor_g3->indel_len[1]);
            //next output tumor specific info
            if (r->gt < 2) fprintf(indel_fh,"%s/%s\t", r->s[r->gt], r->s[r->gt]);
            else fprintf(indel_fh,"%s/%s\t", r->s[0], r->s[1]);
            fprintf(indel_fh,"%d\t%d\t", r->q_cns, r->q_ref);
            fprintf(indel_fh,"%d\t%d\t", gTumor->max_mapQ, n1);
            fprintf(indel_fh,"%d\t%d\t%d\t%d\t", r->cnt1, r->cnt2, r->cnt_ambi, r->cnt_anti);
            fprintf(indel_fh,"%d\t%d\t%d\t%d\t", tumor_g3->min_lk, tumor_g3->lk[0], tumor_g3->lk[1], tumor_g3->lk[2]);
            if (q->gt < 2) fprintf(indel_fh,"%s/%s\t", q->s[q->gt], q->s[q->gt]);
            else fprintf(indel_fh,"%s/%s\t", q->s[0], q->s[1]);
            fprintf(indel_fh,"%d\t%d\t", q->q_cns, q->q_ref);
            fprintf(indel_fh,"%d\t%d\t", gNormal->max_mapQ, n2);
            fprintf(indel_fh,"%d\t%d\t%d\t%d\t", q->cnt1, q->cnt2, q->cnt_ambi, q->cnt_anti);
            fprintf(indel_fh,"%d\t%d\t%d\t%d\t", normal_g3->min_lk, normal_g3->lk[0], normal_g3->lk[1], normal_g3->lk[2]);
            fprintf(indel_fh,"%d\t%d\n", r->libs1, r->libs2);
            //}
            glf3_destroy1(tumor_g3);
            glf3_destroy1(normal_g3);

            sniper_maqindel_ret_destroy(r);
            sniper_maqindel_ret_destroy(q);

        }
        free(gTumor);
        free(gNormal);
        return qPosteriorSum;
    }
    return -1;
}




int main(int argc, char *argv[])
{
    int c;
    char *fn_fa = 0;
    pu_data2_t *d = (pu_data2_t*)calloc(1, sizeof(pu_data2_t));
    d->min_somatic_qual=0;
    d->tid = -1; d->mask = BAM_DEF_MASK; d->mapQ = 0;
    d->c = sniper_maqcns_init();
    d->ido = sniper_maqindel_opt_init();
    while ((c = getopt(argc, argv, "f:T:N:r:I:G:q:Q:")) >= 0) {
        switch (c) {
            case 'f': fn_fa = strdup(optarg); break;
            case 'T': d->c->theta = atof(optarg); break;
            case 'N': d->c->n_hap = atoi(optarg); break;
            case 'r': d->c->het_rate = atoi(optarg); break;
            case 'q': d->mapQ = atoi(optarg); break;
            case 'I': d->ido->q_indel = atoi(optarg); break;
            case 'G': d->ido->r_indel = atof(optarg); break;
            case 'Q': d->min_somatic_qual = atoi(optarg); break;         
            default: fprintf(stderr, "Unrecognizd option '-%c'.\n", c); return 1;
        }
    }
    if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "somaticsniper [options] <tumor.bam> <normal.bam> <snp_output_file> <indel_output_file>\n\n");
        fprintf(stderr, "Options: \n");
        fprintf(stderr, "        -q INT    filtering reads with mapping quality less than INT [%d]\n", d->mapQ);
        fprintf(stderr, "        -Q INT    filtering somatic snp output(NOT INDELS!) with somatic quality less than  INT [15]\n");
        fprintf(stderr, "        -f FILE   reference sequence in the FASTA format\n\n");
        fprintf(stderr, "        -T FLOAT  theta in maq consensus calling model (for -c/-g) [%f]\n", d->c->theta);
        fprintf(stderr, "        -N INT    number of haplotypes in the sample (for -c/-g) [%d]\n", d->c->n_hap);

        fprintf(stderr, "        -r FLOAT  prior of a difference between two haplotypes (for -c/-g) [%f]\n", d->c->het_rate);
        fprintf(stderr, "        -G FLOAT  prior of an indel between two haplotypes (for -c/-g) [%f]\n", d->ido->r_indel);
        fprintf(stderr, "        -I INT    phred prob. of an indel in sequencing/prep. (for -c/-g) [%d]\n", d->ido->q_indel);
        fprintf(stderr, "\n");
        free(fn_fa); free(d);
        return 1;
    }
    if (fn_fa) d->fai = fai_load(fn_fa);
    free(fn_fa);
    sniper_maqcns_prepare(d->c);
    fprintf(stderr,"Preparing to snipe some somatics\n");
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
    FILE* indel_fh = fopen(argv[optind+3], "w");

    if(snp_fh && indel_fh) {
        //NEED TO ADD IN AN ACTUAL FUNCTION NAME HERE
        bam_sspileup_file(fp1, fp2, d->mask, d->mapQ, glf_somatic, d, snp_fh, indel_fh);
    }
    else {
        fprintf(stderr, "Unable to open snp or indel files!!!!!!!!!\n");
        exit(1);
    }

    bam_close(fp1);
    bam_close(fp2);
    kh_destroy(64, d->hash);
    bam_header_destroy(d->h1);
    bam_header_destroy(d->h2);
    if (d->fai) fai_destroy(d->fai);
    sniper_maqcns_destroy(d->c);
    free(d->ido); free(d->ref); free(d);
    return 0;
}


int prob_indel(glf3_t* indel, int prior_prob, int likelihood_flag) {
    int i; 
    int qSum=255;  
    for (i=0; i < 3; i++) {
        int lk_prior= (int)indel->lk[i] + prior_prob;
        qSum = qAdd(lk_prior, qSum);
    }
    return (indel->lk[likelihood_flag] + prior_prob - qSum);
}






