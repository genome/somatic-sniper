#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "dqstats.h"
#include "output_format.h"
#include "somatic_sniper.h"

int get_next_pos(bam_plbuf_t *buf,bamFile fp);

static int qAddTable[1024];
double THETA = 0.001 ;      /* population scaled mutation rate */
static int prior[16][10] ;  /* index over reference base, genotype */

#define qAdd(x,y)  (x + qAddTable[512+y-x])

char **__bam_get_lines(const char *fn, int *_n);
void bam_init_header_hash(bam_header_t *header);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);

int isHom[16] = {0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0} ;
int isHet[16] = {0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0} ;
int glfBase[10] = { 1, 3, 5, 9, 2, 6, 10, 4, 12, 8 } ; /* mapping from 10 genotypes to 4 bit base coding */

void makeSoloPrior (void) {
    int i, b, ref ;
    for (ref = 0 ; ref < 16 ; ++ref) {

        for (i = 0 ; i < 10 ; ++i) { 
            b = glfBase[i] ;
            if (!(b & ~ref))	/* ie b is compatible with ref */
                prior[ref][i] = 0 ;
            else if (b & ref)	/* ie one allele of b is compatible with ref */
                prior[ref][i] = logPhred(THETA) ;
            else if (isHom[b])	/* single mutation homozygote */
                prior[ref][i] = logPhred(0.5*THETA) ;
            else			/* two mutations */
                prior[ref][i] = logPhred(THETA*THETA) ;
        }
    }
}

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
    for (j = 0 ; j < 10 ; ++j) {
        lkResult[j] -= qSum;
        if(lkResult[j] > 255) {
            lkResult[j] = 255;
        }
    }

}

void qAddTableInit (void) {
    int i ;
    for (i = 0 ; i < 1000 ; ++i) { 
        double e = 1 + expPhred(i-512) ;
        qAddTable[i] = logPhred(e) ;
    }
}

int glf_somatic(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE *snp_fh) {
    //hacked copy from function gl3_func behavior to get a g with 10 probabilities to do somatic probability calculation
    pu_data2_t *d = (pu_data2_t*)data;
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
    
    //now we have the filled g1,g2 to compare with code from glfSomatic
    if (rb != 'N' && gTumor->depth > 0 && gNormal->depth > 0) {
        //calculate tumor posteriors

        uint32_t tumor_cns = sniper_maqcns_call(n1, gTumor, d->c);
        uint32_t normal_cns = sniper_maqcns_call(n2, gNormal, d->c);
        int rb4 = bam_nt16_table[rb];
        int tumor_base1 = tumor_cns >> 28;
        int tumor_base2 = tumor_cns >> 24 & 0xf;
        int tumor_score1 = tumor_cns >> 8 & 0xff;
        int tumor_score2 = tumor_cns & 0xff;
/*
        int tumor_rms_mapping = tumor_cns >> 16 & 0xff;
*/

        int normal_base1 = normal_cns >> 28;
        int normal_base2 = normal_cns >> 24 & 0xf;
        int normal_score1 = normal_cns >> 8 & 0xff;
        int normal_score2 = normal_cns & 0xff;
/*
        int normal_rms_mapping = normal_cns >> 16 & 0xff;
*/

        int tumor_snp_q = 0;
        int normal_snp_q = 0;

        if (rb4 != 15 && tumor_base1 != 15 && tumor_base1 != rb4) { // a SNP
            tumor_snp_q = (tumor_base2 == rb4)? tumor_score1 : tumor_score1 + tumor_score2;
            if (tumor_snp_q > 255) tumor_snp_q = 255;

            if (normal_base1 != 15 && normal_base1 != rb4)
            {
                normal_snp_q = (normal_base2 == rb4)? normal_score1 : normal_score1 + normal_score2;
                if (normal_snp_q > 255) normal_snp_q = 255;
            }

            //calculate the posterior probabilities
            calculatePosteriors(gTumor, lkTumor);
            calculatePosteriors(gNormal, lkNormal);

            int j;
            for(j = 0; j < 10; j++) {
                qPosteriorSum = qAdd(qPosteriorSum,(lkTumor[j] + lkNormal[j]));
            }

            if(d->min_somatic_qual <= qPosteriorSum) {
                sniper_output_t out;
                out.seq_name = d->h1->target_name[tid];
                out.pos = pos;
                out.ref_base = rb;
                out.ref_base4 = rb4;

                out.tumor.genotype = tumor_base1;
                out.tumor.consensus_quality = tumor_score1;
                out.tumor.variant_allele_quality = tumor_snp_q;
                out.tumor.somatic_score = qPosteriorSum;
                out.tumor.variant_status = qPosteriorSum > 0 ? SOMATIC : UNKNOWN;
                get_dqstats(pl1, n1, rb4, rb4|tumor_base1, &out.tumor.dqstats);

                out.normal.genotype = normal_base1;
                out.normal.consensus_quality = normal_score1;
                out.normal.variant_allele_quality = normal_snp_q;
                out.normal.somatic_score = -1;
                out.normal.variant_status = UNKNOWN;
                get_dqstats(pl2, n2, rb4, rb4|normal_base1, &out.normal.dqstats);

                d->output_formatter->output_fn(snp_fh, &out);
                fflush(snp_fh);
            }
        }
        free(gTumor);
        free(gNormal);
        return qPosteriorSum;
    }
    return -1;
}
