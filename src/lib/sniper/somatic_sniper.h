#ifndef SOMATIC_SNIPER
#define SOMATIC_SNIPER

#include "sniper_maqcns.h"

#include "sam.h"
#include "kstring.h"
#include "sam.h"
#include "faidx.h"

#define PHRED_CONST (4.343)
#define expPhred(x) (double)exp((double)(-(x))/PHRED_CONST)
#define logPhred(x) (int)((x) < 1 ? (0.5-PHRED_CONST*log(x)) : (-0.5-PHRED_CONST*log(x)))

#define FORMAT_SIMPLE 1 //simple file format output for backwards compatibility 
#define FORMAT_EXTENDED 2 //extended file format output

#ifdef __cplusplus
extern "C" {
#endif

/* types */
typedef struct {
    bam_header_t *h1;
    bam_header_t *h2;
    sniper_maqcns_t *c;
    faidx_t *fai;
    uint32_t format;
    int tid, len, last_pos;
    int mask;
    int mapQ;
    int min_somatic_qual;//for limiting snp calls in somatic sniper
    char *ref;
    glfFile fp; // for glf output only
    //for alternative prior construction
    int use_joint_priors;
    double somatic_mutation_rate;
} pu_data2_t;

typedef int (*bam_sspileup_f)(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE* snp_fh);

/* data */
extern int isHom[16];
extern int isHet[16];
extern int glfBase[10];


/* functions */
void qAddTableInit (void);
void makeSoloPrior (void);
void make_joint_prior (double somatic_rate);
int get_next_pos(bam_plbuf_t *buf,bamFile fp); 
int bam_sspileup_file(bamFile fp1, bamFile fp2, int mask, int thresh, bam_sspileup_f func, void *func_data, FILE *snp_fh);
int glf_somatic(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE *snp_fh);

#ifdef __cplusplus
}
#endif

#endif

