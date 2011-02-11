#ifndef SOMATIC_SNIPER
#define SOMATIC_SNIPER

#include "sniper_maqcns.h"

#include "sam.h"
#include "kstring.h"
#include "sam.h"
#include "faidx.h"

#define expPhred(x) (double)exp((double)(-(x))/4.343)
#define logPhred(x) (int)((x) < 1 ? (0.5-4.343*log(x)) : (-0.5-4.343*log(x)))

#ifdef __cplusplus
extern "C" {
#endif

/* types */
typedef struct {
    bam_header_t *h1;
    bam_header_t *h2;
    sniper_maqcns_t *c;
    sniper_maqindel_opt_t *ido;
    faidx_t *fai;
    uint32_t format;
    int tid, len, last_pos;
    int mask;
    int mapQ;
    int min_somatic_qual;//for limiting snp calls in somatic sniper
    char *ref;
    glfFile fp; // for glf output only
} pu_data2_t;

typedef int (*bam_sspileup_f)(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE* snp_fh, FILE* indel_fh);

/* data */
extern int isHom[16];
extern int isHet[16];
extern int glfBase[10];


/* functions */
void qAddTableInit (void);
void makeSoloPrior (void);
int get_next_pos(bam_plbuf_t *buf,bamFile fp); 
int bam_sspileup_file(bamFile fp1, bamFile fp2, int mask, int thresh, bam_sspileup_f func, void *func_data, FILE *snp_fh, FILE* indel_fh);
int glf_somatic(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE *snp_fh, FILE *indel_fh);

#ifdef __cplusplus
}
#endif

#endif

