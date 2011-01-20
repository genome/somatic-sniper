#ifndef SOMATIC_SNIPER
#define SOMATIC_SNIPER

#define expPhred(x) (double)exp((double)(-(x))/4.343)
#define logPhred(x) (int)((x) < 1 ? (0.5-4.343*log(x)) : (-0.5-4.343*log(x)))

typedef int (*bam_sspileup_f)(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE* snp_fh, FILE* indel_fh);
int get_next_pos(bam_plbuf_t *buf,bamFile fp); 
int bam_sspileup_file(bamFile fp1, bamFile fp2, int mask, int thresh, bam_sspileup_f func, void *func_data, FILE *snp_fh, FILE* indel_fh);

#endif

