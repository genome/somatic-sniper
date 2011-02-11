#ifndef SNIPER_MAQCNS_H
#define SNIPER_MAQCNS_H

#include "sam.h"
#include "glf.h"

#ifdef __cplusplus
extern "C" {
#endif

struct __sbmc_aux_t;

typedef struct {
	float het_rate, theta;
	int n_hap, cap_mapQ;

	float eta, q_r;
	double *fk, *coef;
	double *lhet;
	struct __sbmc_aux_t *aux;
} sniper_maqcns_t;

typedef struct {
	int q_indel;
	float r_indel;
	// hidden parameters, unchangeable from command line
	int mm_penalty, indel_err, ambi_thres;
} sniper_maqindel_opt_t;

typedef struct {
	int indel1, indel2;
	int cnt1, cnt2, cnt_ambi, cnt_anti;
	char *s[2];
	//
	int gt, gl[2];
	int q_cns, q_ref;
    int libs1, libs2;
} sniper_maqindel_ret_t;

	sniper_maqcns_t *sniper_maqcns_init();
	void sniper_maqcns_prepare(sniper_maqcns_t *bm);
	void sniper_maqcns_destroy(sniper_maqcns_t *bm);
	glf1_t *sniper_maqcns_glfgen(int n, const bam_pileup1_t *pl, uint8_t ref_base, sniper_maqcns_t *bm);
	uint32_t sniper_maqcns_call(int n, const glf1_t *g, sniper_maqcns_t *bm);
	// return: cns<<28 | cns2<<24 | mapQ<<16 | cnsQ<<8 | cnsQ2
	uint32_t sniper_glf2cns(const glf1_t *g, int q_r);

	sniper_maqindel_opt_t *sniper_maqindel_opt_init();
	sniper_maqindel_ret_t *sniper_maqindel(int n, int pos, const sniper_maqindel_opt_t *mi, const bam_pileup1_t *pl, const char *ref,
									 int _n_types, int *_types, void * data);
	void sniper_maqindel_ret_destroy(sniper_maqindel_ret_t*);

#ifdef __cplusplus
}
#endif

#endif
