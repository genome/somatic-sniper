#include <math.h>
#include <string.h>
#include "bam.h"
#include "sniper_maqcns.h"
#include "ksort.h"
KSORT_INIT_GENERIC(uint32_t)

typedef struct __sbmc_aux_t {
	int max;
	uint32_t *info;
} sniper_bmc_aux_t;

typedef struct {
	float esum[4], fsum[4];
	uint32_t c[4];
	uint32_t rms_mapQ;
} sniper_glf_call_aux_t;

char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };


/*
  P(<b1,b2>) = \theta \sum_{i=1}^{N-1} 1/i
  P(D|<b1,b2>) = \sum_{k=1}^{N-1} p_k 1/2 [(k/N)^n_2(1-k/N)^n_1 + (k/N)^n1(1-k/N)^n_2]
  p_k = i/k / \sum_{i=1}^{N-1} 1/i
 */
static void sniper_cal_het(sniper_maqcns_t *aa)
{
	int k, n1, n2;
	double sum_harmo; // harmonic sum
	double poly_rate;
	double p1 = 0.0, p3 = 0.0; // just for testing

	free(aa->lhet);
	aa->lhet = (double*)calloc(256 * 256, sizeof(double));
	sum_harmo = 0.0;
	for (k = 1; k <= aa->n_hap - 1; ++k)
		sum_harmo += 1.0 / k;
	for (n1 = 0; n1 < 256; ++n1) {
		for (n2 = 0; n2 < 256; ++n2) {
			long double sum = 0.0;
			double lC = lgamma(n1+n2+1) - lgamma(n1+1) - lgamma(n2+1); // \binom{n1+n2}{n1}
			for (k = 1; k <= aa->n_hap - 1; ++k) {
				double pk = 1.0 / k / sum_harmo;
				double log1 = log((double)k/aa->n_hap);
				double log2 = log(1.0 - (double)k/aa->n_hap);
				sum += pk * 0.5 * (expl(log1*n2) * expl(log2*n1) + expl(log1*n1) * expl(log2*n2));
			}
			aa->lhet[n1<<8|n2] = lC + logl(sum);
			if (n1 == 17 && n2 == 3) p3 = lC + logl(expl(logl(0.5) * 20));
			if (n1 == 19 && n2 == 1) p1 = lC + logl(expl(logl(0.5) * 20));
		}
	}
	poly_rate = aa->het_rate * sum_harmo;
	aa->q_r = -4.343 * log(2.0 * poly_rate / (1.0 - poly_rate));
}

/** initialize the helper structure */
static void sniper_cal_coef(sniper_maqcns_t *aa)
{
	int k, n, q;
	long double sum_a[257], b[256], q_c[256], tmp[256], fk2[256];
	double *lC;

	lC = (double*)calloc(256 * 256, sizeof(double));
	// aa->lhet will be allocated and initialized 
	free(aa->fk); free(aa->coef);
	aa->fk = (double*)calloc(256, sizeof(double));
	aa->coef = (double*)calloc(256*256*64, sizeof(double));
	aa->fk[0] = fk2[0] = 1.0;
	for (n = 1; n != 256; ++n) {
		aa->fk[n] = pow(aa->theta, n) * (1.0 - aa->eta) + aa->eta;
		fk2[n] = aa->fk[n>>1]; // this is an approximation, assuming reads equally likely come from both strands
	}
	for (n = 1; n != 256; ++n)
		for (k = 1; k <= n; ++k)
			lC[n<<8|k] = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
	for (q = 1; q != 64; ++q) {
		double e = pow(10.0, -q/10.0);
		double le = log(e);
		double le1 = log(1.0-e);
		for (n = 1; n != 256; ++n) {
			double *coef = aa->coef + (q<<16|n<<8);
			sum_a[n+1] = 0.0;
			for (k = n; k >= 0; --k) { // a_k = \sum_{i=k}^n C^n_k \epsilon^k (1-\epsilon)^{n-k}
				sum_a[k] = sum_a[k+1] + expl(lC[n<<8|k] + k*le + (n-k)*le1);
				b[k] = sum_a[k+1] / sum_a[k];
				if (b[k] > 0.99) b[k] = 0.99;
			}
			for (k = 0; k != n; ++k) // log(\bar\beta_{nk}(\bar\epsilon)^{f_k})
				q_c[k] = -4.343 * fk2[k] * logl(b[k] / e);
			for (k = 1; k != n; ++k) q_c[k] += q_c[k-1]; // \prod_{i=0}^k c_i
			for (k = 0; k <= n; ++k) { // powl() in 64-bit mode seems broken on my Mac OS X 10.4.9
				tmp[k] = -4.343 * logl(1.0 - expl(fk2[k] * logl(b[k])));
				coef[k] = (k? q_c[k-1] : 0) + tmp[k]; // this is the final c_{nk}
			}
		}
	}
	free(lC);
}

sniper_maqcns_t *sniper_maqcns_init()
{
	sniper_maqcns_t *bm;
	bm = (sniper_maqcns_t*)calloc(1, sizeof(sniper_maqcns_t));
	bm->aux = (sniper_bmc_aux_t*)calloc(1, sizeof(sniper_bmc_aux_t));
	bm->het_rate = 0.001;
	bm->theta = 0.85;
	bm->n_hap = 2;
	bm->eta = 0.03;
	bm->cap_mapQ = 60;
	return bm;
}

void sniper_maqcns_prepare(sniper_maqcns_t *bm)
{
	sniper_cal_coef(bm); sniper_cal_het(bm);
}

void sniper_maqcns_destroy(sniper_maqcns_t *bm)
{
	if (bm == 0) return;
	free(bm->lhet); free(bm->fk); free(bm->coef); free(bm->aux->info);
	free(bm->aux); free(bm);
}

glf1_t *sniper_maqcns_glfgen(int _n, const bam_pileup1_t *pl, uint8_t ref_base, sniper_maqcns_t *bm)
{
	sniper_glf_call_aux_t *b;
	int i, j, k, w[8], c, n;
	glf1_t *g = (glf1_t*)calloc(1, sizeof(glf1_t));
	float p[16], min_p = 1e30;
	uint64_t rms;

	g->ref_base = ref_base;
	if (_n == 0) return g;

	// construct aux array
	if (bm->aux->max < _n) {
		bm->aux->max = _n;
		kroundup32(bm->aux->max);
		bm->aux->info = (uint32_t*)realloc(bm->aux->info, 4 * bm->aux->max);
	}
	for (i = n = 0; i < _n; ++i) {
		const bam_pileup1_t *p = pl + i;
		uint32_t q, x = 0, qq;
		if (p->is_del || (p->b->core.flag&BAM_FUNMAP)) continue;
		q = (uint32_t)bam1_qual(p->b)[p->qpos];
		x |= (uint32_t)bam1_strand(p->b) << 18 | q << 8 | p->b->core.qual;
		if (p->b->core.qual < q) q = p->b->core.qual;
		x |= q << 24;
		qq = bam1_seqi(bam1_seq(p->b), p->qpos);
		q = bam_nt16_nt4_table[qq? qq : ref_base];
		if (!p->is_del && q < 4) x |= 1 << 21 | q << 16;
		bm->aux->info[n++] = x;
	}
	ks_introsort(uint32_t, n, bm->aux->info);
	// generate esum and fsum
	b = (sniper_glf_call_aux_t*)calloc(1, sizeof(sniper_glf_call_aux_t));
	for (k = 0; k != 8; ++k) w[k] = 0;
	rms = 0;
	for (j = n - 1; j >= 0; --j) { // calculate esum and fsum
		uint32_t info = bm->aux->info[j];
		int tmp;
		if (info>>24 < 4 && (info>>8&0x3f) != 0) info = 4<<24 | (info&0xffffff);
		k = info>>16&7;
		if (info>>24 > 0) {
			b->esum[k&3] += bm->fk[w[k]] * (info>>24);
			b->fsum[k&3] += bm->fk[w[k]];
			if (w[k] < 0xff) ++w[k];
			++b->c[k&3];
		}
		tmp = (int)(info&0x7f) < bm->cap_mapQ? (int)(info&0x7f) : bm->cap_mapQ;
		rms += tmp * tmp;
	}
	b->rms_mapQ = (uint8_t)(sqrt((double)rms / n) + .499);
	// rescale ->c[]
	for (j = c = 0; j != 4; ++j) c += b->c[j];
	if (c > 255) {
		for (j = 0; j != 4; ++j) b->c[j] = (int)(254.0 * b->c[j] / c + 0.5);
		for (j = c = 0; j != 4; ++j) c += b->c[j];
	}
	// generate likelihood
	for (j = 0; j != 4; ++j) {
		// homozygous
		float tmp1, tmp3;
		int tmp2, bar_e;
		for (k = 0, tmp1 = tmp3 = 0.0, tmp2 = 0; k != 4; ++k) {
			if (j == k) continue;
			tmp1 += b->esum[k]; tmp2 += b->c[k]; tmp3 += b->fsum[k];
		}
		if (tmp2) {
			bar_e = (int)(tmp1 / tmp3 + 0.5);
			if (bar_e < 4) bar_e = 4; // should not happen
			if (bar_e > 63) bar_e = 63;
			p[j<<2|j] = tmp1 + bm->coef[bar_e<<16|c<<8|tmp2];
		} else p[j<<2|j] = 0.0; // all the bases are j
		// heterozygous
		for (k = j + 1; k < 4; ++k) {
			for (i = 0, tmp2 = 0, tmp1 = tmp3 = 0.0; i != 4; ++i) {
				if (i == j || i == k) continue;
				tmp1 += b->esum[i]; tmp2 += b->c[i]; tmp3 += b->fsum[i];
			}
			if (tmp2) {
				bar_e = (int)(tmp1 / tmp3 + 0.5);
				if (bar_e < 4) bar_e = 4;
				if (bar_e > 63) bar_e = 63;
				p[j<<2|k] = p[k<<2|j] = -4.343 * bm->lhet[b->c[j]<<8|b->c[k]] + tmp1 + bm->coef[bar_e<<16|c<<8|tmp2];
			} else p[j<<2|k] = p[k<<2|j] = -4.343 * bm->lhet[b->c[j]<<8|b->c[k]]; // all the bases are either j or k
		}
		//
		for (k = 0; k != 4; ++k)
			if (p[j<<2|k] < 0.0) p[j<<2|k] = 0.0;
	}

	{ // fix p[k<<2|k]
		float max1, max2, min1, min2;
		int max_k, min_k;
		max_k = min_k = -1;
		max1 = max2 = -1.0; min1 = min2 = 1e30;
		for (k = 0; k < 4; ++k) {
			if (b->esum[k] > max1) {
				max2 = max1; max1 = b->esum[k]; max_k = k;
			} else if (b->esum[k] > max2) max2 = b->esum[k];
		}
		for (k = 0; k < 4; ++k) {
			if (p[k<<2|k] < min1) {
				min2 = min1; min1 = p[k<<2|k]; min_k = k;
			} else if (p[k<<2|k] < min2) min2 = p[k<<2|k];
		}
		if (max1 > max2 && (min_k != max_k || min1 + 1.0 > min2))
			p[max_k<<2|max_k] = min1 > 1.0? min1 - 1.0 : 0.0;
	}

	// convert necessary information to glf1_t
	g->ref_base = ref_base; g->max_mapQ = b->rms_mapQ;
	g->depth = n > 16777215? 16777215 : n;
	for (j = 0; j != 4; ++j)
		for (k = j; k < 4; ++k)
			if (p[j<<2|k] < min_p) min_p = p[j<<2|k];
	g->min_lk = min_p > 255.0? 255 : (int)(min_p + 0.5);
	for (j = c = 0; j != 4; ++j)
		for (k = j; k < 4; ++k)
			g->lk[c++] = p[j<<2|k]-min_p > 255.0? 255 : (int)(p[j<<2|k]-min_p + 0.5);

	free(b);
	return g;
}

uint32_t sniper_glf2cns(const glf1_t *g, int q_r)
{
	int i, j, k, tmp[16], min = 10000, min2 = 10000, min3 = 10000, min_g = -1, min_g2 = -1;
	uint32_t x = 0;
	for (i = k = 0; i < 4; ++i)
		for (j = i; j < 4; ++j) {
			tmp[j<<2|i] = -1;
			tmp[i<<2|j] = g->lk[k++] + (i == j? 0 : q_r);
		}
	for (i = 0; i < 16; ++i) {
		if (tmp[i] < 0) continue;
		if (tmp[i] < min) {
			min3 = min2; min2 = min; min = tmp[i]; min_g2 = min_g; min_g = i;
		} else if (tmp[i] < min2) {
			min3 = min2; min2 = tmp[i]; min_g2 = i;
		} else if (tmp[i] < min3) min3 = tmp[i];
	}
	x = min_g >= 0? (1U<<(min_g>>2&3) | 1U<<(min_g&3)) << 28 : 0xf << 28;
	x |= min_g2 >= 0? (1U<<(min_g2>>2&3) | 1U<<(min_g2&3)) << 24 : 0xf << 24;
	x |= (uint32_t)g->max_mapQ << 16;
	x |= min2 < 10000? (min2 - min < 256? min2 - min : 255) << 8 : 0xff << 8;
	x |= min2 < 10000 && min3 < 10000? (min3 - min2 < 256? min3 - min2 : 255) : 0xff;
	return x;
}

uint32_t sniper_maqcns_call(int n, const glf1_t *g, sniper_maqcns_t *bm)
{
	uint32_t x;
	if (n) {
		x = sniper_glf2cns(g, (int)(bm->q_r + 0.5));
	} else x = 0xfU<<28 | 0xfU<<24;
	return x;
}

