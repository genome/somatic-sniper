#ifndef __mean_qualities_h__
#define __mean_qualities_h__

#include <bam.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void mean_quality_values(
        const bam_pileup1_t* buf,
        int n_reads,
        uint32_t wanted_bases,
        uint32_t mean_baseQ[4],
        uint32_t count_baseQ[4],
        uint32_t mean_mapQ[4],
        uint32_t count_mapQ[4]
        );

void print_mean_quality_values(FILE* fh, int bases, uint32_t values[4]);

void print_base_count(FILE* fh, int bases, uint32_t counts[4]);

#ifdef __cplusplus
}
#endif

#endif /* __mean_qualities_h__ */
