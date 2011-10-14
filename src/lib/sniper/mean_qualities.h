#ifndef __mean_qualities_h__
#define __mean_qualities_h__

#include <bam.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* each member is an array, one element per strand */
typedef struct {
    unsigned supporting_ref[2]; 
    unsigned variants[2]; 
} dp4_t;

void mean_quality_values(
        const bam_pileup1_t* buf,
        int n_reads,
        int ref_base,
        uint32_t wanted_bases,
        uint32_t mean_baseQ[4],
        uint32_t mean_mapQ[4],
        uint32_t base_occ[4],
        dp4_t* dp4
        );

void print_mean_quality_values(FILE* fh, int bases, uint32_t values[4]);

void print_base_count(FILE* fh, int bases, uint32_t counts[4]);

void print_dp4(FILE* fh, const dp4_t* dp4);

#ifdef __cplusplus
}
#endif

#endif /* __mean_qualities_h__ */
