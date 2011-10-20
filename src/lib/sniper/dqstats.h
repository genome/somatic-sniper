#pragma once

/* per-site depth and quality statistics */

#include <bam.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* depth and quality statistics */
typedef struct {
    uint32_t mean_baseQ[4];
    uint32_t mean_mapQ[4];
    uint32_t base_occ[4];
    /* note about dp4 indices:
     * 0: reads supporting reference on forward strand
     * 1: reads supporting reference on reverse strand
     * 2: reads NOT supporting reference on forward strand
     * 3: reads NOT supporting reference on reverse strand
     */
    uint32_t dp4[4]; 
    uint32_t total_depth;
    uint32_t total_mean_mapQ;
} dqstats_t;

void get_dqstats(
        const bam_pileup1_t* buf,
        int n_reads,
        int ref_base,
        uint32_t wanted_bases,
        dqstats_t *dqstats
        );

void print_mean_quality_values(FILE* fh, int bases, const uint32_t values[4]);

void print_base_count(FILE* fh, int bases, const uint32_t counts[4]);

void print_dp4(FILE* fh, const uint32_t dp4[4]);

#ifdef __cplusplus
}
#endif
