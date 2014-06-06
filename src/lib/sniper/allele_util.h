#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/* The functions and macros in this module operate on alleles encoded as
 * integers where A=1, C=2, G=4, T=8, so A/C = 3, and so on.
 */

typedef enum {
    WILDTYPE = 0,
    GERMLINE = 1,
    SOMATIC  = 2,
    LOH      = 3,
    UNKNOWN  = 4
} variant_status_t;

/* what alleles do a and b have in common? */
#define genotype_intersection(a, b) ((a) & (b))

/* returns true if the alleles encoded in a are a proper subset of those in b */
#define genotype_is_proper_subset(a, b) \
    ((b) != (a) && genotype_intersection((a), (b)) == a)

/* returns the alleles encoded in a that do not exist in b (a \ b) */
#define genotype_set_difference(a, b) ((a) & ~(b))

/* alias for genotype_is_proper_subset in domain terms (loss of heterozygostity)
 * for convenience. LOH can be definied in different ways, so try to use
 * genotype_is_proper_subset to avoid ambiguity elsewhere in the code. */
#define is_loh(a, b) genotype_is_proper_subset((a), (b))


/* return the number of alleles encoded in the int 'a' */
int count_alleles(int a);

/* returns 1 if the site should be filtered out because it is LOH from either the tumor or normal */
int should_filter_as_loh(int ref_base, int tumor_genotype, int normal_genotype);

#ifdef __cplusplus
}
#endif /* __cplusplus */
