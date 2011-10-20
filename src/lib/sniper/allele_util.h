#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef enum {
    WILDTYPE = 0,
    GERMLINE = 1,
    SOMATIC  = 2,
    LOH      = 3,
    UNKNOWN  = 4
} variant_status_t;

/* return the number of alleles encoded in the int 'a',
 * where A=1, C=2, G=4, T=8
 */
int count_alleles(int a);

/* returns 1 if the alleles in a represent a loss of heterozygosity event with respect to b */
int is_loh(int a, int b);

#ifdef __cplusplus
}
#endif /* __cplusplus */
