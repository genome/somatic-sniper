#include "allele_util.h"

/* return the number of alleles encoded in the int 'a',
 * where A=1, C=2, G=4, T=8
 */
int count_alleles(int a) {
    return (a & 1) + ((a>>1)&1) + ((a>>2)&1) + ((a>>3)&1);
}

/* returns 1 if the alleles in a represent a loss of heterozygosity event with respect to b */
int is_loh(int a, int b) {
    if ((a&b) != a) /* a contains something b doesn't, can't be LOH */
        return 0;

    return ((a&b) != b);
}
