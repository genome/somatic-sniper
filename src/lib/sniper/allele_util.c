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

/* returns 1 if the site should be filtered out because it is LOH from either the tumor or normal 
 * R T  N
 * A AA AG LOH
 * A GG AG LOH
 * A GG CG LOH
 * A CC CG LOH
 * A AG GG LOH not really but should be filtered
 * A AG AA Not LOH
 * */
int should_filter_as_loh(int ref_base, int tumor_genotype, int normal_genotype) {
    if (is_loh(tumor_genotype, normal_genotype) || (is_loh(normal_genotype, tumor_genotype) && normal_genotype != ref_base)) {
        return 1;
    }
    else {
        return 0;
    }
}
