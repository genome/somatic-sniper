#include "allele_util.h"

/* return the number of alleles encoded in the int 'a',
 * where A=1, C=2, G=4, T=8
 */
int count_alleles(int a) {
    return (a & 1) + ((a>>1)&1) + ((a>>2)&1) + ((a>>3)&1);
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
    return /* normal LOH in tumor */
        genotype_is_proper_subset(tumor_genotype, normal_genotype)
        || /* or gain of reference allele in the tumor */
        (genotype_is_proper_subset(normal_genotype, tumor_genotype)
            && genotype_set_difference(tumor_genotype, normal_genotype) == ref_base);
}
