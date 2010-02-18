#include "prior.h"
#include "sniper_util.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sam.h"

//prior probabilities of transversions and transitions
static double transitionProb = 0.66666667;
static double transversionProb = 0.16666667;

//prior probabilities of germline sites and somatic sites
static int germline_priors[16][10] ;  /* index over reference base, genotype. Stores precalculated prior probabilities for germline assumption */
static int diploid_transition_transversion[10][10]; //genotype transition/transversion probabilities.
static int prior_of_normal_genotype_yielding_somatic[10];

//this calculates whether the 2bit single base encodings are transitions or transversions and returns the appropriate probabilities
double transition_transversion_calc(int a, int b) {
    if(abs(a-b) == 2) {
        return transitionProb;
    }
    else {
        return transversionProb;
    }
}

int prior_for_genotype(int tumor_genotype, int normal_genotype, int ref) {

    if(tumor_genotype == normal_genotype) {
        return germline_priors[ref][tumor_genotype];
    }
    else {
        //factor in the predisposition to somatic mutation here
        return (diploid_transition_transversion[normal_genotype][tumor_genotype]+prior_of_normal_genotype_yielding_somatic[normal_genotype]);
    }

}

//This generates the germline priors based on the constant THETA. 
//TODO add more comments on what the actually means
void initialize_germline_priors (float prior) {
    int i, b, ref ;

    for (ref = 0 ; ref < 16 ; ++ref) {
        for (i = 0 ; i < 10 ; ++i) { 
            b = glfBase[i];
            int r = baseGlf[ref];
            if(r > -1) {

                if (!(b & ~ref))	/* ie b is compatible with ref */
                    germline_priors[ref][i] = 0;
                else if (b & ref)	/* ie one allele of b is compatible with ref */
                    germline_priors[ref][i] = logPhred(prior) + diploid_transition_transversion[r][i];
                else if (isHom(b))	/* single mutation homozygote */
                    germline_priors[ref][i] = logPhred(0.5*prior) + diploid_transition_transversion[r][i];
                else			/* two mutations */
                    germline_priors[ref][i] = logPhred(prior*prior) + diploid_transition_transversion[r][i];
            }
            else {
                germline_priors[ref][i] = 255;  //make non homozygous reference bases really unlikely...FIXME at some point when we care
            }
        }
    }
}

//initializes the somatic priors with a "neutral" mutation model
void initialize_diploid_transition_transversion() {
    int i,j;
    for(i = 0; i < 10; ++i) {
        //for now all normal genotypes are equally likely so
        prior_of_normal_genotype_yielding_somatic[i] = logPhred(0.1);   
        for(j = 0; j < 10; ++j) {
            //initialize based on the number of changes and phasing etc
            if(i==j) {
                diploid_transition_transversion[i][j] = 255;    //no change so no probability of change
            }
            else {
                char a = glfBase[i];
                char b = glfBase[j];
                
                //decompose into alleles, there must be a smarter way to do this
                int a_alleles[2], b_alleles[2], num_alleles_in_a = 0, num_alleles_in_b = 0;
                int bit = 0, a_allele_index = 0, b_allele_index = 0;//which bit we are dealing with
                //expand into an array of values indicating whether or not the bits are set
                for(bit=0; bit < 4; bit++) {
                    int a_bit_set = (a>>bit) & 1;
                    int b_bit_set = (b>>bit) & 1;
                    if(a_bit_set) {
                        a_alleles[a_allele_index] = bit;
                        a_allele_index++;
                        num_alleles_in_a++;
                    }
                    if(b_bit_set) {
                        b_alleles[b_allele_index] = bit;
                        b_allele_index++;
                        num_alleles_in_b++;
                    }
                }

                if(num_alleles_in_a == 1) {
                    if(num_alleles_in_b == 1) {
                        //both are homozygous, probability is the change squared
                        double prob = transition_transversion_calc(a_alleles[0],b_alleles[0]);
                        diploid_transition_transversion[i][j] = logPhred(prob * prob);
                    }
                    else {
                        //a is hom and b is het, the probability is the prob of different alleles * each other
                        int b_allele;
                        double prob = 1.0;
                        for(b_allele = 0; b_allele < 2; b_allele++) {
                            if(b_alleles[b_allele] != a_alleles[0]) {
                                prob *= transition_transversion_calc(a_alleles[0],b_alleles[b_allele]);
                            }
                        }
                        diploid_transition_transversion[i][j] = logPhred(prob);
                    }

                }
                else {
                    if(num_alleles_in_b == 1) {
                        //is b is hom and a is het, the prob is the same as for same situation above 
                        int a_allele;
                        double prob = 1.0;
                        for(a_allele = 0; a_allele < 2; a_allele++) {
                            if(a_alleles[a_allele] != b_alleles[0]) {
                                prob *= transition_transversion_calc(b_alleles[0],a_alleles[a_allele]);
                            }
                        }
                        diploid_transition_transversion[i][j] = logPhred(prob);

                    }
                    else {
                        //if both are het then it is the product of the combination of those alleles divided by two (number of combinations)
                        diploid_transition_transversion[i][j] = logPhred( transition_transversion_calc(a_alleles[0],b_alleles[0]) * transition_transversion_calc(a_alleles[1],b_alleles[1]) + transition_transversion_calc(a_alleles[0],b_alleles[1]) * transition_transversion_calc(a_alleles[1],b_alleles[0])) - logPhred(2);    
                    }
                }
            }
        }
    }
}

void print_transition_tranversion_priors() {
    int i,j;
    for(i = 0; i < 10; ++i) {
        char a = bam_nt16_rev_table[glfBase[i]];
        fprintf(stderr,"\t%c",a);
    }
    fprintf(stderr,"\n");
    for(i = 0; i < 10; ++i) {
        char a = bam_nt16_rev_table[glfBase[i]];
        fprintf(stderr,"%c",a);
        int sum = 255;
        for(j = 0; j < 10; ++j) {
           //fprintf(stderr,"\t%i",diploid_transition_transversion[i][j]);
           fprintf(stderr,"\t%f",expPhred(diploid_transition_transversion[i][j]));
           sum = logAdd(sum,diploid_transition_transversion[i][j]);
        }
        fprintf(stderr,"\tSum: %d\n", sum);
    }
}

void print_germline_priors() {
    int i,j;
    for(j = 0; j < 10; ++j) {
        char b = bam_nt16_rev_table[glfBase[j]];
        fprintf(stderr,"\t%c",b);
    }
    fprintf(stderr,"\n");
    for(i = 0; i < 16; ++i) {
        char a = bam_nt16_rev_table[i];
        int sum = 1000;
        fprintf(stderr,"%c",a);
        for(j = 0; j < 10; ++j) {
           fprintf(stderr,"\t%i",germline_priors[i][j]);
           sum = logAdd(sum, germline_priors[i][j]);
        }
        fprintf(stderr,"\tSum: %d\n", sum);
    }
}
