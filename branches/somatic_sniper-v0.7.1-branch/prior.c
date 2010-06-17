#include "prior.h"
#include "sniper_util.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam.h"

//prior probabilities of transversions and transitions
const static double somaticHomotransversionProb = 0.1111;
const static double somaticHomotransitionProb = 0.4444;


//prior probabilities of transversions and transitions
const static double somaticHettransversionProb = 0.07275;
const static double somaticHettransitionProb = 4.0*0.0725;

//prior probabilities of transversions and transitions
const static double germtransversionProb = 1.0/6.0;
const static double germtransitionProb = 4.0*1.0/6.0;

//prior probabilities of germline sites and somatic sites
static double germline_priors[16][10] ;  /* index over reference base, genotype. Stores precalculated prior probabilities for germline assumption */
static double diploid_transition_transversion[10][10]; //genotype transition/transversion probabilities.
static double prior_of_normal_genotype_yielding_somatic[10];
static double somatic_mutation_priors[3];  //store precalculated prior probability of somatic events

//this calculates whether the 2bit single base encodings are transitions or transversions and returns the appropriate probabilities
double transition_transversion_calc(int a, int b, float transitionProb, float transversionProb) {
    if(abs(a-b) == 2) {
        return transitionProb;
    }
    else {
        return transversionProb;
    }
}

//returns first bit set as its own number
int first_allele(int genotype) {
    if(genotype == 0) {
        return 0;
    }
    int mask = 1;
    while(!(genotype & mask)) {
        mask <<= 1;
    }
    return mask;
}

//calculate the prior of the genotype. 
//Basically should be a little less than 1 for germline
//the somatic mutation rate for single changes
//and the somatic mutation rate squared for double changes
double somatic_prior_for_genotype(int tumor_genotype, int normal_genotype) {
    if(tumor_genotype == normal_genotype) {
        //0 somatic changes
        return somatic_mutation_priors[0];
    }
    if(tumor_genotype & normal_genotype) {
        //they share an allele and therefore the number of changes is 1
        return somatic_mutation_priors[1];
    }
    else {
        //they don't share an allele and there must have been 2 events
        return somatic_mutation_priors[2];
    }
}

void initialize_somatic_priors(double prior) {
    somatic_mutation_priors[1] = logPhred(prior);
    somatic_mutation_priors[2] = somatic_mutation_priors[1];//*2;
    somatic_mutation_priors[0] = logPhred(1 - prior - prior*prior); //hopefully no underflow. I really need to start checking that.
}

double prior_for_genotype(int tumor_genotype, int normal_genotype, int ref) {

    if(tumor_genotype == normal_genotype) {
        return germline_priors[ref][tumor_genotype];
    }
    else {
        //factor in the predisposition to somatic mutation here
        return (diploid_transition_transversion[normal_genotype][tumor_genotype]+prior_of_normal_genotype_yielding_somatic[normal_genotype]);
    }

}

//This generates the germline priors based on the constant THETA. 
//TODO add more comments on what that actually means
void initialize_germline_priors (float prior) {
    int i, b, ref ;

    for (ref = 0 ; ref < 16 ; ++ref) {
        for (i = 0 ; i < 10 ; ++i) { 
            b = glfBase[i];
            int r = glfBase[ref];
            if(isHom(ref)) {

                if (!(b & ~ref))	/* ie b is compatible with ref */
                    germline_priors[ref][i] = 0;
                else if (b & ref) {	/* ie one allele of b is compatible with ref */
                    int a = b^r;    //changed allele
                    germline_priors[ref][i] = logPhred(prior) + logPhred(transition_transversion_calc(r,a,germtransitionProb,germtransversionProb));
                }
                else if (isHom(b))	/* single mutation homozygote */
                    germline_priors[ref][i] = logPhred(0.5*prior) + logPhred(transition_transversion_calc(r,b,germtransitionProb,germtransversionProb));
                else {			/* two mutations */
                    int a = first_allele(b);
                    germline_priors[ref][i] = logPhred(prior*prior) + logPhred(transition_transversion_calc(r,a,germtransitionProb,germtransversionProb)) + logPhred(transition_transversion_calc(r,a^b,germtransitionProb,germtransversionProb));
                }
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
    fprintf(stderr,"Initializing default somatic probabilities based on transition/transversion probabilities\n");
    for(i = 0; i < 10; ++i) {
        //for now all normal genotypes are equally likely so
        prior_of_normal_genotype_yielding_somatic[i] = logPhred(0.1);   
        
        int total_probability = 10000;  //scaling variable to make sure these end up adding up to one
        double transitionProb = 0.0;
        double transversionProb = 0.0;
        if(isHom(glfBase[i])) {
            transversionProb = somaticHomotransversionProb;
            transitionProb = somaticHomotransitionProb;
        }
        else {
            transversionProb = somaticHettransversionProb;
            transitionProb = somaticHettransitionProb;
        }
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
                        double prob = transition_transversion_calc(a_alleles[0],b_alleles[0],transitionProb,transversionProb);
                        diploid_transition_transversion[i][j] = 2 * logPhred(prob);
                    }
                    else {
                        //a is hom and b is het, the probability is the prob of different alleles * each other
                        int b_allele;
                        double prob = 0.0;
                        for(b_allele = 0; b_allele < 2; b_allele++) {
                            if(b_alleles[b_allele] != a_alleles[0]) {
                                prob += logPhred(transition_transversion_calc(a_alleles[0],b_alleles[b_allele],transitionProb,transversionProb));
                            }
                        }
                        diploid_transition_transversion[i][j] = prob;
                    }

                }
                else {
                    if(num_alleles_in_b == 1) {
                        //is b is hom and a is het, the prob is the same as for same situation above 
                        int a_allele;
                        double prob = 0.0;
                        for(a_allele = 0; a_allele < 2; a_allele++) {
                            if(a_alleles[a_allele] != b_alleles[0]) {
                                prob += logPhred(transition_transversion_calc(b_alleles[0],a_alleles[a_allele],transitionProb,transversionProb));
                            }
                        }
                        diploid_transition_transversion[i][j] = prob;

                    }
                    else {
                        //both are het
                        //if they share an allele then it is the change of the variant alleles
                        //do this bitwise
                        if(a & b) {
                            if(a_alleles[0] == b_alleles[0]) {
                                diploid_transition_transversion[i][j] = logPhred(transition_transversion_calc(a_alleles[1],b_alleles[1],transitionProb,transversionProb));
                            }
                            else if(a_alleles[0] == b_alleles[1]) {
                                diploid_transition_transversion[i][j] = logPhred(transition_transversion_calc(a_alleles[1],b_alleles[0],transitionProb,transversionProb));
                            }
                            else if(a_alleles[1] == b_alleles[0]) {
                                diploid_transition_transversion[i][j] = logPhred(transition_transversion_calc(a_alleles[0],b_alleles[1],transitionProb,transversionProb));
                            }
                            else {
                                diploid_transition_transversion[i][j] = logPhred(transition_transversion_calc(a_alleles[0],b_alleles[0],transitionProb,transversionProb));
                            }
                            
                        }
                        else {
                        //if they share no allele then it is the square of the transition
                            diploid_transition_transversion[i][j] = logPhred( transition_transversion_calc(3,1,transitionProb,transversionProb))*2;
                        }
                    }
                }
            }
            total_probability = logAdd(total_probability,diploid_transition_transversion[i][j]);
        }
        //normalize priors
        //for(j = 0; j < 10; ++j) {
        //    diploid_transition_transversion[i][j] -= total_probability;
        //}
            
    }
}

/*
//initializes the somatic priors with a "neutral" mutation model
void initialize_diploid_transition_transversion() {
    int i,j;
    fprintf(stderr,"Initializing default somatic probabilities based on transition/transversion probabilities\n");
    for(i = 0; i < 10; ++i) {
        //for now all normal genotypes are equally likely so
        prior_of_normal_genotype_yielding_somatic[i] = logPhred(0.1);   
        
        int total_probability = 10000;  //scaling variable to make sure these end up adding up to one
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
                        diploid_transition_transversion[i][j] = 2 * logPhred(prob);
                    }
                    else {
                        //a is hom and b is het, the probability is the prob of different alleles * each other
                        int b_allele;
                        int prob = 1.0;
                        for(b_allele = 0; b_allele < 2; b_allele++) {
                            if(b_alleles[b_allele] != a_alleles[0]) {
                                prob += logPhred(transition_transversion_calc(a_alleles[0],b_alleles[b_allele]));
                            }
                        }
                        diploid_transition_transversion[i][j] = prob;
                    }

                }
                else {
                    //if(num_alleles_in_b == 1) {
                    //    //is b is hom and a is het, the prob is the same as for same situation above 
                    //    int a_allele;
                    //    double prob = 1.0;
                    //    for(a_allele = 0; a_allele < 2; a_allele++) {
                    //        if(a_alleles[a_allele] != b_alleles[0]) {
                    //            prob *= transition_transversion_calc(b_alleles[0],a_alleles[a_allele]);
                    //        }
                    //    }
                    //    diploid_transition_transversion[i][j] = logPhred(prob);

                    //}
                    //else {
                        //if a is het or both are het then it is the product of the combination of those alleles divided by two (number of combinations)
                        diploid_transition_transversion[i][j] = logPhred( transition_transversion_calc(a_alleles[0],b_alleles[0]) * transition_transversion_calc(a_alleles[1],b_alleles[1]) + transition_transversion_calc(a_alleles[0],b_alleles[1]) * transition_transversion_calc(a_alleles[1],b_alleles[0])) - logPhred(2);    
                    //}
                }
                total_probability = logAdd(total_probability,diploid_transition_transversion[i][j]);
            }
        }
        //normalize priors
        for(j = 0; j < 10; ++j) {
            diploid_transition_transversion[i][j] -= total_probability;
        }
            
    }
}*/

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
        double sum = 0.0;
        for(j = 0; j < 10; ++j) {
           fprintf(stderr,"\t%f",expPhred(diploid_transition_transversion[i][j]));
           //fprintf(stderr,"\t%f",expPhred(diploid_transition_transversion[i][j]));
           if(i!=j)
           //sum = logAdd(sum,diploid_transition_transversion[i][j]);
           sum += expPhred(diploid_transition_transversion[i][j]);
        }
        fprintf(stderr,"\tSum: %f\n", sum);
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
        double sum = 1000.0;
        fprintf(stderr,"%c",a);
        for(j = 0; j < 10; ++j) {
           fprintf(stderr,"\t%f",germline_priors[i][j]);
           sum = logAdd(sum, germline_priors[i][j]);
        }
        fprintf(stderr,"\tSum: %f\n", sum);
    }
}

//loads the prior probabilities from a matrix file
int load_priors_from_file(FILE *file, double prior) {
        //check header
        char header[20];
        char *expected_header = "\tA\tM\tR\tW\tC\tS\tY\tG\tK\tT";
        int result = fscanf(file,"%20c ",header);
        header[20]='\0';    //null terminate the string
        if(result == 0 || result == EOF || strcmp(expected_header,header) != 0) {
            fprintf(stderr, "Improperly formatted header line for the somatic probabilities\n");
            return 0;
        }
        //load somatic probabilities
        int i=0;
        
        for(i = 0; i<10; i++) {
            result = 0;
            char base;
            double raw_probabilities[10];
            result = fscanf(file,"%c\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf ",&base,&raw_probabilities[0],&raw_probabilities[1],&raw_probabilities[2],&raw_probabilities[3],&raw_probabilities[4],&raw_probabilities[5],&raw_probabilities[6],&raw_probabilities[7],&raw_probabilities[8],&raw_probabilities[9]);
            if(result != 11) {
                fprintf(stderr, "Improperly formatted probability line for the somatic probabilities\n");
                return 0;
            }
            int j = 0;
            for(j=0; j < 10; j++) {
                diploid_transition_transversion[i][j] = logPhred(raw_probabilities[j]);
            }
        } 
        //initalize germline priors'
        initialize_germline_priors(prior);
        return 1;
}
