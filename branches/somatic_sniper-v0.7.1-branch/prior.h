#ifndef PRIOR_H
#define PRIOR_H
#include <stdio.h>

//loads the prior probabilities from a matrix file
int load_priors_from_file(FILE *file, double prior);

//returns the appropriate prior probability for a given reference base, and tumor and normal genotypes
double prior_for_genotype(int tumor_genotype, int normal_genotype, int ref);

//This generates the germline priors based on the prior probability of a het mutation. 
//TODO add more comments on what the actually means
void initialize_germline_priors (float prior);

//This generates the somatic priors based on transition/transversion probabilities    
void initialize_diploid_transition_transversion();

//This generates the somatic priors based on transition/transversion probabilities    
void initialize_somatic_priors(double prior);

double somatic_prior_for_genotype(int tumor_genotype, int normal_genotype);


//Print out the somatic priors
void print_transition_tranversion_priors();
void print_germline_priors();

#endif
