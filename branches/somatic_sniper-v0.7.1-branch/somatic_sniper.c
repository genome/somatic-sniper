#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include "sam.h"
#include "faidx.h"
#include "sniper_maqcns.h"
#include "khash.h"
#include "kstring.h"
#include "somatic_sniper.h"
#include "prior.h"
#include "sniper_util.h"

//TODO not sure what these are for yet
typedef int *indel_list_t;
KHASH_MAP_INIT_INT64(64, indel_list_t)

//These duplicate what was in the BAM code and are largely useless at this point I think    
#define BAM_PLF_SIMPLE     0x01
#define BAM_PLF_CNS        0x02
#define BAM_PLF_INDEL_ONLY 0x04
#define BAM_PLF_GLF        0x08
#define BAM_PLF_VAR_ONLY   0x10

//Macros to improve clarity    
#define HOMO_INDEL1 0
#define HOMO_INDEL2 1
#define HET_INDEL   2
    

static int qAddTable[1024]; //This is a predefined table of adds in Phred space. It is filled upon execution


//do phred equivalent of log(x + y) where both x and y are already in log space
#define qAdd(x,y)  (x - qAddTable[512+y-x])




//a big giant ass heng li type struct for storing variables that need to get passed to the callback
typedef struct {
    bam_header_t *tumor_header;   //probably tumor. This is a terrible name
    bam_header_t *normal_header;   //probably normal TODO this is a terrible name too
    sniper_maqcns_t *c; //options etc used by the snp calling model
    sniper_maqindel_opt_t *ido; //options etc used by indel caller
    faidx_t *fai;               //reference sequence index
    khash_t(64) *hash;          //TODO figure out wtf this is for
    uint32_t format;            //largely or entirely irrelevant
    int tid, len, last_pos;     //info about where in the reference the program is at
    int mask;                   //mask to determine what reads to allow
    int mapQ;                   //minimum mapping quality value TODO Check this is true
    int min_somatic_qual;   //for limiting snp calls in somatic sniper
    char *ref;              //the reference seqeunce in characters for the current chromosome
    
    //somatic specific
    float somatic_rate; //rate of somatic point mutations

    //new model for storing reads from both in relatively static memory
    bam_pileup1_t *pl3;  //for storing both tumor and normal alignments
    int max_n3; //for knowing how many tumor and normal alignments we have space for

    
} pu_data2_t;

//TODO determine what this is actually doing 
int prob_indel(glf3_t* indel, int prior_prob, int likelihood_flag);




sniper_maqindel_ret_t *sniper_somaticindelscore(int n, int pos, const sniper_maqindel_opt_t *mi, const bam_pileup1_t *pl, const char *ref, sniper_maqindel_ret_t *tumor_indel);

//I assume these are actually from libbam but weren't declared in a header.
char **__bam_get_lines(const char *fn, int *_n);
void bam_init_header_hash(bam_header_t *header);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);


void print_glf(glf1_t *g) {
    fprintf(stderr, "\tmin_lk: %d\n",g->min_lk);
    int i = 0;

    for(i=0; i<10; i++) {
        fprintf(stderr, "\t%c: %d\n",bam_nt16_rev_table[glfBase[i]],g->lk[i]);
    }

}


//Initialize the precalculated phred space add table
//TODO try to explain wtf this is doing better
static void qAddTableInit (void)
{
    int i ;
    for (i = 0 ; i < 1000 ; ++i)
    { double e = 1 + expPhred(i-512) ;
        qAddTable[i] = -logPhred(e) ;
    }
}

glf3_t *sniper_maqindel2glf(sniper_maqindel_ret_t *r, int n) {
    glf3_t *g3;
    g3 = glf3_init1();

    int het = 3 * n, min;
    min = het;
    if (min > r->gl[0]) min = r->gl[0];
    if (min > r->gl[1]) min = r->gl[1];
    g3->ref_base = 0;
    g3->rtype = GLF3_RTYPE_INDEL;
    memset(g3->lk, 0, 10);
    g3->lk[0] = r->gl[0] - min < 255? r->gl[0] - min : 255;
    g3->lk[1] = r->gl[1] - min < 255? r->gl[1] - min : 255;
    g3->lk[2] = het - min < 255? het - min : 255;
    g3->offset = 0;
    g3->indel_len[0] = r->indel1;
    g3->indel_len[1] = r->indel2;
    g3->min_lk = min < 255? min : 255;
    g3->max_len = (abs(r->indel1) > abs(r->indel2)? abs(r->indel1) : abs(r->indel2)) + 1;
    g3->indel_seq[0] = strdup(r->s[0]+1);
    g3->indel_seq[1] = strdup(r->s[1]+1);
    return g3;
}



static int glf_somatic(uint32_t tid, uint32_t pos, int n1, int n2, const bam_pileup1_t *pl1, const bam_pileup1_t *pl2, void *data, FILE *snp_fh, FILE *indel_fh) {
    //hacked copy from function gl3_func behavior to get a g with 10 probabilities to do somatic probability calculation    
    pu_data2_t *d = (pu_data2_t*)data;

    //fix up accumulated alignments from tumor AND normal
    /*if((n1 + n2) > d->max_n3) {
        //allocate space to store the reads
        d->pl3 = realloc(d->pl3, sizeof(bam_pileup1_t) * (n1 + n2));
        if(!d->pl3) {
            fprintf(stderr, "Couldn't allocate memory for tumor/normal combined reads\n");
            exit(1);
        }
        else {
            d->max_n3 = n1 + n2;
        }
    }

    //copy over tumor and normal reads
    memcpy(d->pl3,pl1,sizeof(bam_pileup1_t)*n1);
    memcpy(d->pl3 + n1, pl2, sizeof(bam_pileup1_t)*n2); //I hope to god this works
      */      

    sniper_maqindel_ret_t *r = 0;
    if (d->fai && (int)tid != d->tid) {
        free(d->ref);
        d->ref = fai_fetch(d->fai, d->tumor_header->target_name[tid], &d->len);
        d->tid = tid;
    }
    int rb = (d->ref && (int)pos < d->len)? d->ref[pos] : 'N';
    int rb4 = bam_nt16_table[rb];

    int lkSomatic[10][10];

    int qPosteriorSum = 255;
    int qSomatic = 255;
    int qProbabilityData = 255;

    glf1_t *gTumor =sniper_maqcns_glfgen(n1, pl1, bam_nt16_table[rb], d->c);
    //fprintf(stderr, "Tumor likelihood for %d:\n",pos+1);
    //print_glf(gTumor);
    glf1_t *gNormal =sniper_maqcns_glfgen(n2, pl2, bam_nt16_table[rb], d->c);
    //fprintf(stderr, "Normal likelihood for %d:\n",pos+1);
    //print_glf(gNormal);
    //glf1_t *gAll = sniper_maqcns_glfgen(n1 + n2, d->pl3, bam_nt16_table[rb], d->c);
    //fprintf(stderr, "All likelihood for %d:\n",pos+1);
    //print_glf(gAll);
        
    //now we have the filled g1,g2 to compare with code from larsons's glfSomatic
    
    if (rb != 'N' && gTumor->depth > 0 && gNormal->depth > 0) {
        //calculate posteriors
            int tumor, normal;
            int min_lk_tumor_genotype = -1, min_lk_normal_genotype = -1, min_joint_lk = 1000;
            for(tumor = 0; tumor < 10; tumor++) {
                for(normal = 0; normal < 10; normal++) {
                    lkSomatic[tumor][normal] = (gTumor->lk[tumor] + gTumor->min_lk) + (gNormal->lk[normal] + gNormal->min_lk) + prior_for_genotype(tumor,normal,rb4) + somatic_prior_for_genotype(tumor,normal);
                    /*                    if(lkSomatic[tumor][normal] > 255) {
                        lkSomatic[tumor][normal] = 255;
                    }
                    */
                    qProbabilityData = logAdd(lkSomatic[tumor][normal],qProbabilityData );
                    if(lkSomatic[tumor][normal] < min_joint_lk) {
                        min_joint_lk = lkSomatic[tumor][normal];
                        min_lk_normal_genotype = normal;
                        min_lk_tumor_genotype = tumor;
                    }
                }
            }
            for(tumor = 0; tumor < 10; tumor++) {
                for(normal = 0; normal < 10; normal++) {
                    lkSomatic[tumor][normal] -= qProbabilityData;
                    //cap low likelihoods
                    //if(lkSomatic[tumor][normal] > 255) {
                    //    lkSomatic[tumor][normal] = 255;
                    //}   
                    if(tumor != min_lk_tumor_genotype && normal != min_lk_normal_genotype)
                        qPosteriorSum = logAdd(lkSomatic[tumor][normal],qPosteriorSum);
                    if(tumor == normal) {
                        qSomatic = logAdd(lkSomatic[tumor][normal], qSomatic);
                    }
                }
            }

            // int result = logAdd(0,-qPosteriorSum);
            if(d->min_somatic_qual <= qSomatic) {  
                fprintf(snp_fh, "%s\t%d\t%c\t%c\t%c\t%d\t%d\t%d\t%d\n",d->tumor_header->target_name[tid], pos + 1 , rb, bam_nt16_rev_table[glfBase[min_lk_tumor_genotype]], bam_nt16_rev_table[glfBase[min_lk_normal_genotype]], qSomatic, qPosteriorSum, n1, n2);
                //original caller compatible format
                //fprintf(snp_fh, "%s\t%d\t%c\t%c\t%d\t%d\t%d\t%d\n",d->tumor_header->target_name[tid], pos + 1 , rb, bam_nt16_rev_table[glfBase[min_lk_tumor_genotype]], qSomatic, qPosteriorSum, n1, n2);
            }
            
        r = sniper_maqindel(n1, pos, d->ido, pl1, d->ref, 0,0, d->tumor_header);
        if (r) {
            sniper_maqindel_ret_t *q=0;
            glf3_t *tumor_g3=0, *normal_g3=0;
            tumor_g3=sniper_maqindel2glf(r, n1);
            q=sniper_somaticindelscore(n2, pos, d->ido, pl2, d->ref,r);
            q->libs1=0;
            q->libs2=0;
            normal_g3=sniper_maqindel2glf(q, n2);
            int prior1 = tumor_g3->indel_len[0] ? 43 : 0;
            int prior2 = tumor_g3->indel_len[1] ? 43 : 0;
            int prior3 = prior1 && prior2 ? 80 : 40;
            int hetindeltumor = prob_indel(tumor_g3, prior3, HET_INDEL);
            int hom2indeltumor = prob_indel(tumor_g3, prior2, HOMO_INDEL2);
            int hom1indeltumor = prob_indel(tumor_g3, prior1, HOMO_INDEL1);
            int hetindelnormal = prob_indel(normal_g3, prior3, HET_INDEL);
            int hom2indelnormal = prob_indel(normal_g3, prior2, HOMO_INDEL2);
            int hom1indelnormal = prob_indel(normal_g3, prior1, HOMO_INDEL1);
            qPosteriorSum = 255;
            int het_sum = (hetindeltumor+hetindelnormal) > 255 ? 255 : hetindeltumor+hetindelnormal; 
            int hom1_sum = (hom1indeltumor+hom1indelnormal) > 255 ? 255 : hom1indeltumor+hom1indelnormal;  
            int hom2_sum = (hom2indeltumor+hom2indelnormal) > 255 ? 255 : hom2indeltumor+hom2indelnormal;  
            //fprintf(stdout,"het_sum:%d\thom1_sum%d\t:hom2_sum%d\n", het_sum, hom1_sum, hom2_sum);
            qPosteriorSum = logAdd(qPosteriorSum, het_sum);
            qPosteriorSum = logAdd(qPosteriorSum, hom1_sum);
            qPosteriorSum = logAdd(qPosteriorSum, hom2_sum);
            //print general info on the site in common among both tumor and normal
            fprintf(indel_fh, "%s\t%d\t*\t%d\t%s\t%s\t%d\t%d\t",d->tumor_header->target_name[tid], pos + 1,qPosteriorSum, tumor_g3->indel_len[0] ? tumor_g3->indel_seq[0] : "*", tumor_g3->indel_len[1] ? tumor_g3->indel_seq[1] : "*", tumor_g3->indel_len[0], tumor_g3->indel_len[1]);
            //next output tumor specific info
            if (r->gt < 2) fprintf(indel_fh,"%s/%s\t", r->s[r->gt], r->s[r->gt]);
            else fprintf(indel_fh,"%s/%s\t", r->s[0], r->s[1]);
            fprintf(indel_fh,"%d\t%d\t", r->q_cns, r->q_ref);
            fprintf(indel_fh,"%d\t%d\t", gTumor->max_mapQ, n1);
            fprintf(indel_fh,"%d\t%d\t%d\t%d\t", r->cnt1, r->cnt2, r->cnt_ambi, r->cnt_anti);
            fprintf(indel_fh,"%d\t%d\t%d\t%d\t", tumor_g3->min_lk, tumor_g3->lk[0], tumor_g3->lk[1], tumor_g3->lk[2]);
            if (q->gt < 2) fprintf(indel_fh,"%s/%s\t", q->s[q->gt], q->s[q->gt]);
            else fprintf(indel_fh,"%s/%s\t", q->s[0], q->s[1]);
            fprintf(indel_fh,"%d\t%d\t", q->q_cns, q->q_ref);
            fprintf(indel_fh,"%d\t%d\t", gNormal->max_mapQ, n2);
            fprintf(indel_fh,"%d\t%d\t%d\t%d\t", q->cnt1, q->cnt2, q->cnt_ambi, q->cnt_anti);
            fprintf(indel_fh,"%d\t%d\t%d\t%d\t", normal_g3->min_lk, normal_g3->lk[0], normal_g3->lk[1], normal_g3->lk[2]);
            fprintf(indel_fh,"%d\t%d\n", r->libs1, r->libs2);
            glf3_destroy1(tumor_g3);
            glf3_destroy1(normal_g3);

            sniper_maqindel_ret_destroy(r);
            sniper_maqindel_ret_destroy(q);

        }
        free(gTumor);
        free(gNormal);
        return qPosteriorSum;
    }
    return -1;
}




int main(int argc, char *argv[])
{
    int c;  //For processing command line options
    char *fn_fa = NULL;    //pointer to the reference sequence filename
    char *prior_prob = NULL; //pointer to the prior somatic probabilities filename
    
    pu_data2_t *d = (pu_data2_t*)calloc(1, sizeof(pu_data2_t)); //This stores file info and variables related to the calculating of likelihoods 
    
    //initialize some basic filter variables
    d->min_somatic_qual = 0;  //by default there is no minimum somatic quality
    d->tid = -1; //stores the current chromosome
    d->mask = BAM_DEF_MASK; 
    d->mapQ = 0;
    d->somatic_rate = 0.00001; //default to 1e-6. TODO test some other possible values
    d->c = sniper_maqcns_init();
    d->ido = sniper_maqindel_opt_init();
    d->pl3 = NULL;
    d->max_n3 = 0;
    
    while ((c = getopt(argc, argv, "f:T:N:r:I:G:q:Q:s:p:")) >= 0) {
        switch (c) {
            case 'f': fn_fa = strdup(optarg); break;    //reference file name.
            case 'T': d->c->theta = atof(optarg); break;    //error correlation correction
            case 'N': d->c->n_hap = atoi(optarg); break;    //number of haplotypes
            case 'r': d->c->het_rate = atoi(optarg); break; //het rate of the population
            case 'q': d->mapQ = atoi(optarg); break;        //minimum mapping quality
            case 'I': d->ido->q_indel = atoi(optarg); break;    //probability of indel in prep
            case 'G': d->ido->r_indel = atof(optarg); break;    //probability of a indel between haplotypes
            case 'Q': d->min_somatic_qual = atoi(optarg); break;//minimum somatic quality to report of SNVs         
            case 's': d->somatic_rate = atof(optarg); break;   //probability of observing a somatic point mutuation         
            case 'p': prior_prob = strdup(optarg); break;   //somatic genotype probability file         
            default: fprintf(stderr, "Unrecognizd option '-%c'.\n", c); return 1;
        }
    }
    if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "somaticsniper [options] -f <ref.fasta> <tumor.bam> <normal.bam> <snp_output_file> <indel_output_file>\n\n");
        fprintf(stderr, "Required Option: \n");
        fprintf(stderr, "        -f FILE   REQUIRED reference sequence in the FASTA format\n\n");
        fprintf(stderr, "Options: \n");
        fprintf(stderr, "        -q INT    filtering reads with mapping quality less than INT [%d]\n", d->mapQ);
        fprintf(stderr, "        -Q INT    filtering somatic snp output(NOT INDELS!) with somatic quality less than  INT [15]\n");
        fprintf(stderr, "        -T FLOAT  theta in maq consensus calling model [%f]\n", d->c->theta);
        fprintf(stderr, "        -N INT    number of haplotypes in the sample [%d]\n", d->c->n_hap);
        fprintf(stderr, "        -r FLOAT  prior of a difference between two haplotypes in the population [%f]\n", d->c->het_rate);
        fprintf(stderr, "        -G FLOAT  prior of an indel between two haplotypes [%f]\n", d->ido->r_indel);
        fprintf(stderr, "        -I INT    phred prob. of an indel in sequencing/prep [%d]\n", d->ido->q_indel);
        fprintf(stderr, "        -s FLOAT  prior of a somatic point mutation occuring [%f]\n", d->somatic_rate);
        fprintf(stderr, "        -p FILE   somatic genotype probability file\n");
        fprintf(stderr, "\n");
        free(fn_fa); sniper_maqcns_destroy(d->c); free(d->ido); free(d);
        if(prior_prob) {
            free(prior_prob);
        }
        return 1;
    }

    //Handle the reading of the reference sequence and alert the user that it is a required option. TODO discuss fixing this with others.
    if (fn_fa) {
        d->fai = fai_load(fn_fa);
    }
    else {
        fprintf(stderr, "You MUST specify a reference sequence. It isn't optional.\n");
        sniper_maqcns_destroy(d->c);
        free(d->ido); 
        free(d);
        return 1; 
    }
    free(fn_fa);

    //Prepare to run the caller
    sniper_maqcns_prepare(d->c); //This precalculates some tables.
    qAddTableInit();    //Initialize a precalculated table of phred space additions ie when x and y are log space and you need log(x + y)
    if(prior_prob) {
        FILE *somatic_prior_prob = fopen(prior_prob,"r");
        if(somatic_prior_prob == NULL) {
            //TODO check errno here
            fprintf(stderr, "Error opening somatic probability file %s.\n",prior_prob); 
            //TODO report a more coherent error message taking into account the errno value
            free(prior_prob);
            sniper_maqcns_destroy(d->c);
            free(d->ido); 
            free(d);
            return 1;
        }
        else {
            free(prior_prob);
            if(load_priors_from_file(somatic_prior_prob,d->c->het_rate) == 0) {
                //there was some sort of error
                sniper_maqcns_destroy(d->c);
                free(d->ido); 
                free(d);
                fclose(somatic_prior_prob);
                return 1;
            }
        }

    }
    else {
        initialize_diploid_transition_transversion();
        initialize_germline_priors(d->c->het_rate);
        initialize_somatic_priors(d->somatic_rate);
    }
        
    print_transition_tranversion_priors();
    fprintf(stderr,"\n\n");

    print_germline_priors();
    fprintf(stderr,"\n\n");
    
    fprintf(stderr,"Preparing to snipe some somatics\n");
    bamFile tumor_fp, normal_fp;
    tumor_fp = bam_open(argv[optind], "r"); 
    if(tumor_fp) {
        fprintf(stderr, "Tumor bam is %s\n", argv[optind]);
        d->tumor_header = bam_header_read(tumor_fp);    //FIXME check an error code here
        sam_header_parse_rg(d->tumor_header);           //FIXME check an error code here
    }
    else {
        fprintf(stderr, "Unable to open %s as tumor bam file\n",argv[optind]);
        if (d->fai) fai_destroy(d->fai);
        sniper_maqcns_destroy(d->c);
        free(d->ido); free(d->ref); free(d);
        return 1;    
    }

    //open the normal bam file and parse the header
    normal_fp = bam_open(argv[optind+1], "r");
    if(normal_fp) {
        fprintf(stderr, "Normal bam is %s\n", argv[optind+1]);
        d->normal_header = bam_header_read(normal_fp);  //FIXME check and error here
        sam_header_parse_rg(d->normal_header);          //FIXME check for an error here
    }
    else {
        fprintf(stderr, "Unable to open %s as tumor bam file\n",argv[optind]);
        bam_close(tumor_fp);
        bam_header_destroy(d->tumor_header);
        if (d->fai) fai_destroy(d->fai);
        sniper_maqcns_destroy(d->c);
        free(d->ido); free(d->ref); free(d);
        return 1;    
    }
    


    FILE* snp_fh = fopen(argv[optind+2], "w");
    FILE* indel_fh = fopen(argv[optind+3], "w");

    if(snp_fh && indel_fh) {
        //NEED TO ADD IN AN ACTUAL FUNCTION NAME HERE
        bam_sspileup_file(tumor_fp, normal_fp, d->mask, d->mapQ, glf_somatic, d, snp_fh, indel_fh);
    }
    else {
        fprintf(stderr, "Unable to open snp or indel files!!!!!!!!!\n");
        exit(1);
    }

    bam_close(tumor_fp);
    bam_close(normal_fp);
    kh_destroy(64, d->hash);
    bam_header_destroy(d->tumor_header);
    bam_header_destroy(d->normal_header);
    if (d->fai) fai_destroy(d->fai);
    sniper_maqcns_destroy(d->c);
    free(d->ido); free(d->ref); free(d);
    fclose(snp_fh);
    fclose(indel_fh);
    return 0;
}


int prob_indel(glf3_t* indel, int prior_prob, int likelihood_flag) {
    int i; 
    int qSum=255;  
    for (i=0; i < 3; i++) {
        int lk_prior= (int)indel->lk[i] + prior_prob;
        qSum = logAdd(lk_prior, qSum);
    }
    return (indel->lk[likelihood_flag] + prior_prob - qSum);
}






