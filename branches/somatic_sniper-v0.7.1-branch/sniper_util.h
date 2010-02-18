#ifndef SNIPER_UTIL_H
#define SNIPER_UTIL_H

//Utility arrays for converting between base encodings
/* glf genotype order is: AA/AC/AG/AT/CC/CG/CT/GG/GT/TT, or AMRWCSYGKT in IUPAC */
static int glfBase[10] = { 1, 3, 5, 9, 2, 6, 10, 4, 12, 8 } ; /* mapping from 10 genotypes to 4 bit base coding */
static int baseGlf[16] = { -1, 0, 4, 1, 7, 2, 5, -1, 9, 3, 6, -1, 8, -1, -1, -1 } ; /* inverse of glfBase */

//Macros related to bam
#define BAM_REF_BASE       0  //The 4bit encoding of a base matching the reference
#define BAM_N_BASE         15 //The 4bit encoding of the N base               

#define isHom(x) ((x) != BAM_REF_BASE && !((x) & ((x) - 1))) //Test to see if the 4bit encoded nucleotide is Homozygous (a power of two). Got this off the web.
#define isHet(y) ((y) != BAM_REF_BASE && (y) != BAM_N_BASE && ((y) & ((y) - 1)))    //Test to see if a 4bit encoded nucleotide is Heterozygous.

#define nt4_to_nt2(x) ( isHom(x) ? log2(x): (-1) )

/* Functions for converting to and from Phred space */
#define expPhred(x) (double)exp((double)(-(x))/4.343)
#define logPhred(x) (int)((x) < 1 ? (0.5-4.343*log(x)) : (-0.5-4.343*log(x)))

//some math functions
int logAdd(int a, int c);

#endif
