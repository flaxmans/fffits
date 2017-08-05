

/*  ************* INCLUDES FOR SPECIAL DATA TYPES ***************** */
#include <gsl_rng.h>            // gnu scientific library
#include <gsl_randist.h>        // gnu scientific library


/*  ************* DEFAULTS FOR GLOBALS ***************** */
// NOTE: these are overridden by anything in parameters.ini.txt
#define nSITES_DEFAULT 1000000 // imaginary length of DNA we are watching evolve
#define nLINKAGE_GROUPS_DEFAULT 2 // # of independently assorting units
#define nGENERATIONS_DEFAULT 10
#define nPOPULATIONS_DEFAULT 3 // # islands/metapopulations/demes
#define MU_DEFAULT 1.0E-9 // per base mutation rate per meiosis per generation
#define MIGRATION_RATE_DEFAULT 0.01
#define K_DEFAULT 1000.0 // carrying capacity
#define INCLUDE_SELECTION_DEFAULT 1
#define MEAN_S_DEFAULT 0.01
#define FIXED_POP_SIZE_DEFAULT 0
#define MAX_POP_GROWTH_RATE_DEFAULT 0.01
#define RECOMBINATION_RATE_PER_KB_DEFAULT 0.001
#define PROBABILITY_SITE_SELECTED_DEFAULT 0.01
#define TIME_SERIES_SAMPLE_FREQ_DEFAULT 1000


/*  ************* SITE CLASSIFICATIONS AND MAGIC NUMBERS ***************** */
#define nSITE_CLASSES 4
#define SITE_CLASS_NEUTRAL 0
#define SITE_CLASS_BGS 1 // background selection
#define SITE_CLASS_POS 2 // positive selection
#define SITE_CLASS_DIV 3 // divergent selection
#define ALLELE_CODE_ANCESTRAL 0
#define ALLELE_CODE_DERIVED 1
#define PLOIDY 2 // diploid
#define ENVT_TYPE_GRADIENT 0 // this is default in fffits.c
#define ENVT_TYPE_MOSAIC 1
#define ENVT_TYPE_INVARIANT 2
#define ENVT_MAX_DEFAULT 1.0 // max environmental value -- like a crude "niche"
#define ENVT_MIN_DEFAULT -1.0 // min environmental value
#define CODOMINANCE 0.5 // codominant effects of alleles in diploid genotypes
#define FITNESS_MODEL_ADDITIVE 0
#define FITNESS_MODEL_MULTIPLICATIVE 1 // default in fffits.c
#define LOCUS_STATUS_INACTIVE 0
#define LOCUS_STATUS_VARIABLE_IN_PARENTS 1
#define LOCUS_STATUS_VARIABLE_PLUS_MUT 2
#define LOCUS_STATUS_NEW_MUT_ONLY 3
#define LOCUS_STATUS_TRACKED_IN_PARENTS 4


/*  ************* FUNCTION DECLARATIONS ***************** */
// initialization functions
void initializeRNG(unsigned int RNG_SEED);
int initializationSteps( int argc, char *argv[], char *progname );
double randExp(double meanValue);
unsigned readInParametersFromFile(void);
void setUpDataFiles(void);
void setUpGenome(void);
void setUpInitialAlleleFrequencies(double *expectedFreq, unsigned long long int *SFScounts);
void setUpPopulations(void);
void usage(char *progname);
void wrongParametersIniOption(char *expected, char *previous, char *found);
// pop gen calculations:
void calcAlleleCountsFreqsFSTandSFS(double *alleleFreqsByPop, double *globalFreqs);
void calculateAndPrintPi(double *alleleFreqsByPop, double *globalFreqs);
void calculateAndPrintDXY(double *alleleFreqsByPop);
void calculateAlleleCountsByPop(long int focalSite, long int *alleleCountsByPopulation);
void calculateAlleleFreqsByPop(long int *alleleCountsByPopulation, double *alleleFreqsByPop);
void calculateFST( double *FSTarray, double *alleleFreqsByPop, double *globalFreqs );
void dataRecording(void);
void writeAbundances(void);


/*  ************* GLOBAL VARIABLES ***************** */
extern _Bool VERBOSE;
extern long int nGENERATIONS, nSelectedSites;
extern const gsl_rng_type *rngType;		/* generator type */
extern gsl_rng *rngState;				/* rng instance */
extern int nPOPULATIONS, nLINKAGE_GROUPS;
extern long int N, *abundances;						// population size in total
extern unsigned long long int nSITES, blockSizes[2];	// number of sites in genome, memory size for genome data
extern double MU;						// per base mutation rate
extern unsigned long long int nTrackedSitesInParents;
extern short int *sitesStatuses, *genotypes0, *genotypes1, *gts, currentBlock;
extern int *siteClassifications, *linkageGroupMembership, *locations;
extern double *alleleFrequencies, *K_VALUES, *selectionCoefficients;
extern unsigned long long int *parentalTrackedSiteIndexes, *siteIndexes, *alleleCounts; // pointers for memory blocks for sites in genome
extern double PROBABILITY_SITE_DIV, PROBABILITY_SITE_POS, PROBABILITY_SITE_BGS, PROBABILITY_SITE_NEUTRAL;
extern double MEAN_S_BGS, MEAN_S_DIV, MEAN_S_POS;
extern int ENVIRONMENT_TYPE; // defaults for how selection works; magic numbers defined below
extern double *environmentGradient, ENVT_MAX, ENVT_MIN;
extern int nDEMOGRAPHIC_CHANGES;
extern long int *DEMOGRAPHIC_CHANGE_TIMES;
extern _Bool FIXED_POP_SIZE;
extern double MAX_POP_GROWTH_RATE;
extern int nMIGRATION_CHANGES;
extern long int *MIGRATION_CHANGE_TIMES;
extern double *M_VALUES;
extern double RECOMBINATION_RATE_PER_KB;
extern _Bool INCLUDE_SELECTION;
extern int FITNESS_MODEL;
extern long int TIME_SERIES_SAMPLE_FREQ;
extern double *migRatePt, *KvalPt;
extern double GENOME_MU;
extern _Bool TEST_MODE;
extern long int t;

// global file pointers
extern FILE *dataFile_alleleFreqTS;
extern FILE *dataFile_SFS_TS;
extern FILE *dataFile_segSiteTS;
extern FILE *dataFile_derivedFixationTS;
extern FILE *dataFile_abundances;
extern FILE *dataFile_PiAndDXY;





