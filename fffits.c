/* fast forward-in-time simulator for prediction of population genetic metrics and ABC inference */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <string.h>

#include <gsl_rng.h>            // gnu scientific library //
#include <gsl_randist.h>        // gnu scientific library //
const gsl_rng_type *rngType;    /* generator type */
gsl_rng *rngState;              /* rng instance */



const char *version = "ffits1.0.0";

// defaults for globals
#define nSITES_DEFAULT 1000000 // imaginary length of DNA we are watching evolve
#define nLINKAGE_GROUPS_DEFAULT 2 // # of independently assorting units
#define nGENERATIONS_DEFAULT 10
#define nPOPULATIONS_DEFAULT 3 // # islands/metapopulations/demes
#define MU_DEFAULT 1E-9 // per base mutation rate per meiosis per generation
#define MIGRATION_RATE_DEFAULT 0.01
#define K_DEFAULT 1000.0
#define INCLUDE_SELECTION_DEFAULT 1
#define MEAN_S_DEFAULT 0.01
#define FIXED_POP_SIZE_DEFAULT 0
#define MAX_POP_GROWTH_RATE_DEFAULT 0.01
#define RECOMBINATION_RATE_PER_KB_DEFAULT 0.001
#define PROBABILITY_SITE_SELECTED_DEFAULT 0.01

// globals for cmd line and or parameter file options
unsigned long long int nSITES = nSITES_DEFAULT;
long int nGENERATIONS = nGENERATIONS_DEFAULT;
int nPOPULATIONS = nPOPULATIONS_DEFAULT;
double MU = MU_DEFAULT;
int nLINKAGE_GROUPS = nLINKAGE_GROUPS_DEFAULT;
int nDEMOGRAPHIC_CHANGES;
long int *DEMOGRAPHIC_CHANGE_TIMES, *MIGRATION_CHANGE_TIMES;
int nMIGRATION_CHANGES;
_Bool INCLUDE_SELECTION = INCLUDE_SELECTION_DEFAULT, FIXED_POP_SIZE = FIXED_POP_SIZE_DEFAULT;
double MEAN_S_BGS = MEAN_S_DEFAULT, MEAN_S_DIV = MEAN_S_DEFAULT, MEAN_S_POS = MEAN_S_DEFAULT;
double MAX_POP_GROWTH_RATE = MAX_POP_GROWTH_RATE_DEFAULT;
double RECOMBINATION_RATE_PER_KB = RECOMBINATION_RATE_PER_KB_DEFAULT;
double PROBABILITY_SITE_DIV = PROBABILITY_SITE_SELECTED_DEFAULT;
double PROBABILITY_SITE_POS = PROBABILITY_SITE_SELECTED_DEFAULT;
double PROBABILITY_SITE_BGS = PROBABILITY_SITE_SELECTED_DEFAULT;
double PROBABILITY_SITE_NEUTRAL;
int ENVIRONMENT_TYPE = 0, FITNESS_MODEL = 0; // defaults for how selection works; magic numbers defined below
double *environmentGradient, ENVT_MAX = 1.0, ENVT_MIN = -1.0;
long int nSelectedSites;

short int *genotypes0, *genotypes1, *gts; // pointer to memory blocks for individual genotypes
int *locations; // locations in discrete space of individuals
unsigned long long int *variableSiteIndexes, *siteIndexes, *SFScounts, *alleleCounts; // pointers for memory blocks for sites in genome
double *alleleFrequencies, *migrationRates, *K_VALUES, *M_VALUES, *selectionCoefficients;
_Bool *siteIsVariable;
short int currentBlock = 0; // keep track of which block is in use
unsigned long long int blockSizes[2];
long int t; // time counter
long int N; // current total population size
long int *abundances; // size of each population
unsigned long long int nVariableSites = 0;
_Bool VERBOSE = 0;
int currentMigrationPeriod = 0;
int currentDemographyPeriod = 0;
double *migRatePt, *KvalPt;

// site classifications and magic numbers
int *siteClassifications;
#define nSITE_TYPES 4
#define SITE_CLASS_NEUTRAL 0
#define SITE_CLASS_BGS 1 // background selection
#define SITE_CLASS_POS 2 // positive selection
#define SITE_CLASS_DIV 3 // divergent selection
#define ALLELE_CODE_ANCESTRAL 0
#define ALLELE_CODE_DERIVED 1
#define PLOIDY 2 // diploid
#define ENVT_TYPE_GRADIENT 0
#define ENVT_TYPE_MOSAIC 1
#define ENVT_TYPE_INVARIANT 2
#define CODOMINANCE 0.5 // codominant effects of alleles in diploid genotypes
#define FITNESS_MODEL_ADDITIVE 0
#define FITNESS_MODEL_MULTIPLICATIVE 1


// function declarations
long int calculateNumOffspring(int pop);
short int * checkMemoryBlocks(int totalOffspring);
void computeFitness(double *fitnessValues);
void finalTasks(void);
void initializeRNG(unsigned int RNG_SEED);
void migration(void);
void printParametersToFiles(unsigned RNG_SEED);
double randExp(double meanValue);
unsigned readInParametersFromFile(void);
void reproduction(void);
void setUpGenome(void);
void setUpInitialAlleleFrequencies(double *expectedFreq);
void setUpPopulations(void);
void usage(char *progname);
//void viabilitySelection(void);
void wrongParametersIniOption(char *expected, char *previous, char *found);



int main(int argc, char *argv[])
{
    unsigned int RNG_SEED;
    int ch;
    char *progname = argv[0];
    
    // read in optional command line arguments
    while ((ch = getopt(argc, argv, "V?")) != -1) {
        switch (ch) {
            case 'V':
                VERBOSE = 1;
                break;
            case '?':
            default:
                usage(progname);
                exit(-1);
        }
    }
    
    // initialization steps
    RNG_SEED = readInParametersFromFile();
    initializeRNG(RNG_SEED);
    setUpGenome();
    setUpPopulations();
    
    for ( t = 1; t <= nGENERATIONS; t++ ) {

        migration();
        
        reproduction(); // includes fecundity selection
        
    }
    
    finalTasks();
    return 0;
}


long int calculateNumOffspring(int pop)
{
    long int numOffspring = 0, focalTime;
    double expected, currentPop, k;
    
    // check to see if now is a time to change demographic parameters
    if ( nDEMOGRAPHIC_CHANGES > 0 ) {
        focalTime = DEMOGRAPHIC_CHANGE_TIMES[currentDemographyPeriod];
        if ( t == focalTime ) {
            KvalPt += nPOPULATIONS;
            currentDemographyPeriod++;
        }
    }
    
    if ( FIXED_POP_SIZE )
        numOffspring = (long int) *(KvalPt + pop);
    else {
        k = *(KvalPt + pop);
        currentPop = ((double) *(abundances + pop));
        expected = currentPop + (currentPop * MAX_POP_GROWTH_RATE * (k - currentPop)/k); // logistic equation
        numOffspring = gsl_ran_poisson( rngState, expected );
    }
    
    return numOffspring;
}



short int * checkMemoryBlocks(int totalOffspring)
{
    short int *offGT;
    int offspringBlock;
    unsigned long long int neededSize;
    
    if ( currentBlock )
        offspringBlock = 0; // current parents use 1, so offspring will use 0
    else
        offspringBlock = 1; // current 0, offspring 1
    
    neededSize = PLOIDY * totalOffspring * nVariableSites * (sizeof(short int));
    
    // check if we need bigger memory blocks and allocate if needed
    if ( blockSizes[offspringBlock] < neededSize ) {
        if ( offspringBlock ) {
            free( genotypes1 );
            genotypes1 = (short int *) malloc( neededSize );
            offGT = genotypes1;
        }
        else {
            free( genotypes0 );
            genotypes0 = (short int *) malloc( neededSize );
            offGT = genotypes0;
        }
        blockSizes[offspringBlock] = neededSize;
    }
    
    return offGT;
}


void computeFitness(double *fitnessValues)
{
    unsigned long long int i, j, selectedSiteIndexesMaster[nSelectedSites];
    unsigned long long int selectedSiteIndexesLocal[nSelectedSites], *ullpt;
    double cf_scoeffs[nSelectedSites], *dpt, cf_s, cf_gradient; // prefix cf_ denoting local to this function
    long int nFound;
    int cf_siteClass, cf_siteClasses[nSelectedSites], *locpt;
    short int *sipt, gtsum;
    
    dpt = fitnessValues;
    for ( i = 0; i < N; i++ ) {
        *dpt = 1.0;
        dpt++;
    }
    
    nFound = 0;
    ullpt = variableSiteIndexes;
    for ( i = 0; i < nVariableSites; i++ ) {
        j = *ullpt; // focal site master index
        cf_siteClass = *(siteClassifications + j);
        if ( cf_siteClass > SITE_CLASS_NEUTRAL ) {
            cf_siteClasses[nFound] = cf_siteClass;
            selectedSiteIndexesLocal[nFound] = i;
            selectedSiteIndexesMaster[nFound] = j;
            cf_scoeffs[nFound] = *(selectionCoefficients + j);
            nFound++;
        }
        ullpt++;
    }
    
    if ( VERBOSE ) {
        printf("\nLocal and Master selected site indexes, site type codes, and selection coefficients:\n");
        for ( i = 0; i < nSelectedSites; i++ )
            printf("\t%llu\t%llu\t%i\t%f\n", selectedSiteIndexesLocal[i], selectedSiteIndexesMaster[i], cf_siteClasses[i], cf_scoeffs[i]);
    }
    
    for ( i = 0; i < nSelectedSites; i++ ) {
        sipt = gts + (PLOIDY * selectedSiteIndexesLocal[i]); // first individual, selected site i
        dpt = fitnessValues; // first individual's fitness score
        cf_siteClass = cf_siteClasses[i];
        cf_s = cf_scoeffs[i];
        for ( j = 0; j < N; j++ ) {
            locpt = locations;
            gtsum = *sipt + *(sipt + 1); // number of derived alleles at this locus
            if ( (cf_siteClass == SITE_CLASS_POS || cf_siteClass == SITE_CLASS_BGS) && gtsum ) { // + and - selection handled same
                if ( FITNESS_MODEL == FITNESS_MODEL_ADDITIVE )
                    *dpt += ((double) gtsum) * cf_s;
                else
                    *dpt *= 1.0 + ((double) gtsum) * cf_s;
            }
            else if ( cf_siteClass == SITE_CLASS_DIV ) { // divergent selection requires special handling.
                cf_gradient = *( environmentGradient + *locpt );
                if ( gtsum == 2 ) { // homozygous derived
                    if ( FITNESS_MODEL == FITNESS_MODEL_ADDITIVE )
                        *dpt += 2.0 * cf_s * cf_gradient;
                    else
                        *dpt *= 1.0 + (2.0 * cf_s * cf_gradient);
                }
                else if ( gtsum == 0 ) { // homozygous ancestral
                    if ( FITNESS_MODEL == FITNESS_MODEL_ADDITIVE )
                        *dpt += -2.0 * cf_s * cf_gradient;
                    else
                        *dpt *= 1.0 - (2.0 * cf_s * cf_gradient);
                }
                else if ( gtsum != 1 ) {
                    fprintf(stderr, "\nError in computeFitness():\n\tgtsum (= %i) out of bounds.\n", gtsum);
                    exit(-1);
                    // heterozygotes are assumed intermediate everywhere, so no fitness change
                }
            }
            else if ( cf_siteClasses == SITE_CLASS_NEUTRAL ) {
                fprintf(stderr, "\nError in computeFitness():\n\tneutral site made it into calcs.\n");
                exit(-1);
            }
            
            sipt += PLOIDY * nVariableSites; // advance to locus i in next individual's genotype
            dpt++; // advance to next individual's fitness
            locpt++; // advance to next individual's location
        }
    }
    
    if ( t == 1 ) {
        FILE *cf_fitvals;
        cf_fitvals = fopen("InitialFitnessValues.csv", "w");
        fprintf(cf_fitvals, "individual,location,fitness\n");
        for ( i = 0; i < N; i++ )
            fprintf(cf_fitvals, "%llu,%i,%f\n", i, *(locations + i), *(fitnessValues + i));
        fclose(cf_fitvals);
    }
    dpt = fitnessValues;
    for ( i = 0; i < N; i++ ) {
        if ( *dpt < 0.0 )
            *dpt = 0.0;
        dpt++;
    }
}


void finalTasks(void)
{
    free(genotypes0);
    free(genotypes1);
    free(variableSiteIndexes);
    free(siteIndexes);
    free(siteIsVariable);
    free(alleleFrequencies);
    free(alleleCounts);
    free(SFScounts);
    free(siteClassifications);
    free(DEMOGRAPHIC_CHANGE_TIMES);
    free(K_VALUES);
    free(MIGRATION_CHANGE_TIMES);
    free(M_VALUES);
    free(abundances);
    free(selectionCoefficients);
    free(locations);
    free(environmentGradient);
    
    printf("\nPlease check code to make sure all malloc() are freed\n");
}


void initializeRNG(unsigned int RNG_SEED)
{
    int i;

    // initialize random number generator
    // T and r are global variables for the RNG state and calls
    gsl_rng_env_setup(); // set up the environment variables for the RNG
    rngType = gsl_rng_mt19937; // Mersenne Twister// for default you can say T = gsl_rng_default;
    rngState = gsl_rng_alloc(rngType);
    gsl_rng_set(rngState, RNG_SEED);
    
//    for ( i = 0; i < 10; i++ )
//        printf("%f\n", gsl_rng_uniform(rngState) );
}



void migration(void)
{
    long int i, j, counter, focalTime, nhere, nmoving, mover;
    double cumulativeM, numberLine[nPOPULATIONS], dum, *dpt;
    long int choiceVector[N], thoseMoving[N], immigrants[nPOPULATIONS], emigrants[nPOPULATIONS];
    int chosen, dumi;
    
    for ( i = 0; i < N; i++ )
        choiceVector[i] = i;
    for ( i = 0; i < nPOPULATIONS; i++ ) {
        immigrants[i] = 0;
        emigrants[i] = 0;
    }
    
    // check to see if now is a time to change migration rates
    if ( nMIGRATION_CHANGES > 0 ) {
        focalTime = MIGRATION_CHANGE_TIMES[currentMigrationPeriod];
        if ( t == focalTime ) {
            migRatePt += ( nPOPULATIONS * nPOPULATIONS );
            currentMigrationPeriod++;
        }
    }
    
    counter = 0;
    for ( i = 0; i < nPOPULATIONS; i++ ) {
        nhere = abundances[i]; // abundance in subpopulation i
        if ( nhere > 0 ) {
            cumulativeM = *(migRatePt + (i * nPOPULATIONS)); // first entry in ith "row" of current migration matrix
            numberLine[0] = cumulativeM;
            for ( j = 1; j < nPOPULATIONS; j++ ) {
                cumulativeM += (*(migRatePt + (i * nPOPULATIONS) + j));
                numberLine[j] = cumulativeM;
            }
            for ( j = 0; j < nPOPULATIONS; j++ )
                numberLine[j] = numberLine[j] / cumulativeM; // normalize
            
            if ( cumulativeM > 0.0 ) {
                
                // determine how many migrate: stochastic
                nmoving = gsl_ran_binomial(rngState, cumulativeM, nhere);
                
                if ( nmoving > 0 ) {
                    // determine which ones migrate: stochastic
                    // int gsl_ran_choose (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size)
                    gsl_ran_choose( rngState, thoseMoving, nmoving, choiceVector, nhere, sizeof(long int) );
                    
                    for ( j = 0; j < nmoving; j++ ) {
                        dum = gsl_rng_uniform( rngState );
                        dpt = numberLine;
                        dumi = 0;
                        while ( (*dpt) < dum && dumi < nPOPULATIONS) {
                            dpt++;
                            dumi++;
                        }
                        
                        if ( i == dumi ) {
                            fprintf(stderr, "\nError in migration():\n\tmigrator moving to same patch it is already in.\n\tNeed to fix algorithm\n");
                            exit(-1);
                        }
                        
                        mover = thoseMoving[j];
                        if ( *(locations + counter + mover) != i ) {
                            fprintf(stderr, "\nError in migration():\n\tYour mover isn't where you think she is.\n\tNeed to fix algorithm\n");
                            exit(-1);
                        }
                        *(locations + counter + mover) = dumi;
                        
                        immigrants[dumi] = immigrants[dumi] + 1;
                        emigrants[i] = emigrants[i] + 1;
                        
                    }
                    
                }
            }
        }
        counter += nhere;
    }
    
    if ( VERBOSE ) {
        printf("\nAbundances before migration:\n");
        for ( i = 0; i < nPOPULATIONS; i++ )
            printf("\t%li", abundances[i]);
        printf("\nEmigrants:\n");
        for ( i = 0; i < nPOPULATIONS; i++ )
            printf("\t%li", emigrants[i]);
        printf("\nImmigrants:\n");
        for ( i = 0; i < nPOPULATIONS; i++ )
            printf("\t%li", immigrants[i]);
    }
    
    for ( i = 0; i < nPOPULATIONS; i++ ) {
        abundances[i] = abundances[i] + immigrants[i] - emigrants[i];
    }
    
    if ( VERBOSE ) {
        printf("\nAbundances after migration:\n");
        for ( i = 0; i < nPOPULATIONS; i++ )
            printf("\t%li", abundances[i]);
    }
    
    
    //printf("\nWarning: migration() not written yet!\n");
    //exit(0);
}


void printParametersToFiles(unsigned RNG_SEED)
{
    
}


double randExp(double meanValue)
{
    return ( (log(1.0 - gsl_rng_uniform(rngState))) * (-meanValue) );
}


unsigned readInParametersFromFile(void)
{
    char c, option[80], option2[80];
    long int INITIAL_N;
    int i, j, k, dumi, temp;
    unsigned RNG_SEED = 1;
    double value;
    FILE *pfile;
    _Bool nPopSet = 0, migrationSet = 0, demographySet = 0;
    
    pfile = fopen("parameters.ini.txt", "r");
    if ( pfile == NULL ) {
        fprintf(stderr, "\nError in readInParametersFromFile():\n\tparameters.ini.txt not found!\n");
        exit(-1);
    }
    
    // the action of reading.  Order of things in parameters.ini.txt file matters!
    while ( !feof(pfile) ) {
        
        fscanf(pfile, "%s", option);
        
        if ( !strcmp( option, "RNG_SEED" ) ) { // set RNG seed value
            fscanf(pfile, "%i", &RNG_SEED);
            if ( VERBOSE )
                printf("Found RNG_SEED (%s) = %u\n", option, RNG_SEED);
        }
        else if ( !strcmp( option, "nGENERATIONS" ) ) { // set RNG seed value
            fscanf(pfile, "%li", &nGENERATIONS);
            if ( VERBOSE )
                printf("Found nGENERATIONS (%s) = %li\n", option, nGENERATIONS);
        }
        else if ( !strcmp( option, "nPOPULATIONS" ) ) { // set number of populations
            fscanf(pfile, "%i", &nPOPULATIONS);
            nPopSet = 1;
            if ( nPOPULATIONS < 0 ) {
                fprintf(stderr, "\nError in readInParametersFromFile():\n\t");
                fprintf(stderr, "nPOPULATIONS (%i) < 0\n\t", nPOPULATIONS);
                fprintf(stderr, "*** Exiting *** \n");
                exit(-1);
            }
            if ( VERBOSE )
                printf("Found nPOPULATIONS (%s) = %i\n", option, nPOPULATIONS);
        }
        else if ( !strcmp( option, "nDEMOGRAPHIC_CHANGES" ) ) { // set demography
            if ( !nPopSet ) {
                fprintf(stderr, "\nWarning in readInParametersFromFile():\n\t");
                fprintf(stderr, "nPOPULATIONS should be set in parameters.ini.txt before nDEMOGRAPHIC_CHANGES\n\t");
            }
            
            demographySet = 1;
            
            // number of demographic events:
            fscanf(pfile, "%i", &nDEMOGRAPHIC_CHANGES);
            if ( VERBOSE)
                printf("Found nDEMOGRAPHIC_CHANGES (%s) = %i\n", option, nDEMOGRAPHIC_CHANGES);
            
            // times of demographic events:
            fscanf(pfile, "%s", option);
            if ( strcmp( option, "DEMOGRAPHIC_CHANGE_TIMES" ) ) {
                wrongParametersIniOption( "DEMOGRAPHIC_CHANGE_TIMES", "nDEMOGRAPHIC_CHANGES", option );
            }
            else {
                if ( nDEMOGRAPHIC_CHANGES > 0 ) {
                    DEMOGRAPHIC_CHANGE_TIMES = (long int *) malloc( nDEMOGRAPHIC_CHANGES * sizeof(long int));
                    for ( i = 0; i < nDEMOGRAPHIC_CHANGES; i++ ) {
                        fscanf(pfile, "%li", (DEMOGRAPHIC_CHANGE_TIMES + i));
                    }
                }
                else {
                    DEMOGRAPHIC_CHANGE_TIMES = (long int *) malloc( sizeof(long int) );
                    *DEMOGRAPHIC_CHANGE_TIMES = 0;
                }
            }
            
            // K values for demographic event periods:
            fscanf(pfile, "%s", option);
            if ( strcmp( option, "K_VALUES" ) ) {
                wrongParametersIniOption( "K_VALUES", "DEMOGRAPHIC_CHANGE_TIMES", option );
            }
            else {
                INITIAL_N = 0;
                K_VALUES = (double *) malloc( (nDEMOGRAPHIC_CHANGES + 1) * nPOPULATIONS * sizeof(double) );
                for ( i = 0; i < ((nDEMOGRAPHIC_CHANGES + 1) * nPOPULATIONS); i++ ) {
                    fscanf(pfile, "%lf", (K_VALUES + i));
                    if ( i < nPOPULATIONS )
                        INITIAL_N += K_VALUES[i];
                    if ( *(K_VALUES+i) < 0.0 ) {
                        fprintf(stderr, "\nError in readInParametersFromFile():\n\t");
                        fprintf(stderr, "K_VALUES should be non-negative, but instead found %iith value = %f.\n\tPlease fix parameters.ini.txt\n", (i+1), *(K_VALUES+i));
                        fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
                        exit(-1);
                    }
                }
            }
            // test check
            if ( VERBOSE ) {
                printf("Found DEMOGRAPHIC_CHANGE_TIMES: ");
                for ( i = 0; i < nDEMOGRAPHIC_CHANGES; i++ ) {
                    printf(" %li", DEMOGRAPHIC_CHANGE_TIMES[i]);
                }
                printf("\nINITIAL_N = %li\nFound the following K_VALUES:\n", INITIAL_N);
                for ( i = 0; i < ((nDEMOGRAPHIC_CHANGES + 1) * nPOPULATIONS); i++ ) {
                    printf("\t%f", K_VALUES[i]);
                    if ( i % nPOPULATIONS == nPOPULATIONS - 1 )
                        printf("\n");
                }
            }
            

        }
        else if ( !strcmp( option, "FIXED_POP_SIZE"  ) ) { // whether or not ppulation sizes are fixed
            fscanf(pfile, "%i", &temp);
            FIXED_POP_SIZE = temp;
            if ( VERBOSE )
                printf("Found FIXED_POP_SIZE (%s) = %i\n", option, FIXED_POP_SIZE);
        }
        else if ( !strcmp( option, "MAX_POP_GROWTH_RATE"  ) ) { // maximum population growth rate in logistic equation
            fscanf(pfile, "%lf", &MAX_POP_GROWTH_RATE);
            if ( FIXED_POP_SIZE ) {
                fprintf(stderr, "\nWarning in readInParametersFromFile():\n\t");
                fprintf(stderr, "Since population size is fixed, MAX_POP_GROWTH_RATE won't be used\n\t");
            }
            else if ( VERBOSE )
                printf("Found MAX_POP_GROWTH_RATE (%s) = %f\n", option, MAX_POP_GROWTH_RATE);
        }
        else if ( !strcmp( option, "nMIGRATION_CHANGES" ) ) { // set up migration rates
            if ( !nPopSet ) {
                fprintf(stderr, "\nWarning in readInParametersFromFile():\n\t");
                fprintf(stderr, "nPOPULATIONS should be set in parameters.ini.txt before nMIGRATION_CHANGES\n\t");
            }
            migrationSet = 1;
            
            // number of migration change events:
            fscanf(pfile, "%i", &nMIGRATION_CHANGES);
            if ( VERBOSE )
                printf("Found nMIGRATION_CHANGES (%s) = %i\n", option, nMIGRATION_CHANGES);
            
            // times of migration change events:
            fscanf(pfile, "%s", option);
            if ( strcmp( option, "MIGRATION_CHANGE_TIMES" ) ) {
                wrongParametersIniOption( "MIGRATION_CHANGE_TIMES", "nMIGRATION_CHANGES", option);
            }
            else {
                if ( nMIGRATION_CHANGES > 0 ) {
                    MIGRATION_CHANGE_TIMES = (long int *) malloc( nMIGRATION_CHANGES * sizeof(long int));
                    for ( i = 0; i < nMIGRATION_CHANGES; i++ ) {
                        fscanf(pfile, "%li", (MIGRATION_CHANGE_TIMES + i));
                    }
                }
                else {
                    MIGRATION_CHANGE_TIMES = (long int *) malloc( sizeof(long int) );
                    *MIGRATION_CHANGE_TIMES = 0;
                }
            }
            
            // migration rate values for periods.  expecting nPOPULATIONS x nPOPULATIONS entries (like a matrix) for each period
            fscanf(pfile, "%s", option);
            if ( strcmp( option, "M_VALUES" ) ) {
                wrongParametersIniOption( "M_VALUES", "MIGRATION_CHANGE_TIMES", option);
            }
            else {
                M_VALUES = (double *) malloc( nPOPULATIONS * nPOPULATIONS * (nMIGRATION_CHANGES + 1) * sizeof(double) );
                dumi = 0;
                for ( i = 0; i < (nMIGRATION_CHANGES+1); i++ ) {
                    for ( j = 0; j < nPOPULATIONS; j++ ) {
                        for ( k = 0; k < nPOPULATIONS; k++ ) {
                            fscanf(pfile, "%lf", (M_VALUES + dumi));
                            value =  *(M_VALUES + dumi);
                            if ( value < 0.0 || value > 1.0 ) {
                                fprintf(stderr, "\nError in readInParametersFromFile():\n\t");
                                fprintf(stderr, "found migration rate value out of bounds [0,1] = %f\n\t", value);
                                fprintf(stderr, "Please fix parameters.ini.txt\n");
                                fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
                                exit(-1);
                            }
                            else if ( j == k && value > 0.0 ) {
                                fprintf(stderr, "\nError in readInParametersFromFile():\n\t");
                                fprintf(stderr, "found non-zero migration rate on diagonal = %f\n\t", value);
                                fprintf(stderr, "Please fix parameters.ini.txt\n");
                                fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
                                exit(-1);
                            }
                            dumi++;
                        }
                    }
                }
            }
            // test check
            if ( VERBOSE ) {
                printf("Found MIGRATION_CHANGE_TIMES: ");
                for ( i = 0; i < nMIGRATION_CHANGES; i++ ) {
                    printf(" %li", MIGRATION_CHANGE_TIMES[i]);
                }
                printf("\nFound M_VALUES:\n");
                dumi = 0;
                for ( i = 0; i < (nMIGRATION_CHANGES+1); i++ ) {
                    for ( j = 0; j < nPOPULATIONS; j++ ) {
                        for ( k = 0; k < nPOPULATIONS; k++ ) {
                            printf("\t%f", M_VALUES[dumi++]);
                        }
                        printf("\n");
                    }
                    printf("\n");
                }
            }
            
            
        }
        else if ( !strcmp( option, "nSITES"  ) ) { // total number of sites in genome
            fscanf(pfile, "%llu", &nSITES);
            if ( nSITES < 1 ) {
                fprintf(stderr, "\nError in readInParametersFromFile():\n\tnSITES (= %llu) should be > 1\n", nSITES);
                fprintf(stderr, "Please fix parameters.ini.txt\n");
                fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
                exit(-1);
            }
            // test check
            if ( VERBOSE )
                printf("Found nSITES (%s) = %llu\n", option, nSITES);
        }
        else if ( !strcmp( option, "nLINKAGE_GROUPS"  ) ) { // total number of sites in genome
            fscanf(pfile, "%i", &nLINKAGE_GROUPS);
            if ( nLINKAGE_GROUPS < 1 ) {
                fprintf(stderr, "\nError in readInParametersFromFile():\n\nLINKAGE_GROUPS (= %i) should be > 1\n", nLINKAGE_GROUPS);
                fprintf(stderr, "Please fix parameters.ini.txt\n");
                fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
                exit(-1);
            }
            // test check
            if ( VERBOSE )
                printf("Found nLINKAGE_GROUPS (%s) = %i\n", option, nLINKAGE_GROUPS);
        }
        else if ( !strcmp( option, "RECOMBINATION_RATE_PER_KB"  ) ) { // background selection mean coefficient
            fscanf(pfile, "%lf", &RECOMBINATION_RATE_PER_KB);
            if ( VERBOSE )
                printf("Found RECOMBINATION_RATE_PER_KB (%s) = %f\n", option, RECOMBINATION_RATE_PER_KB);
        }
        else if ( !strcmp( option, "INCLUDE_SELECTION"  ) ) { // using selection
            fscanf(pfile, "%i", &temp);
            INCLUDE_SELECTION = temp;
            if ( VERBOSE )
                printf("Found INCLUDE_SELECTION (%s) = %i\n", option, INCLUDE_SELECTION);
        }
        else if ( !strcmp( option, "MEAN_S_BGS"  ) ) { // background selection mean coefficient
            fscanf(pfile, "%lf", &MEAN_S_BGS);
            if ( VERBOSE )
                printf("Found MEAN_S_BGS (%s) = %f\n", option, MEAN_S_BGS);
        }
        else if ( !strcmp( option, "MEAN_S_POS"  ) ) { // positive selection mean coefficient
            fscanf(pfile, "%lf", &MEAN_S_POS);
            if ( VERBOSE )
                printf("Found MEAN_S_POS (%s) = %f\n", option, MEAN_S_POS);
        }
        else if ( !strcmp( option, "MEAN_S_DIV"  ) ) { // divergent selection mean coefficient
            fscanf(pfile, "%lf", &MEAN_S_DIV);
            if ( VERBOSE )
                printf("Found MEAN_S_DIV (%s) = %f\n", option, MEAN_S_DIV);
        }
        else if ( !strcmp( option, "PROBABILITY_SITE_BGS"  ) ) { // divergent selection mean coefficient
            fscanf(pfile, "%lf", &PROBABILITY_SITE_BGS);
            if ( VERBOSE )
                printf("Found PROBABILITY_SITE_BGS (%s) = %f\n", option, PROBABILITY_SITE_BGS);
        }
        else if ( !strcmp( option, "PROBABILITY_SITE_POS"  ) ) { // divergent selection mean coefficient
            fscanf(pfile, "%lf", &PROBABILITY_SITE_POS);
            if ( VERBOSE )
                printf("Found PROBABILITY_SITE_POS (%s) = %f\n", option, PROBABILITY_SITE_POS);
        }
        else if ( !strcmp( option, "PROBABILITY_SITE_DIV"  ) ) { // divergent selection mean coefficient
            fscanf(pfile, "%lf", &PROBABILITY_SITE_DIV);
            if ( VERBOSE )
                printf("Found PROBABILITY_SITE_DIV (%s) = %f\n", option, PROBABILITY_SITE_DIV);
        }
        else if ( !strcmp( option, "ENVIRONMENT_TYPE"  ) ) { // divergent selection mean coefficient
            fscanf(pfile, "%s", option2);
            if ( !strcmp( option2, "GRADIENT" ) )
                ENVIRONMENT_TYPE = ENVT_TYPE_GRADIENT;
            else if ( !strcmp( option2, "MOSAIC" ) )
                ENVIRONMENT_TYPE = ENVT_TYPE_MOSAIC;
            else if ( !strcmp( option2, "INVARIANT" ) )
                ENVIRONMENT_TYPE = ENVT_TYPE_INVARIANT;
            else {
                fprintf(stderr, "\nError in readInParametersFromFile():\n\tENVIRONMENT_TYPE option (%s) not recognized.", option2);
                fprintf(stderr, "\n\tUsable options are 'GRADIENT', 'MOSAIC', and 'INVARIANT'.\n\tPlease fix parameters.ini.txt\n");
                fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
                exit(-1);
            }
            if ( VERBOSE )
                printf("Found ENVIRONMENT_TYPE (%s) = %s = %i\n", option, option2, ENVIRONMENT_TYPE);
        }
        else if ( !strcmp( option, "ENVT_MIN"  ) ) { // divergent selection mean coefficient
            fscanf(pfile, "%lf", &ENVT_MIN);
            if ( VERBOSE )
                printf("Found ENVT_MIN (%s) = %f\n", option, ENVT_MIN);
        }
        else if ( !strcmp( option, "ENVT_MAX"  ) ) { // divergent selection mean coefficient
            fscanf(pfile, "%lf", &ENVT_MAX);
            if ( VERBOSE )
                printf("Found ENVT_MAX (%s) = %f\n", option, ENVT_MAX);
        }
        else if ( !strcmp( option, "FITNESS_MODEL"  ) ) { // divergent selection mean coefficient
            fscanf(pfile, "%s", option2);
            if ( !strcmp( option2, "ADDITIVE" ) )
                FITNESS_MODEL = FITNESS_MODEL_ADDITIVE;
            else if ( !strcmp( option2, "MULTIPLICATIVE" ) )
                FITNESS_MODEL = FITNESS_MODEL_MULTIPLICATIVE;
            else {
                fprintf(stderr, "\nError in readInParametersFromFile():\n\tFITNESS_MODEL option (%s) not recognized.", option2);
                fprintf(stderr, "\n\tUsable options are 'ADDITIVE' and 'MULTIPLICATIVE'.\n\tPlease fix parameters.ini.txt\n");
                fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
                exit(-1);
            }
            if ( VERBOSE )
                printf("Found ENVIRONMENT_TYPE (%s) = %s = %i\n", option, option2, FITNESS_MODEL);
        }

        
        option[0] = '\0'; // reset to avoid double setting last option
        option2[0] = '\0';
    }
    
    fclose(pfile);
    
    PROBABILITY_SITE_NEUTRAL = 1.0 - (PROBABILITY_SITE_BGS + PROBABILITY_SITE_DIV + PROBABILITY_SITE_POS);
    if ( PROBABILITY_SITE_NEUTRAL < 0.0 || PROBABILITY_SITE_BGS < 0.0 || PROBABILITY_SITE_POS < 0.0 || PROBABILITY_SITE_DIV < 0.0 ) {
        fprintf(stderr, "\nError in readInParametersFromFile():\n\t");
        fprintf(stderr, "One or more PROBABILITY_SITE_... variables out of bounds\n\t");
        fprintf(stderr, "Please fix parameters.ini.txt\n");
        fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
        exit(-1);
    }
    
    if ( !demographySet ) {
        INITIAL_N = 0;
        nDEMOGRAPHIC_CHANGES = 0;
        DEMOGRAPHIC_CHANGE_TIMES = (long int *) malloc( sizeof(long int) );
        *DEMOGRAPHIC_CHANGE_TIMES = 0;
        K_VALUES = (double *) malloc( nPOPULATIONS * sizeof(double) );
        for ( i = 0; i < nPOPULATIONS; i++ ) {
            K_VALUES[i] = K_DEFAULT;
            INITIAL_N += K_DEFAULT;
        }
        if ( VERBOSE )
            printf("\nINITIAL_N = %li\n", INITIAL_N);
    }
    if ( !migrationSet ) {
        nMIGRATION_CHANGES = 0;
        MIGRATION_CHANGE_TIMES = (long int *) malloc( sizeof( long int ) );
        *MIGRATION_CHANGE_TIMES = 0;
        M_VALUES = (double *) malloc( nPOPULATIONS * nPOPULATIONS * sizeof(double) );
        k = 0;
        for ( i = 0; i < nPOPULATIONS; i++ ) {
            for ( j = 0; j < nPOPULATIONS; j++ ) {
                if ( i == j )
                    M_VALUES[k] = 0.0;
                else
                    M_VALUES[k] = MIGRATION_RATE_DEFAULT;
                k++;
            }
        }
    }
    KvalPt = K_VALUES;
    migRatePt = M_VALUES;
    
    
    /* printf("\nWarning: not done writing readInParametersFromFile().\nNeed to write code for cases when paramters are left out.\nAnd probably other stuff too.\n\n"); */
    
    //exit(0);
    N = INITIAL_N;
    printParametersToFiles(RNG_SEED);
    
    return RNG_SEED;
}



void reproduction(void)
{
    int pop;
    long int i, j, noffspring[nPOPULATIONS], totalOffspring;
    short int *offspringGTs, *sipt;
    
    if ( INCLUDE_SELECTION ) {
        double fitnessValues[N];
        computeFitness(fitnessValues);
    }
    
    totalOffspring = 0;
    for ( pop = 0; pop < nPOPULATIONS; pop++ ) {
        noffspring[pop] = calculateNumOffspring(pop);
        totalOffspring += noffspring[pop];
    }
    
    if ( VERBOSE ) {
        printf("\nNumber of offspring in each deme:\n");
        for ( i = 0; i < nPOPULATIONS; i++ )
            printf("\t%li", noffspring[i]);
        printf("\n");
    }

    // now, get a memory block to use to store the offspring genotypes
    offspringGTs = checkMemoryBlocks(totalOffspring);
    
    
    
    printf("\nWarning: reproduction() not finished yet!\n");
    exit(0);
}


void setUpGenome(void)
{
    double theta, expectedSegSites, *expectedFreq, *dpt, value;
    unsigned long long int foo, *ullpt;
    long int i;
    int dumi;
    char str[80];
    _Bool firstOne;
    
    // set up standing neutral variation using classical pop. gen. expectations
    // theta is the poulation mutation parameter from coalescent theory, specifically
    // here from Wakeley's book, p. 92, "... twice the expected number of mutations
    // introduced into the population each generation ... "
    foo = N * PLOIDY * nSITES;
    theta = 2.0 * MU * ((double) foo);
    
    // now the expected number of segregating sites, Wakeley, p. 97, eqn. 4.7:
    // note for diploid, we have 2N sequences in the population, so the sum is from 1 to 2N-1;
    // the total number of expected segregating sites is the sum of the expected number of sites with each possible
    // number of copies of a derived allele
    expectedFreq = (double *) malloc( PLOIDY * N * sizeof(double) );
    *expectedFreq = 0.0;
    dpt = expectedFreq + 1; // start with expected number of sites having 1 copy (not with zero copies)
    expectedSegSites = 0.0;
    for ( i = 1; i < (PLOIDY * N); i++ ) {
        expectedSegSites += 1.0 / ((double) i); // summation equation; below sum is multiplied by theta
        *dpt = (*(dpt - 1)) + theta / ((double) i); // equation from Wakeley's book (based on Fu 1995); build cumulative vector
        dpt++;
    }
    expectedSegSites = expectedSegSites * theta;
    nVariableSites = (unsigned long long int) (expectedSegSites + 0.5);
    
    
    // now to choose the sites that are variable and their frequencies
    siteIndexes = (unsigned long long int *) malloc( nSITES * sizeof(unsigned long long int) );
    variableSiteIndexes = (unsigned long long int *) malloc( nVariableSites * sizeof(unsigned long long int) );
    ullpt = siteIndexes;
    for ( foo = 0; foo < nSITES; foo++ ) {
        *ullpt = foo;
        ullpt++;
    }
    
    // use gsl ran choose to pick sites at random to be standing neutral variation
    gsl_ran_choose( rngState, variableSiteIndexes, nVariableSites, siteIndexes, nSITES, sizeof(unsigned long long int) );
    // variable site indexes now stores the ordered indexes of sites that should have segregating variants
    // siteIndexes just stores the indexes [0..nSITES-1], which is useful for random choices like this one, so we keep it
    siteIsVariable = (_Bool *) malloc( nSITES * sizeof(_Bool) );
    siteClassifications = (int *) malloc( nSITES * sizeof(int) );
    memset( siteIsVariable, 0, (nSITES * sizeof(_Bool)) );
    for ( foo = 0; foo < nVariableSites; foo++ )
        *(siteIsVariable + (variableSiteIndexes[foo])) = 1;
    // siteIsVariable stores the same information as variableSiteIndexes, but in a different form; both are useful
    
    // now randomly assign frequencies
    alleleFrequencies = (double *) malloc( nSITES * sizeof(double)); // frequencies (0 to 1) by site
    alleleCounts = (unsigned long long int *) malloc( nSITES * sizeof(unsigned long long int));
    memset( alleleCounts, 0, nSITES * sizeof(unsigned long long int) );
    SFScounts = (unsigned long long int *) malloc( PLOIDY * N * sizeof(unsigned long long int) );
    setUpInitialAlleleFrequencies(expectedFreq);
    
    // pick allele frequencies using the Ewens sampling formula from Wakeley's book.

    
    if ( VERBOSE )
        printf("\ntheta = %f, segSites = %llu\n", theta, nVariableSites);

    
    // assign site types
    double siteProbs[nSITE_TYPES];
    unsigned long long int siteTypeCounts[nSITE_TYPES];
    memset( siteTypeCounts, 0, (nSITE_TYPES * sizeof(unsigned long long int)) );
    selectionCoefficients = (double *) malloc( nSITES * sizeof(double) );
    memset( selectionCoefficients, 0, (nSITES * sizeof(double)) );
    siteProbs[0] = PROBABILITY_SITE_NEUTRAL;
    siteProbs[1] = siteProbs[0] + PROBABILITY_SITE_BGS;
    siteProbs[2] = siteProbs[1] + PROBABILITY_SITE_POS;
    siteProbs[3] = 1.0;
    for ( foo = 0; foo < nSITES; foo ++ ) {
        value = gsl_rng_uniform(rngState);
        i = -1;
        do {
            i++;
        } while ( value > siteProbs[i] );
        if ( i == 0 ) {
            *(siteClassifications + foo) = SITE_CLASS_NEUTRAL;
            *(selectionCoefficients + foo) = 0.0;
        }
        else if ( i == 1 ) {
            *(siteClassifications + foo) = SITE_CLASS_BGS;
            *(selectionCoefficients + foo) = -(randExp(MEAN_S_BGS)) * CODOMINANCE;
        }
        else if ( i == 2 ) {
            *(siteClassifications + foo) = SITE_CLASS_POS;
            *(selectionCoefficients + foo) = randExp(MEAN_S_POS) * CODOMINANCE;
        }
        else if ( i == 3 ) {
            *(siteClassifications + foo) = SITE_CLASS_DIV;
            *(selectionCoefficients + foo) = randExp(MEAN_S_DIV) * CODOMINANCE;
        }
        else {
            fprintf(stderr, "\nError in setUpGenome():\n\tbogus site type (%li)\n\t*** Exiting ***\n\n", i);
            exit(-1);
        }
        siteTypeCounts[i] = siteTypeCounts[i] + 1;
    }
    
    foo = 0;
    if ( VERBOSE )
        printf("\nSite type counts:\n");
    for ( i = 0; i < nSITE_TYPES; i++ ) {
        if ( VERBOSE )
            printf("\t%llu", siteTypeCounts[i]);
        foo += siteTypeCounts[i];
    }
    if ( VERBOSE )
        printf("\n");
    if ( foo != nSITES ) {
        fprintf(stderr, "\nError in setUpGenome():\n\tsite total (%llu) != nSITES (%llu)\n\t*** Exiting ***\n\n", foo, nSITES);
        exit(-1);
    }
    
    // save data on selected sites to files
    FILE *siteDesignations, *sd2;
    siteDesignations = fopen("SiteClassIndexes.R", "w"); // R format
    sd2 = fopen("SiteClassIndexes.m", "w"); // matlab format
    for ( i = 1; i < nSITE_TYPES; i++ ) {
        if ( i == SITE_CLASS_BGS ) {
            dumi = SITE_CLASS_BGS;
            strcpy( str, "BGS_SITE_INDEXES" );
        }
        else if ( i == SITE_CLASS_POS ) {
            dumi = SITE_CLASS_POS;
            strcpy( str, "POS_SITE_INDEXES" );
        }
        else if ( i == SITE_CLASS_DIV ) {
            dumi = SITE_CLASS_DIV;
            strcpy( str, "DIV_SITE_INDEXES" );
        }
        else {
            fprintf(stderr, "\nError in setUpGenome():\n\tsite designation out of bounds\n");
            exit(-1);
        }
        fprintf( siteDesignations, "%s <- c(", str );
        fprintf( sd2, "%s = [", str );
        firstOne = 1;
        for ( foo = 0; foo < nSITES; foo++ ) {
            if ( *(siteClassifications + foo) == dumi ) {
                if ( firstOne ) { // so no trailing comma or space on end
                    fprintf( siteDesignations, "%llu", foo );
                    fprintf( sd2, "%llu", foo );
                    firstOne = 0;
                }
                else {
                    fprintf( siteDesignations, ",%llu", foo );
                    fprintf( sd2, " %llu", foo );
                }
            }
        }
        fprintf( siteDesignations, ");\n\n" );
        fprintf( sd2, "];\n\n" );
    }
    fclose(siteDesignations);
    fclose(sd2);
    
    nSelectedSites = 0;
    FILE *initialFreqs;
    unsigned long long int focalSiteIndex;
    int SiteClassCode;
    initialFreqs = fopen("InitialAlleleFreqs.csv", "w");
    fprintf(initialFreqs, "SiteIndex,DerivedAlleleCount,DerivedAlleleFreq,SelectionCoefficient,SiteClassCode,SiteClassName\n");
    for ( foo = 0; foo < nVariableSites; foo++ ) {
        focalSiteIndex = *(variableSiteIndexes + foo);
        
        if ( *(siteClassifications + focalSiteIndex) > SITE_CLASS_NEUTRAL )
            nSelectedSites++;

        
        fprintf(initialFreqs, "%llu,%llu,%E,%E", focalSiteIndex, alleleCounts[focalSiteIndex], alleleFrequencies[focalSiteIndex], selectionCoefficients[focalSiteIndex] );
        SiteClassCode = siteClassifications[focalSiteIndex];
        fprintf(initialFreqs, ",%i", SiteClassCode);
        if ( SiteClassCode == SITE_CLASS_BGS )
            fprintf(initialFreqs, ",SITE_CLASS_BGS\n");
        else if ( SiteClassCode == SITE_CLASS_POS )
            fprintf(initialFreqs, ",SITE_CLASS_POS\n");
        else if ( SiteClassCode == SITE_CLASS_DIV )
            fprintf(initialFreqs, ",SITE_CLASS_DIV\n");
        else if ( SiteClassCode == SITE_CLASS_NEUTRAL )
            fprintf(initialFreqs, ",SITE_CLASS_NEUTRAL\n");
        else {
            fprintf(stderr, "\nError in setUpGenome:\n\tSite class (%i) not found\n", SiteClassCode);
            exit(-1);
        }
    }
    fclose(initialFreqs);
    
    if ( VERBOSE )
        printf("\nnSelectedSites = %li\n", nSelectedSites);
    
    initialFreqs = fopen("InitialSFS.csv", "w");
    fprintf(initialFreqs, "AlleleCount,NumberOfSites\n");
    for ( foo = 0; foo < (PLOIDY * N); foo++ ) {
        if ( SFScounts[foo] ) {
            fprintf(initialFreqs, "%llu,%llu\n", foo, SFScounts[foo]);
        }
    }
    fclose(initialFreqs);
    
    
    
    
    free(expectedFreq);
    
    //printf("\nWarning: setUpGenome() not done yet.  Still need to set selection coefficients.\n");
    //exit(0);
}


void setUpInitialAlleleFrequencies(double *expectedFreq)
{
    unsigned long long int i, copies, focalSiteIndex;
    double dum, *dpt, twoN;
    _Bool looking;
    
    twoN = ((double) PLOIDY) * ((double) N);
    
    // make proper scaling of expectedFreq
    dpt = expectedFreq + 1;
    dum = *(expectedFreq + (PLOIDY * N) - 1); // cumulative sum
    for ( i = 1; i < ((PLOIDY * N) - 1); i++ ) {
        *dpt = (*dpt) / dum; // rescale so total sums to one
        dpt++;
//        if ( i < 20 )
//            printf("%E\t%E\n", *(dpt-2), *(dpt-1));
    }
    *(expectedFreq + (PLOIDY * N) - 1) = 1.1; // top out last entry
    
    
    for ( i = 0; i < nVariableSites; i++ ) {
        
        // get index of focal site
        focalSiteIndex = *(variableSiteIndexes + i);
        
        // draw random uniform number to choose frequency
        dum = gsl_rng_uniform(rngState);
        looking = 1;
        copies = 1;
        dpt = expectedFreq + 1;
        // run through cumulative probability vector of copy number find chosen one
        while ( looking && (copies < ((PLOIDY * N) - 1)) ) {
            if ( dum < (*dpt) )
                looking = 0;
            else {
                copies++;
                dpt++;
            }
        }
        
        // store results
        *(alleleFrequencies + focalSiteIndex) = ((double) copies) / twoN;
        *(alleleCounts + focalSiteIndex) = copies;
        *(SFScounts + copies) = (*(SFScounts + copies)) + 1;
    }
}



void setUpPopulations(void)
{
    unsigned long long int neededSize, focalSiteIndex;
    int i;
    unsigned long long int j, k, index, counter, siteAlleleCount;
    unsigned long long int choiceVector[(PLOIDY * N)], allelesToSwitch[(PLOIDY * N)];
    short int *sipt;
    double stepSize;
    
    // allocate space for genotypes and keep track of how big they are in memory
    neededSize = PLOIDY * N * sizeof(short int) * nVariableSites;
    genotypes0 = (short int *) malloc( neededSize );
    genotypes1 = (short int *) malloc( neededSize );
    currentBlock = 0;
    gts = genotypes0;
    blockSizes[0] = neededSize;
    blockSizes[1] = neededSize;
    memset(genotypes0, ALLELE_CODE_ANCESTRAL, neededSize);
    memset(genotypes1, ALLELE_CODE_ANCESTRAL, neededSize);
    
    // make array of abundances in subpopulations
    abundances = (long int *) malloc( nPOPULATIONS * sizeof( long int ) );
    for ( i = 0; i < nPOPULATIONS; i++ )
        abundances[i] = K_VALUES[i]; // initialize with populations at current carrying capacity
    
    // set locations of individuals
    locations = (int *) malloc( N * sizeof(int) );
    counter = 0;
    for ( i = 0; i < nPOPULATIONS; i++ ) {
        for ( j = 0; j < abundances[i]; j++ ) {
            locations[counter] = i;
            counter++;
        }
    }
    
    // make a choice vector for choosing which alleles in which genotypes to change from zeros to ones
    counter = 0;
    for ( j = 0; j < N; j++ ) { // individuals' genotypes
        for ( i = 0; i < PLOIDY; i++ ) { // the two copies each individual has
            /* the sequence of choices needs to be: 0,1, PLOIDY * nVariableSites, PLOIDY*nVariableSites + 1, 2*PLOIDY*nVariableSites, 2*PLOIDY*nVariableSites + 1, etc... */
            choiceVector[counter] = (j * PLOIDY * nVariableSites) + i;
            counter++;
        }
    }
    
    for ( j = 0; j < nVariableSites; j++ ) {
        focalSiteIndex = variableSiteIndexes[j];
        siteAlleleCount = alleleCounts[focalSiteIndex];
        if ( siteAlleleCount <= 0 || siteAlleleCount >= (PLOIDY * N) ) {
            printf("\nError in setUpPopulations():\n\tsiteAlleleCount (= %llu) out of bounds.\n", siteAlleleCount);
            exit(-1);
        }
            
        // int gsl_ran_choose (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size)
        gsl_ran_choose( rngState, allelesToSwitch, siteAlleleCount, choiceVector, (PLOIDY * N), sizeof(unsigned long long int) );
        
        for ( k = 0; k < siteAlleleCount; k++ ) {
            // index into genotypes array
            // genotypes array: each individual = PLOIDY * nVariableSites consecutive entries
            // two adjacent entries are the two entries for the two alleles an individual has at that locus
            // need to add "1" allele at jth locus for each of selected genotypes that get a derived allele copy
            // hence, index of jth locus in kth genotype copy is:
            index = (PLOIDY * j) + allelesToSwitch[k]; // kth genotype copy of jth locus to add allele to
            *(gts + index) = ALLELE_CODE_DERIVED; // derived allele
        }
        
    }
    
    
    FILE *initgts;
    initgts = fopen("InitialGenotypes.csv", "w");
    fprintf(initgts, "Locus0Copy0,Locus0Copy1");
    for ( k = 1; k < nVariableSites; k++ ) {
        for ( i = 0; i < PLOIDY; i++ ) {
            fprintf(initgts, ",Locus%lluCopy%i", k, i);
        }
    }
    fprintf(initgts,"\n");
    
    sipt = gts;
    for ( k = 0; k < N; k++ ) {
        fprintf(initgts, "%i", *sipt);
        sipt++;
        for ( j = 1; j < (PLOIDY * nVariableSites); j++ ) {
            fprintf(initgts, ",%i", *sipt);
            sipt++;
        }
        fprintf(initgts, "\n");
    }
    fclose(initgts);
    
    environmentGradient = (double *) malloc( nPOPULATIONS * sizeof(double) );
    if ( ENVIRONMENT_TYPE == ENVT_TYPE_GRADIENT )
        stepSize = (ENVT_MAX - ENVT_MIN) / (((double) nPOPULATIONS) - 1.0);
    for ( i = 0; i < nPOPULATIONS; i++ ) {
        if ( ENVIRONMENT_TYPE == ENVT_TYPE_GRADIENT )
            environmentGradient[i] = ENVT_MIN + (((double) i) * stepSize);
        else if ( ENVIRONMENT_TYPE == ENVT_TYPE_MOSAIC ) {
            if ( i % 2 )
                environmentGradient[i] = ENVT_MAX;
            else
                environmentGradient[i] = ENVT_MIN;
        }
        else if ( ENVIRONMENT_TYPE == ENVT_TYPE_INVARIANT )
            environmentGradient[i] = 0.0;
        else {
            fprintf(stderr, "\nError in setUpPopulations():\n\tENVIRONMENT_TYPE (%i) unrecognized\n", ENVIRONMENT_TYPE);
            exit(-1);
        }
    }
    if ( VERBOSE ) {
        printf("\nenvironmentGradient:\n");
        for ( i = 0; i < nPOPULATIONS; i++ )
            printf("\t%f", environmentGradient[i]);
        printf("\n");
    }
    
}



void usage(char *progname)
{
    fprintf(stdout, "\n%s:\n\tOptions:\n", progname);
    
    fprintf(stdout, "\nNote that limited options are available here on the command line.\nFor extensive flexible settings, put settings in parameters.ini.txt.\nOptions available on command line are as follows:\n");
    
    fprintf(stdout, "\n\t-V\tVerbose: print lots of human readable messages to\n\t\tstdout related to run status.\n");
}



//void viabilitySelection(void)
//{
//    long int i, j;
//    double fitnessValues[N];
//    
//    memset( fitnessValues, 0, (N * sizeof(double)) );
//    computeFitness(fitnessValues);
//    
//    
//    
//    printf("\nWarning: viabilitySelection() not written yet!\n");
//    exit(0);
//}


void wrongParametersIniOption(char *expected, char *previous, char *found)
{
    fprintf(stderr, "\nError in readInParametersFromFile():\n\t");
    fprintf(stderr, "Expecting '%s' after '%s'\n\t", expected, previous);
    fprintf(stderr, "in parameters.ini.txt, but instead found: '%s'.\n\tPlease fix parameters.ini.txt\n", found);
    fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
    exit(-1);
}






