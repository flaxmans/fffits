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
#include "fffits.h"

#include <gsl_rng.h>            // gnu scientific library //
#include <gsl_randist.h>        // gnu scientific library //
const gsl_rng_type *rngType;    /* generator type */
gsl_rng *rngState;              /* rng instance */



const char *version = "fffits1.0.0";

// globals for cmd line and or parameter file options
unsigned long long int nSITES = nSITES_DEFAULT;
long int nGENERATIONS = nGENERATIONS_DEFAULT;
int nPOPULATIONS = nPOPULATIONS_DEFAULT;
double MU = MU_DEFAULT, GENOME_MU;
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
long int TIME_SERIES_SAMPLE_FREQ = TIME_SERIES_SAMPLE_FREQ_DEFAULT;

short int *genotypes0, *genotypes1, *gts; // pointer to memory blocks for individual genotypes
int *locations, *linkageGroupMembership; // locations in discrete space of individuals
unsigned long long int *parentalTrackedSiteIndexes, *siteIndexes, *alleleCounts; // pointers for memory blocks for sites in genome
double *alleleFrequencies, *migrationRates, *K_VALUES, *M_VALUES, *selectionCoefficients;
short int *sitesStatuses; // codes for locus's current status (see below for codes)
short int currentBlock = 0; // keep track of which block is in use
unsigned long long int blockSizes[2];
long int t; // time counter
long int N; // current total population size
long int *abundances; // size of each population
unsigned long long int nTrackedSitesInParents = 0;
_Bool VERBOSE = 0;
int currentMigrationPeriod = 0;
int currentDemographyPeriod = 0;
double *migRatePt, *KvalPt;
FILE *dataFile_alleleFreqTS, *dataFile_alleleFreqTSbyPop, *dataFile_SFS_TS, *dataFile_segSiteTS;
FILE *dataFile_derivedFixationTS;
int *siteClassifications; // site classifications


// function declarations
long int calculateNumOffspring(int pop);
void calculatePopGenMetrics(void);
short int * checkMemoryBlocks(long int totalOffspring, long int nSitesInOffspring);
void chooseParents(long int *mommy, long int *daddy, double *dpt, long int nparents);
void chooseParentsAtRandom(long int *mommy, long int *daddy, long int *randomNumberLine, long int nhere);
void computeFitness(double *fitnessValues);
long int figureOutOffspringGenomeSites( unsigned long long int *offsp_SiteIndexes, short int *offsp_lociStates, long int nNewMutations, unsigned long long int *mutatedLoci );
void finalTasks(unsigned RNG_SEED);
void makeCumulativeFitnessNumLines(double *fitnessValues, double *fitnessNumLines, long int *individualsInDeme);
void makeDemesIndexes(long int *individualsInDeme);
void makeOneOffspring(long int momIndex, long int dadIndex, short int *offGTpt, long int nSitesInOffspring, unsigned long long int *offsp_SiteIndexes, short int *offsp_lociStates);
void migration(void);
void printParametersToFiles(unsigned RNG_SEED);
void putInMutations( short int *offspringGTs, short int *offsp_lociStates, unsigned long long int *offsp_SiteIndexes, long int nNewMutations, long int totalOffspring, long int nSitesInOffspring );
void reproduction(void);
//void viabilitySelection(void);



int main(int argc, char *argv[])
{
    unsigned int RNG_SEED;
	char *progname = argv[0];
	
	RNG_SEED = initializationSteps(argc, argv, progname);
    
    for ( t = 1; t <= nGENERATIONS; t++ ) {

        migration();
        
        reproduction(); // includes fecundity selection
        
        if ( N == 0 ) {
            // all extinct
            fprintf(stdout, "\nAll extinct at end of generation %li.\n", t);
            finalTasks(RNG_SEED);
            exit(1);
        }
        
        if ( (t % TIME_SERIES_SAMPLE_FREQ == 0) || (t == nGENERATIONS) )
            calculatePopGenMetrics();
    }
    
    finalTasks(RNG_SEED);
    return 0;
}


long int calculateNumOffspring(int pop)
{
    long int numOffspring = 0, focalTime, nhere;
    double expected, currentPop, k;
    
    // check to see if now is a time to change demographic parameters
    if ( nDEMOGRAPHIC_CHANGES > 0 ) {
        focalTime = DEMOGRAPHIC_CHANGE_TIMES[currentDemographyPeriod];
        if ( t == focalTime ) {
            KvalPt += nPOPULATIONS;
            currentDemographyPeriod++;
        }
    }
    
    nhere = *(abundances + pop);
    if ( nhere > 1 ) {
        if ( FIXED_POP_SIZE )
            numOffspring = (long int) *(KvalPt + pop);
        else {
            k = *(KvalPt + pop);
            currentPop = ((double) nhere);
            expected = currentPop + (currentPop * MAX_POP_GROWTH_RATE * (k - currentPop)/k); // logistic equation
            numOffspring = gsl_ran_poisson( rngState, expected );
        }
    }
    else
        numOffspring = 0;
    
    return numOffspring;
}


void calculatePopGenMetrics(void)
{
    long int i, j, alleleCountsByPopulation[(nPOPULATIONS * nTrackedSitesInParents)];
    unsigned long long int SFScountsByPopulation[(nPOPULATIONS * 2 * N)], nDivSites = 0, nPosSites = 0;
    unsigned long long int siteMasterIndex, nSegSites = 0, nNeutralSites = 0, nBGsites = 0, SFScounts[(PLOIDY * N)];
    int pop;
    short int *sipt, siteClass;
    long int nhere, count;
    double oneOver2N;
    
    oneOver2N = 1.0 / (((double) PLOIDY) * ((double) N));
    
    for ( count = 0; count < (PLOIDY * N); count++ )
        SFScounts[count] = 0;
    
    for ( i = 0; i < nTrackedSitesInParents; i++ ) {
        sipt = gts + (PLOIDY * i);
        
        // overall counts and frequencies
        siteMasterIndex = *(parentalTrackedSiteIndexes + i);
        if ( *(sitesStatuses + siteMasterIndex) == LOCUS_STATUS_VARIABLE_IN_PARENTS ) {
            // site types
            nSegSites++;
            siteClass = *(siteClassifications + siteMasterIndex);
            if ( siteClass == SITE_CLASS_NEUTRAL )
                nNeutralSites++;
            else if ( siteClass == SITE_CLASS_BGS )
                nBGsites++;
            else if ( siteClass == SITE_CLASS_DIV )
                nDivSites++;
            else if ( siteClass == SITE_CLASS_POS )
                nPosSites++;
            else {
                fprintf(stderr, "\nError in calculatePopGenMetrics():\n siteClass = %i not recognized\n", siteClass);
                exit(-1);
            }
            
            // global allele frequencies and SFS
            count = *(alleleCounts + siteMasterIndex);
            
            if ( count <= 0 ) {
                fprintf(stderr, "\nError in calculatePopGenMetrics():\n\tcount (%li) <= 0 for siteMasterIndex = %llu\n", count, siteMasterIndex);
                exit(-1);
            }
            fprintf(dataFile_alleleFreqTS, "%li,%llu,%i,%li,%E,%E,%i\n", t, siteMasterIndex, *(linkageGroupMembership + siteMasterIndex), count, (((double) count) * oneOver2N), *(selectionCoefficients + siteMasterIndex), *(siteClassifications + siteMasterIndex) );
            
            *(SFScounts + count) += 1;
        }
        
        for ( j = 0; j < nPOPULATIONS; j++ ) {
            
        }
    }
    
    // print segregating site counts
    fprintf(dataFile_segSiteTS, "%li,%llu,%llu,%llu,%llu,%llu\n", t, nSegSites, nNeutralSites, nBGsites, nPosSites, nDivSites);
    
    // print the SFS
    for ( i = 1; i < (PLOIDY * N); i++ ) {
        if ( *(SFScounts + i) ) {
            fprintf(dataFile_SFS_TS, "%li,%li,%llu\n", t, i, *(SFScounts + i));
        }
    }
    
    for ( i = 0; i < nPOPULATIONS; i++ ) {
        
    }
    
}



short int * checkMemoryBlocks(long int totalOffspring, long int nSitesInOffspring)
{
    int offspringBlock;
    unsigned long long int neededSize, makeSize;
    
    if ( currentBlock )
        offspringBlock = 0; // current parents use 1, so offspring will use 0
    else
        offspringBlock = 1; // current 0, offspring 1
    
    neededSize = PLOIDY * totalOffspring * nSitesInOffspring * (sizeof(short int));
    
    // check if we need bigger memory blocks and allocate if needed
    if ( blockSizes[offspringBlock] < neededSize ) {
        makeSize = PLOIDY * totalOffspring * (nSitesInOffspring + 10) * (sizeof(short int));
        if ( offspringBlock ) {
            free( genotypes1 );
            genotypes1 = (short int *) malloc( makeSize );
        }
        else {
            free( genotypes0 );
            genotypes0 = (short int *) malloc( makeSize );
        }
        blockSizes[offspringBlock] = makeSize;
    }
    
    if ( offspringBlock )
        return genotypes1;
    else
        return genotypes0;
}


void chooseParents(long int *mommy, long int *daddy, double *fnlpt, long int nparents)
{
    long int i, j;
    double dum, *dpt;

    dpt = fnlpt;
    // dpt is the pointer to first entry in the the fitness number line for this patch

    // choose mom
    dum = gsl_rng_uniform(rngState);
    i = 0;
    while ( dum > *dpt ) {
        i++;
        dpt++;
    }
    *mommy = i;
    
    // choose dad; make sure not the same as mom
    do {
        dum = gsl_rng_uniform(rngState);
        j = 0;
        dpt = fnlpt;
        while ( dum > *dpt ) {
            dpt++;
            j++;
        }
    } while ( j == i );
    *daddy = j;
    
//    if ( VERBOSE ) {
//        printf("\nMommy = %li, daddy = %li\n", *mommy, *daddy);
//    }

}


void chooseParentsAtRandom(long int *mommy, long int *daddy, long int *randomNumberLine, long int nhere)
{
    long int chosenOnes[2];
    
    gsl_ran_choose( rngState, chosenOnes, 2, randomNumberLine, nhere, sizeof(long int) );
    
    *mommy = chosenOnes[0];
    *daddy = chosenOnes[1];
}


void computeFitness(double *fitnessValues)
{
    unsigned long long int i, j, selectedSiteIndexesMaster[nSelectedSites];
    unsigned long long int selectedSiteIndexesLocal[nSelectedSites], *ullpt;
    double cf_scoeffs[nSelectedSites], *dpt, cf_s, cf_gradient; // prefix cf_ denoting local to this function
    long int nFound;
    int cf_siteClass, cf_siteClasses[nSelectedSites], *locpt;
    short int *sipt, gtsum;
//    unsigned long long int *dbpt1, *dbpt2; // pointers for debugging with lldb
//    double *dbpt3;
//    int *dbpt4;
    
    dpt = fitnessValues;
    for ( i = 0; i < N; i++ ) {
        *dpt = 1.0;
        dpt++;
    }
    
//    dbpt1 = selectedSiteIndexesMaster;
//    dbpt2 = selectedSiteIndexesLocal;
//    dbpt3 = cf_scoeffs;
//    dbpt4 = cf_siteClasses;
    
    nFound = 0;
    ullpt = parentalTrackedSiteIndexes;
    for ( i = 0; i < nTrackedSitesInParents; i++ ) {
        j = *ullpt; // focal site master index
        cf_siteClass = *(siteClassifications + j);
        if ( cf_siteClass != SITE_CLASS_NEUTRAL ) {
            cf_siteClasses[nFound] = cf_siteClass;
            selectedSiteIndexesLocal[nFound] = i;
            selectedSiteIndexesMaster[nFound] = j;
            cf_scoeffs[nFound] = *(selectionCoefficients + j);
            nFound++;
        }
        ullpt++;
    }
    
    if ( VERBOSE && t < 10 ) {
        printf("\nLocal and Master selected site indexes, site type codes, and selection coefficients:\n");
        for ( i = 0; i < nSelectedSites; i++ )
            printf("\t%llu\t%llu\t%i\t%f\n", selectedSiteIndexesLocal[i], selectedSiteIndexesMaster[i], cf_siteClasses[i], cf_scoeffs[i]);
    }
    
    for ( i = 0; i < nSelectedSites; i++ ) {
        sipt = gts + (PLOIDY * selectedSiteIndexesLocal[i]); // first individual, selected site i
        dpt = fitnessValues; // first individual's fitness score
        cf_siteClass = cf_siteClasses[i];
        cf_s = cf_scoeffs[i];
        locpt = locations;
        for ( j = 0; j < N; j++ ) {
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
            else if ( cf_siteClass == SITE_CLASS_NEUTRAL ) {
                fprintf(stderr, "\nError in computeFitness():\n\tneutral site made it into calcs.\n");
                exit(-1);
            }
            
            sipt += PLOIDY * nTrackedSitesInParents; // advance to locus i in next individual's genotype
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


long int figureOutOffspringGenomeSites( unsigned long long int *offsp_SiteIndexes, short int *offsp_lociStates, long int nNewMutations, unsigned long long int *mutatedLoci )
{
    long int nSitesInOffspring = 0, newMutationSiteCount = 0, parentalSiteCount = 0, nRepeatMutations = 0;
    unsigned long long int locus, *newLocusIndexPt, *variableSitePt, *offsp_pt;
    long int i, totalSitesDone = 0, retainedSites = 0, nSitesToDo, nBrandNewMutations = 0;
    short int *offsp_ls;

    newLocusIndexPt = mutatedLoci;
    variableSitePt = parentalTrackedSiteIndexes;
    offsp_pt = offsp_SiteIndexes;
    offsp_ls = offsp_lociStates;
    
    nSitesToDo = nNewMutations + nTrackedSitesInParents;
    while ( totalSitesDone < nSitesToDo ) {
        if ( (parentalSiteCount < nTrackedSitesInParents) && (newMutationSiteCount < nNewMutations) ) {
            // there are remaining sites to fill in for both types
            if ( *newLocusIndexPt > *variableSitePt ) {
                // add in the parental locus
                locus = *variableSitePt;
                if ( *(sitesStatuses + locus) == LOCUS_STATUS_VARIABLE_IN_PARENTS ) {
                    // should be added in, else ignored
                    *offsp_pt = locus; // record the locus
                    offsp_pt++; // advance offspring locus pointer
                    
                    *offsp_ls = LOCUS_STATUS_VARIABLE_IN_PARENTS; // record the status
                    offsp_ls++; // advance the status pointer
                    
                    nSitesInOffspring++; // increment the total offspring site counter
                    retainedSites++;
                }
                parentalSiteCount++; // one more parental site taken care of
                variableSitePt++; // increment the variable site pointer
                
                totalSitesDone++;
            }
            else if ( *newLocusIndexPt < *variableSitePt ) {
                // add in the new locus
                *offsp_pt = *newLocusIndexPt; // record the index
                offsp_pt++; // advance the index pointer
                
                *offsp_ls = LOCUS_STATUS_NEW_MUT_ONLY; // record the status
                offsp_ls++; // advance the status pointer
                
                newLocusIndexPt++; // advance new locus index pointer
                
                nSitesInOffspring++; // increment the total site counter
                
                newMutationSiteCount++;
                
                totalSitesDone++;
                
                nBrandNewMutations++;
            }
            else {
                // the mutation is at an already used site
                if ( *newLocusIndexPt != *variableSitePt ) {
                    fprintf(stderr, "\nError in figureOutOffspringGenomeSites():\n\t*newLocusIndexPt (%llu) != *variableSitePt (%llu)\n", *newLocusIndexPt, *variableSitePt);
                    exit(-1);
                }
                *offsp_pt = *newLocusIndexPt; // record the index
                offsp_pt++; // advance the index pointer
                
                *offsp_ls = LOCUS_STATUS_VARIABLE_PLUS_MUT; // record the status
                offsp_ls++; // advance the status pointer
                
                newLocusIndexPt++; // advance BOTH parental AND mutations pointers and their respective counts
                newMutationSiteCount++;
                variableSitePt++;
                parentalSiteCount++;
                
                nSitesInOffspring++; // increment count of sites in offspring
                nRepeatMutations++;
                
                
                totalSitesDone += 2; // took care of a parental site AND a mutation site at the same time
            }
        }
        else if ( parentalSiteCount < nTrackedSitesInParents && newMutationSiteCount == nNewMutations ) {
            // done with new, just need to add parentals
            locus = *variableSitePt;
            if ( *(sitesStatuses + locus) == LOCUS_STATUS_VARIABLE_IN_PARENTS ) {
                // should be added in, else ignored
                *offsp_pt = locus; // record the locus
                offsp_pt++; // advance offspring locus pointer
                
                *offsp_ls = LOCUS_STATUS_VARIABLE_IN_PARENTS; // record the status
                offsp_ls++; // advance the status pointer
                
                nSitesInOffspring++; // increment the total offspring site counter
                retainedSites++;
            }
            parentalSiteCount++; // one more parental site taken care of
            variableSitePt++; // increment the variable site pointer
            
            totalSitesDone++;
        }
        else if ( parentalSiteCount == nTrackedSitesInParents && newMutationSiteCount < nNewMutations ) {
            // done with parentals, just need to add new
            *offsp_pt = *newLocusIndexPt; // record the index
            offsp_pt++; // advance the index pointer
            
            *offsp_ls = LOCUS_STATUS_NEW_MUT_ONLY; // record the status
            offsp_ls++; // advance the status pointer
            
            newLocusIndexPt++; // advance new locus index pointer
            
            nSitesInOffspring++; // increment the total site counter
            
            newMutationSiteCount++;
            
            totalSitesDone++;
            
            nBrandNewMutations++;
        }
        else {
            fprintf(stderr, "\nError in figureOutOffspringGenomeSites():\n\tcounts aren't working how you think\n\t");
            fprintf(stderr, "parental count = %li\tnewMutationSiteCount = %li\n\tnTrackedSitesInParents = %llu, nNewMutations = %li\n", parentalSiteCount, newMutationSiteCount, nTrackedSitesInParents, nNewMutations);
            exit(-1);
        }
    }
    
    // error checking
    if ( nNewMutations != newMutationSiteCount ) {
        fprintf(stderr, "\nError in figureOutOffspringGenomeSites():\n\tnNewMutations (%li) != newMutationSiteCount (%li)\n", nNewMutations, newMutationSiteCount);
        exit(-1);
    }
    if ( newMutationSiteCount != (nRepeatMutations + nBrandNewMutations) ) {
        fprintf(stderr, "\nError in figureOutOffspringGenomeSites():\n\tnewMutationSiteCount (%li) != (nRepeatMutations (%li) + nBrandNewMutations (%li))\n", newMutationSiteCount, nRepeatMutations, nBrandNewMutations);
        exit(-1);
    }
    if ( totalSitesDone != (parentalSiteCount + newMutationSiteCount) ) {
        fprintf(stderr, "\nError in figureOutOffspringGenomeSites():\n\ttotalSitesDone (%li) != (parentalSiteCount (%li) + newMutationSiteCount (%li))\n", totalSitesDone, parentalSiteCount, newMutationSiteCount);
        exit(-1);
    }
    
    if ( VERBOSE && t < 10 ) {
        printf("\nParental loci and status codes:\n\tnumber\tindex\tcode\n");
        for ( i = 0; i < nTrackedSitesInParents; i++ ) {
            locus = *(parentalTrackedSiteIndexes + i);
            printf("\t%li\t%llu\t%i\n", i, locus, sitesStatuses[locus]);
        }
        printf("\nOffspring loci and status codes:\n\tnumber\tindex\tcode\n");
        for ( i = 0; i < nSitesInOffspring; i++ ) {
            locus = *(offsp_SiteIndexes + i);
            printf("\t%li\t%llu\t%i\n", i, locus, offsp_lociStates[i]);
        }
        for ( i = 1; i < nSitesInOffspring; i++ ) {
            if ( *(offsp_SiteIndexes + i) <= *(offsp_SiteIndexes + i - 1) ) {
                fprintf(stderr, "\nError in figureOutOffspringGenomeSites():\n\t*(offsp_SiteIndexes + i) (%llu) <= *(offsp_SiteIndexes + i - 1) (%llu)\n", *(offsp_SiteIndexes + i), *(offsp_SiteIndexes + i - 1));
                exit(-1);
            }
        }
        
        printf("\nnSitesInOffspring = %li\nNew sites from mutations = %li, retained sites = %li, repeat mutations = %li\n", nSitesInOffspring, newMutationSiteCount, retainedSites, nRepeatMutations);
    }
    
    return nSitesInOffspring;
}


void finalTasks(unsigned RNG_SEED)
{
    printParametersToFiles(RNG_SEED);
    
    fclose(dataFile_alleleFreqTS);
    fclose(dataFile_alleleFreqTSbyPop);
    fclose(dataFile_SFS_TS);
    fclose(dataFile_segSiteTS);
    fclose(dataFile_derivedFixationTS);
    
    free(genotypes0);
    free(genotypes1);
    free(parentalTrackedSiteIndexes);
    free(siteIndexes);
    free(sitesStatuses);
    free(alleleFrequencies);
    free(alleleCounts);
    free(siteClassifications);
    free(DEMOGRAPHIC_CHANGE_TIMES);
    free(K_VALUES);
    free(MIGRATION_CHANGE_TIMES);
    free(M_VALUES);
    free(abundances);
    free(selectionCoefficients);
    free(locations);
    free(environmentGradient);
    free(linkageGroupMembership);
}


void makeCumulativeFitnessNumLines(double *fitnessValues, double *fitnessNumLines, long int *individualsInDeme)
{
    int pop;
    long int i, j, nhere, *lipt, masterIndex, count;
    double fitSum, *dptfnl, fitval, *startpt;
    
    dptfnl = fitnessNumLines; // ordered for each deme
    lipt = individualsInDeme; // list of individual indexes in deme
    count = 0;
    for ( pop = 0; pop < nPOPULATIONS; pop++ ) { // go population by population
        nhere = abundances[pop]; // get total number here
        if ( nhere > 0 ) {
            startpt = dptfnl; // save starting point in array
            masterIndex = *lipt; // individual's index
            if ( VERBOSE ) {
                if ( pop != *(locations + masterIndex) ) {
                    fprintf(stderr, "\nError in makeCumulativeFitnessNumLines():\n\tlocation (%i) doesn't match pop (%i)!\n", *(locations + masterIndex), pop);
                    exit(-1);
                }
            }
            
            fitval = *(fitnessValues + masterIndex);  // individual's fitness
            fitSum = fitval; // first individual in deme: set cumulative sum start point
            *dptfnl = fitval; // save value in fitness number line
            lipt++; // advance individual indexing pointer
            dptfnl++; // advance number line pointer
            
            for ( i = 1; i < nhere; i++ ) { // start with 1 because 0 done above
                masterIndex = *lipt; // individual's index
                fitval = *(fitnessValues + masterIndex); // individual's fitness
                *dptfnl = *(dptfnl - 1) + fitval; // make number line cumulative;
                fitSum += fitval; // store total sum
                lipt++; // advance individual's index pointer
                dptfnl++; // advance fitness number line pointer
            }
            fitSum = 1.0 / fitSum; // now done with this deme's raw numbers; time to normalize
            for ( i = 0; i < nhere; i++ ) {
                *startpt *= fitSum;
                startpt++;
            }
        }
    }
    
    if ( VERBOSE ) {
        printf("\n");
        dptfnl = fitnessNumLines;
        for ( pop = 0; pop < nPOPULATIONS; pop++ ) {
            if ( abundances[pop] > 0 ) {
                printf("Deme %i, first fit = %f, last fit = %f\n", pop, *dptfnl, *(dptfnl + abundances[pop] - 1));
                dptfnl++;
                for ( i = 1; i < abundances[pop]; i++ ) {
                    if ( *dptfnl < *(dptfnl - 1) ) {
                        printf("\nError in makeCumulativeFitnessNumLines():\n\tDecreasing values in fitness vector:\n\t%f\t%f\n", *dptfnl, *(dptfnl - 1) );
                        exit(-1);
                    }
                    dptfnl++;
                }
            }
            else
                printf("No individuals in Deme %i\n", pop);
        }
    }
    
}


void makeDemesIndexes(long int *individualsInDeme)
{
    // this function makes a vector of consecutive indexes in each deme
    // this is done because -- at times that this is called -- the individuals are unsorted due to migration
    // with high migration rates, maybe a sort would actually be more efficient???
    
    long int i, j, demePosition[nPOPULATIONS], spot;
    int *ipt, loc;
    
    demePosition[0] = 0;
    for ( i = 1; i < nPOPULATIONS; i++ )
        demePosition[i] = demePosition[(i-1)] + abundances[(i-1)];
    
    
    ipt = locations;
    
    for ( i = 0; i < N; i++ ) {
        loc = *ipt;
        spot = demePosition[loc];
        *(individualsInDeme + spot) = i;
        demePosition[loc] = demePosition[loc] + 1;
        
        ipt++;
    }
    
    if ( VERBOSE ) {
        FILE *testMakeIndexes;
        testMakeIndexes = fopen("TestMakeDemesIndexes.txt","w");
        for ( i = 0; i < N; i++ ) {
            j = individualsInDeme[i];
            fprintf(testMakeIndexes, "%li,%i\n", j, locations[j]);
        }
        fclose(testMakeIndexes);
    }
}


void makeOneOffspring(long int momIndex, long int dadIndex, short int *offGTpt, long int nSitesInOffspring, unsigned long long int *offsp_SiteIndexes, short int *offsp_lociStates)
{
    int i, currentLinkageGroup;
    unsigned long long int j, focalSite, *offsp_SIpt, *parentalLocusIndexes;
    short int *sipt, *parentPoint, *offsp_ls;
    int chromosome;
    double meanRecombDistance = 1000.0 / RECOMBINATION_RATE_PER_KB;
    unsigned long long int nextRecombinationSpot;
    long int parentalLocusCounter;
    
    if ( VERBOSE ) {
        if ( *(locations + momIndex) != *(locations + dadIndex) ) {
            printf("\nError in makeOneOffspring:\n\tIndexes are off giving bogus mom or dad\n");
            exit(-1);
        }
    }
    
    for ( i = 0; i < 2; i++ ) {
        
        sipt = offGTpt + i; // which haploid set in offspring
        
        
        parentalLocusCounter = 0;
        parentalLocusIndexes = parentalTrackedSiteIndexes;
        if ( nTrackedSitesInParents ) {
            while ( *(sitesStatuses + *parentalLocusIndexes) != LOCUS_STATUS_VARIABLE_IN_PARENTS ) {
                parentalLocusCounter++;
                parentalLocusIndexes++;
            }
        }
        
        // which parent
        if ( i == 0 )
            parentPoint = gts + (PLOIDY * nTrackedSitesInParents * momIndex ) + (PLOIDY * parentalLocusCounter);
        else
            parentPoint = gts + (PLOIDY * nTrackedSitesInParents * dadIndex ) + (PLOIDY * parentalLocusCounter);
        
        
        
        // which chromosome to start with in the parent
        if ( gsl_rng_uniform( rngState ) < 0.5 ) {
            chromosome = 0;
        }
        else {
            chromosome = 1;
            parentPoint++; // advance 1 to "second" choromosome
        }
        currentLinkageGroup = 0;

        offsp_ls = offsp_lociStates;
        offsp_SIpt = offsp_SiteIndexes;
        nextRecombinationSpot = (unsigned long long int) randExp( meanRecombDistance );
        for ( j = 0; j < nSitesInOffspring; j++ ) {
            focalSite = *offsp_SIpt;
            
            // now copy alleles from parents to offspring
            if ( *offsp_ls == LOCUS_STATUS_VARIABLE_IN_PARENTS || *offsp_ls == LOCUS_STATUS_VARIABLE_PLUS_MUT ) {
                
                // first handle recombination and independent assortment
                if ( *(linkageGroupMembership + focalSite) != currentLinkageGroup ) {
                    // new "chromosome"; need to implement independent assortment
                    if ( gsl_rng_uniform(rngState) < 0.5 ) {
                        // flip to other chromosome in parent
                        if ( chromosome ) {
                            parentPoint--; // move back one
                            chromosome = 0;
                        }
                        else {
                            parentPoint++; // move forward one
                            chromosome = 1;
                        }
                    }
                    if ( *(linkageGroupMembership + focalSite) < currentLinkageGroup ) {
                        fprintf(stderr, "\nError in makeOneOffspring():\n\t*(linkageGroupMembership + focalSite) (%i) < currentLinkageGroup (%i)\n", *(linkageGroupMembership + focalSite), currentLinkageGroup );
                        exit(-1);
                    }
                    currentLinkageGroup = *(linkageGroupMembership + focalSite);
                    
                    nextRecombinationSpot = focalSite + ((unsigned long long int) randExp( meanRecombDistance ));
                }
                else if ( focalSite > nextRecombinationSpot ) {
                    // implement recombination within a linkage group
                    do {
                        // do-while loop allows for multiple recombination events
                        if ( chromosome ) {
                            parentPoint--; // move back one
                            chromosome = 0;
                        }
                        else {
                            parentPoint++; // move forward one
                            chromosome = 1;
                        }
                        nextRecombinationSpot += ((unsigned long long int) randExp( meanRecombDistance ));
                    } while ( focalSite > nextRecombinationSpot );
                }
                
                *sipt = *parentPoint;
                parentalLocusCounter++;  // one more parental locus taken care of
                if ( j < (nSitesInOffspring - 1) && parentalLocusCounter < nTrackedSitesInParents ) {
                    parentalLocusIndexes++; // advance to index of next parental site/locus
                    parentPoint += PLOIDY; // advance to next site in genome
                    while ( *(sitesStatuses + *parentalLocusIndexes) != LOCUS_STATUS_VARIABLE_IN_PARENTS && parentalLocusCounter < nTrackedSitesInParents ) {
                        parentalLocusCounter++; // count prior locus toward those that are "done"
                        // the while part helps skip over loci that will no longer be tracked
                        if ( parentalLocusCounter < nTrackedSitesInParents ) {
                            // if statement prevents seg fault on last While check
                            parentalLocusIndexes++; // advance to index of next parental site/locus
                            parentPoint += PLOIDY; // advance to next site in genome
                        }
                    }
                }
            }
            else if ( *offsp_ls == LOCUS_STATUS_NEW_MUT_ONLY ) {
                *sipt = ALLELE_CODE_ANCESTRAL; // just put in the locus as a placeholder; mutations added later
                // do NOT advance any parent pointers or counters
            }
            else {
                fprintf(stderr, "\nError in makeOneOffspring():\n\t*offsp_ls = %i code not recognized.\n", *offsp_ls);
                exit(-1);
            }
            
            if ( *sipt )
                *(alleleCounts + focalSite) += *sipt;
            
            sipt += PLOIDY;
            offsp_SIpt++;
            offsp_ls++;
        }
    }
    
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
    // make an R-sourceable script of parameter values
    int i, j;
    FILE *rfile;
    double *dpt, THETA;
    unsigned long long int foo;
    
    rfile = fopen("metadataAndParameters.R","w");
    fprintf(rfile, "RNG_SEED <- %u\n", RNG_SEED);
    fprintf(rfile, "TIME_SERIES_SAMPLE_FREQ <- %li\n", TIME_SERIES_SAMPLE_FREQ);
    fprintf(rfile, "nGENERATIONS <- %li\n", nGENERATIONS);
    fprintf(rfile, "nPOPULATIONS <- %i\n", nPOPULATIONS);
    fprintf(rfile, "\n");
    
    // demography
    fprintf(rfile, "# demography\nnDEMOGRAPHIC_CHANGES <- %i\n", nDEMOGRAPHIC_CHANGES);
    if ( nDEMOGRAPHIC_CHANGES > 0 ) {
        if ( nDEMOGRAPHIC_CHANGES == 1 )
            fprintf(rfile, "DEMOGRAPHIC_CHANGE_TIMES <- %li\n", *DEMOGRAPHIC_CHANGE_TIMES);
        else {
            fprintf(rfile, "DEMOGRAPHIC_CHANGE_TIMES <- c(%li", *DEMOGRAPHIC_CHANGE_TIMES);
            for ( i = 1; i < nDEMOGRAPHIC_CHANGES; i++ )
                fprintf(rfile, ",%li", *(DEMOGRAPHIC_CHANGE_TIMES + i));
            fprintf(rfile, ")\n");
        }
    }
    else
        fprintf(rfile, "DEMOGRAPHIC_CHANGE_TIMES <- -1\n");
    
    fprintf(rfile, "K_VALUES <- c(%E", *K_VALUES);
    for ( i = 1; i < (nPOPULATIONS * (nDEMOGRAPHIC_CHANGES+1)); i++ ) {
        fprintf(rfile, ",%E", *(K_VALUES+i));
    }
    fprintf(rfile, ")\n");
    fprintf(rfile, "K_VALUES <- matrix(data = K_VALUES, nrow = nDEMOGRAPHIC_CHANGES+1, ncol = nPOPULATIONS, byrow = TRUE)\n");
    fprintf(rfile, "FIXED_POP_SIZE <- %i\n", FIXED_POP_SIZE);
    fprintf(rfile, "MAX_POP_GROWTH_RATE <- %E\n", MAX_POP_GROWTH_RATE);
    fprintf(rfile, "\n");
    
    // migration parameters
    fprintf(rfile, "# migration parameters\nnMIGRATION_CHANGES <- %i\n", nMIGRATION_CHANGES);
    if ( nMIGRATION_CHANGES > 0 ) {
        if ( nMIGRATION_CHANGES == 1 )
            fprintf(rfile, "MIGRATION_CHANGE_TIMES <- %li\n", *MIGRATION_CHANGE_TIMES);
        else {
            fprintf(rfile, "MIGRATION_CHANGE_TIMES <- c(%li", *MIGRATION_CHANGE_TIMES);
            for ( i = 1; i < nMIGRATION_CHANGES; i++ )
                fprintf(rfile, ",%li", *(MIGRATION_CHANGE_TIMES + i));
            fprintf(rfile, ")\n");
        }
    }
    else
        fprintf(rfile, "MIGRATION_CHANGE_TIMES <- -1\n");
    dpt = M_VALUES;
    for ( i = 1; i <= (nMIGRATION_CHANGES + 1); i++ ) {
        fprintf(rfile, "M_VALUES_period_%i <- c(%E", i, *dpt);
        dpt++;
        for ( j = 1; j < (nPOPULATIONS * nPOPULATIONS); j++ ) {
            fprintf(rfile, ",%E", *dpt);
            dpt++;
        }
        fprintf(rfile, ")\n");
        fprintf(rfile, "M_VALUES_period_%i <- matrix(data = M_VALUES_period_%i, nrow = nPOPULATIONS, ncol = nPOPULATIONS, byrow = TRUE)\n", i, i);
    }
    fprintf(rfile, "\n");
    
    // genome
    fprintf(rfile, "# genome parameters\nnSITES <- %llu\n", nSITES);
    fprintf(rfile, "nLINKAGE_GROUPS <- %i\n", nLINKAGE_GROUPS);
    fprintf(rfile, "PLOIDY <- %i\n", PLOIDY);
    fprintf(rfile, "RECOMBINATION_RATE_PER_KB <- %E\n", RECOMBINATION_RATE_PER_KB);
    fprintf(rfile, "MU <- %E\n", MU);
    fprintf(rfile, "GENOME_MU <- %E\n", GENOME_MU);
    foo = N * PLOIDY * nSITES;
    THETA = 2.0 * MU * ((double) foo);
    fprintf(rfile, "THETA <- %E\n", THETA);
    fprintf(rfile, "\n");
    
    // natural selection
    fprintf(rfile, "# natural selection\nINCLUDE_SELECTION <- %i\n", INCLUDE_SELECTION);
    fprintf(rfile, "MEAN_S_BGS <- %E\n", MEAN_S_BGS);
    fprintf(rfile, "MEAN_S_POS <- %E\n", MEAN_S_POS);
    fprintf(rfile, "MEAN_S_DIV <- %E\n", MEAN_S_DIV);
    fprintf(rfile, "PROBABILITY_SITE_BGS <- %E\n", PROBABILITY_SITE_BGS);
    fprintf(rfile, "PROBABILITY_SITE_POS <- %E\n", PROBABILITY_SITE_POS);
    fprintf(rfile, "PROBABILITY_SITE_DIV <- %E\n", PROBABILITY_SITE_DIV);
    fprintf(rfile, "ENVIRONMENT_TYPE <- %i\n", ENVIRONMENT_TYPE);
    fprintf(rfile, "ENVT_MIN <- %E\n", ENVT_MIN);
    fprintf(rfile, "ENVT_MAX <- %E\n", ENVT_MAX);
    fprintf(rfile, "environmentGradient <- c(%E", *environmentGradient);
    for ( i = 1; i < nPOPULATIONS; i++ )
        fprintf(rfile, ",%E", *(environmentGradient + i));
    fprintf(rfile, ")\n");
    fprintf(rfile, "FITNESS_MODEL <- %i\n", FITNESS_MODEL);
    fprintf(rfile, "\n");
    
    // codes and magic numbers
    fprintf(rfile, "# codes and magic numbers\nnSITE_CLASSES <- %i\n", nSITE_CLASSES);
    fprintf(rfile, "SITE_CLASS_NEUTRAL <- %i\n", SITE_CLASS_NEUTRAL);
    fprintf(rfile, "SITE_CLASS_BGS <- %i\n", SITE_CLASS_BGS);
    fprintf(rfile, "SITE_CLASS_POS <- %i\n", SITE_CLASS_POS);
    fprintf(rfile, "SITE_CLASS_DIV <- %i\n", SITE_CLASS_DIV);
    fprintf(rfile, "ALLELE_CODE_ANCESTRAL <- %i\n", ALLELE_CODE_ANCESTRAL);
    fprintf(rfile, "ALLELE_CODE_DERIVED <- %i\n", ALLELE_CODE_DERIVED);
    fprintf(rfile, "ENVT_TYPE_GRADIENT <- %i\n", ENVT_TYPE_GRADIENT);
    fprintf(rfile, "ENVT_TYPE_MOSAIC <- %i\n", ENVT_TYPE_MOSAIC);
    fprintf(rfile, "ENVT_TYPE_INVARIANT <- %i\n", ENVT_TYPE_INVARIANT);
    
    fprintf(rfile, "FITNESS_MODEL_ADDITIVE <- %i\n", FITNESS_MODEL_ADDITIVE);
    fprintf(rfile, "FITNESS_MODEL_MULTIPLICATIVE <- %i\n", FITNESS_MODEL_MULTIPLICATIVE);
    fprintf(rfile, "LOCUS_STATUS_INACTIVE <- %i\n", LOCUS_STATUS_INACTIVE);
    fprintf(rfile, "LOCUS_STATUS_VARIABLE_IN_PARENTS <- %i\n", LOCUS_STATUS_VARIABLE_IN_PARENTS);
    fprintf(rfile, "LOCUS_STATUS_VARIABLE_PLUS_MUT <- %i\n", LOCUS_STATUS_VARIABLE_PLUS_MUT);
    fprintf(rfile, "LOCUS_STATUS_NEW_MUT_ONLY <- %i\n", LOCUS_STATUS_NEW_MUT_ONLY);
    fprintf(rfile, "LOCUS_STATUS_TRACKED_IN_PARENTS <- %i\n", LOCUS_STATUS_TRACKED_IN_PARENTS);
    fprintf(rfile, "CODOMINANCE <- %E\n", CODOMINANCE);
    fprintf(rfile, "\n");
    
    // states at end of run
    fprintf(rfile, "#states of variables at the end of the run\nN <- %li\n", N);
    fprintf(rfile, "nSelectedSites <- %li\n", nSelectedSites);
    fprintf(rfile, "nTrackedSitesInParents <- %llu\n", nTrackedSitesInParents);
    fprintf(rfile, "abundances <- c(%li", abundances[0]);
    for ( i = 1; i < nPOPULATIONS; i++ )
        fprintf(rfile, ",%li", abundances[i]);
    fprintf(rfile, ")\n");
    
    fprintf(rfile, "\n");
    
    fclose(rfile);
}


void putInMutations( short int *offspringGTs, short int *offsp_lociStates, unsigned long long int *offsp_SiteIndexes, long int nNewMutations, long int totalOffspring, long int nSitesInOffspring )
{
    long int i, individual, totalMutsAdded = 0;
    short int *spot, *offsp_ls;
    unsigned long long int locus;
    
    // only considering a possibility of two different alleles at each site

    offsp_ls = offsp_lociStates;
    for ( i = 0; i < nSitesInOffspring; i++ ) {
        if ( *offsp_ls == LOCUS_STATUS_NEW_MUT_ONLY || *offsp_ls == LOCUS_STATUS_VARIABLE_PLUS_MUT ) {
            individual = gsl_rng_uniform_int( rngState, totalOffspring );
            //unsigned long int gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n)
            //this function returns a random integer from 0 to n-1
            spot = offspringGTs + ( PLOIDY * nSitesInOffspring * individual ) + ( i * PLOIDY );
            // spot = site in genome array
            
            // choose one of the two alleles at random
            if ( gsl_rng_uniform( rngState ) < 0.5 )
                spot++;
            
            locus = *(offsp_SiteIndexes + i);
            
            // if it is a derived allele, mutate it back to the ancestral state
            if ( *spot == ALLELE_CODE_DERIVED ) {
                *spot = ALLELE_CODE_ANCESTRAL;
                *(alleleCounts + locus) -= 1;
            }
            else {
                // otherwise, insert a derived allele
                *spot = ALLELE_CODE_DERIVED;
                *(alleleCounts + locus) += 1;
            }
            
            // double check allele counts
            if ( *(alleleCounts + locus) > (PLOIDY * totalOffspring) ) {
                fprintf(stderr, "\nError in putInMutations():\n\t*(alleleCounts + locus) (%llu) out of bounds\n", *(alleleCounts + locus));
                exit(-1);
            }
            
            totalMutsAdded++;
        }
        offsp_ls++;
    }
    
    if ( totalMutsAdded != nNewMutations ) {
        fprintf( stderr, "\nError in putInMutations():\n\ttotalMutsAdded (%li) != nNewMutations (%li)\n", totalMutsAdded, nNewMutations);
        exit(-1);
    }
    
}


double randExp(double meanValue)
{
    return ( (log(1.0 - gsl_rng_uniform(rngState))) * (-meanValue) );
}


unsigned readInParametersFromFile(void)
{
    char c, option[80], option2[80];
    long int INITIAL_N;
    unsigned long long int foo;
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
        else if ( !strcmp( option, "MU"  ) ) { // background selection mean coefficient
            fscanf(pfile, "%lf", &MU);
            if ( VERBOSE )
                printf("Found MU (%s) = %E\n", option, MU);
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
        else if ( !strcmp( option, "TIME_SERIES_SAMPLE_FREQ"  ) ) { // divergent selection mean coefficient
            fscanf(pfile, "%li", &TIME_SERIES_SAMPLE_FREQ);
            if ( VERBOSE )
                printf("Found TIME_SERIES_SAMPLE_FREQ (%s) = %li\n", option, TIME_SERIES_SAMPLE_FREQ);
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
    GENOME_MU = MU * ((double) PLOIDY) * ((double) nSITES);
    
    return RNG_SEED;
}



void reproduction(void)
{
    int pop, *locpt, *offspringLocations, nDerivedFixations = 0;
    long int i, j, noffspring[nPOPULATIONS], totalOffspring, mommy, daddy, noff, nhere;
    long int momIndex, dadIndex, driftLosses, actuallyVariable;
    short int *offspringGTs, *sipt;
    long int individualsInDeme[N], nNewMutations, *lipt, randomNumberLine[N], nSitesInOffspring, maxSitesInOffspring;
    double mutRate;
    double fitnessValues[N], fitnessNumLines[N], *dpt;
    unsigned long long int locus, *ullipt;
    unsigned long long int *offsp_SiteIndexes;
    
    if ( INCLUDE_SELECTION )
        computeFitness(fitnessValues);
    else {
        for ( i = 0; i < N; i++ ) {
            randomNumberLine[i] = i;
        }
    }
    
    totalOffspring = 0;
    for ( pop = 0; pop < nPOPULATIONS; pop++ ) {
        noffspring[pop] = calculateNumOffspring(pop);
        totalOffspring += noffspring[pop];
    }

    offspringLocations = (int *) malloc( totalOffspring * sizeof(int) );
    
    if ( VERBOSE ) {
        printf("\nNumber of offspring in each deme:\n");
        for ( i = 0; i < nPOPULATIONS; i++ )
            printf("\t%li", noffspring[i]);
        printf("\n");
    }

    // figure out how many mutations and which sites, including possibility for repeat and back mutation
    nNewMutations = gsl_ran_poisson( rngState, (GENOME_MU * ((double) totalOffspring)) );
    unsigned long long int mutatedLoci[nNewMutations];
    gsl_ran_choose( rngState, mutatedLoci, nNewMutations, siteIndexes, nSITES, sizeof(unsigned long long int) );
    maxSitesInOffspring = nTrackedSitesInParents + nNewMutations;
    short int offsp_lociStates[maxSitesInOffspring];
    offsp_SiteIndexes = (unsigned long long int *) malloc( maxSitesInOffspring * sizeof( unsigned long long int ) );
    
    // determine which sites will need to be tracked among offspring
    nSitesInOffspring = figureOutOffspringGenomeSites( offsp_SiteIndexes, offsp_lociStates, nNewMutations, mutatedLoci );
    
    // debug
    //printf("\nreproduction: time, nTrackedSitesInParents, nSitesInOffspring, nNewMutations, N\n\t%li\t%llu\t%li\t%li\t%li", t, nTrackedSitesInParents, nSitesInOffspring, nNewMutations, N);
    
    // make convenient arrays for choosing parents
    makeDemesIndexes(individualsInDeme); // this is instead of doing a sort; maybe doing a sort will turn out to be better???
    if ( INCLUDE_SELECTION )
        makeCumulativeFitnessNumLines( fitnessValues, fitnessNumLines, individualsInDeme );
    
    // now, get a memory block to use to store the offspring genotypes
    offspringGTs = checkMemoryBlocks(totalOffspring, nSitesInOffspring);
    // choose parents and make offspring
    sipt = offspringGTs; // pointer to beginning of offspring genotype array
    dpt = fitnessNumLines; // pointer to first number line
    locpt = offspringLocations; // pointer to location of first individual
    lipt = individualsInDeme;
    memset( alleleCounts, 0, nSITES * sizeof(unsigned long long int) );
    
    // some debug code
    /*
    if ( t >= 10 ) {
        printf("\nTracked sites in parents (%llu) @ t = %li:\n", nTrackedSitesInParents, t);
        for ( i = 0; i < nTrackedSitesInParents; i++ )
            printf("%llu\t", parentalTrackedSiteIndexes[i]);
        
        printf("\n\nTSites in offspring (%li):\n", nSitesInOffspring);
        for ( i = 0; i < nSitesInOffspring; i++ )
            printf("%llu\t", offsp_SiteIndexes[i]);
        printf("\n\n");
    }
     */
    
    
    
    for ( pop = 0; pop < nPOPULATIONS; pop++ ) {
        noff = noffspring[pop];
        nhere = abundances[pop];
        if ( noff > 0 ) {
            for ( i = 0; i < noff; i++ ) {
                if ( INCLUDE_SELECTION )
                    chooseParents(&mommy, &daddy, dpt, nhere);
                else
                    chooseParentsAtRandom(&mommy, &daddy, randomNumberLine, nhere);
                momIndex = *(lipt + mommy);
                dadIndex = *(lipt + daddy);
                makeOneOffspring(momIndex, dadIndex, sipt, nSitesInOffspring, offsp_SiteIndexes, offsp_lociStates);
                *locpt = pop;
                locpt++;
                sipt += (PLOIDY * nSitesInOffspring);
            }
        }
        if ( nhere > 0 ) {
            dpt += nhere;
            lipt += nhere;
        }
    }
    
    putInMutations( offspringGTs, offsp_lociStates, offsp_SiteIndexes, nNewMutations, totalOffspring, nSitesInOffspring );
    
    // switch the genotype pointer to the offspring's
    gts = offspringGTs;
    if ( currentBlock )
        currentBlock = 0;
    else
        currentBlock = 1;
    for ( i = 0; i < nPOPULATIONS; i++ ) {
        abundances[i] = noffspring[i];
    }
    N = totalOffspring;
    
    ullipt = offsp_SiteIndexes;
    memset( sitesStatuses, LOCUS_STATUS_INACTIVE, nSITES * sizeof(short int) );
    nSelectedSites = 0;
    for ( i = 0; i < nSitesInOffspring; i++ ) {
        locus = *ullipt;
        *(sitesStatuses + locus) = LOCUS_STATUS_TRACKED_IN_PARENTS;
        if ( *(siteClassifications + locus) > SITE_CLASS_NEUTRAL )
            nSelectedSites++;
        ullipt++;
    }
    
    ullipt = alleleCounts;
    sipt = sitesStatuses;
    driftLosses = 0;
    actuallyVariable = 0;
    for ( locus = 0; locus < nSITES; locus++ ) {
        if ( *ullipt > 0 && *ullipt < ( PLOIDY * N ) ) {
            *sipt = LOCUS_STATUS_VARIABLE_IN_PARENTS;
            actuallyVariable++;
        }
        else if ( *ullipt == (PLOIDY * N) ) {
            fprintf( dataFile_derivedFixationTS, "%li,%llu,%i,%E\n", t, locus, *(siteClassifications + locus), *(selectionCoefficients + locus) );
            nDerivedFixations++;
        }
        else if ( *sipt == LOCUS_STATUS_TRACKED_IN_PARENTS ) {
            if ( VERBOSE )
                printf("\nLocus %llu lost in time step %li\n", locus, t);
            driftLosses++;
        }
        sipt++;
        ullipt++;
    }
    
    if ( (driftLosses + actuallyVariable + nDerivedFixations) != nSitesInOffspring ) {
        fprintf(stderr, "\nError in reproduction():\n\t(driftLosses (%li) + actuallyVariable (%li)) + nDerivedFixations (%i) != nSitesInOffspring (%li) )\n", driftLosses, actuallyVariable, nDerivedFixations, nSitesInOffspring);
        exit(-1);
    }
    
    nTrackedSitesInParents = nSitesInOffspring;
    
    free(locations);
    locations = offspringLocations;
    
    free(parentalTrackedSiteIndexes);
    parentalTrackedSiteIndexes = offsp_SiteIndexes;
    
    // need to switch over/record nTrackedSitesInParents and variable site indexes
    // memset all isVariableSite to zero
    // then reset those variable ones to 1
    // also copy over offspring locations
    
    //printf("\nWarning: reproduction() not finished yet!\n");
    
    //exit(0);
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







