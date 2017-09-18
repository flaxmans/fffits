/* fast forward-in-time simulator for prediction of population genetic metrics and 
 hopefully some ABC inference, but maybe that's a pipe dream? */

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
unsigned long int nSITES = nSITES_DEFAULT;
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
int ENVIRONMENT_TYPE = ENVT_TYPE_GRADIENT;
double ENVT_MAX = ENVT_MAX_DEFAULT;
double ENVT_MIN = ENVT_MIN_DEFAULT;
int FITNESS_MODEL = FITNESS_MODEL_MULTIPLICATIVE; // defaults for how selection works; see defines in fffits.h
long int TIME_SERIES_SAMPLE_FREQ = TIME_SERIES_SAMPLE_FREQ_DEFAULT;
_Bool VERBOSE = 0;

long int nSelectedSites;
double *environmentGradient;
short int *genotypes0, *genotypes1, *gts; // pointer to memory blocks for individual genotypes
int *locations, *linkageGroupMembership; // locations in discrete space of individuals
unsigned long int *parentalTrackedSiteIndexes, *alleleCounts; // pointers for memory blocks for sites in genome
double *alleleFrequencies, *migrationRates, *K_VALUES, *M_VALUES, *selectionCoefficients;
short int *sitesStatuses; // codes for locus's current status (see below for codes)
short int currentBlock = 0; // keep track of which block is in use
unsigned long int blockSizes[2];
long int t; // time counter
long int N; // current total population size
long int *abundances; // size of each population
unsigned long int nTrackedSitesInParents = 0;
int currentMigrationPeriod = 0;
int currentDemographyPeriod = 0;
double *migRatePt, *KvalPt;
FILE *dataFile_alleleFreqTS, *dataFile_SFS_TS, *dataFile_segSiteTS; // *dataFile_alleleFreqTSbyPop,
FILE *dataFile_derivedFixationTS, *dataFile_abundances, *dataFile_PiAndDXY;
int *siteClassifications; // site classifications
unsigned long int totalMutationsForRun = 0;
unsigned long int expIncrement, halfIncr;
double *expLookupTable,expMedianVal;
unsigned long int expLookupMedian;


// function declarations
long int arraySearch(double *fnlpt, long int nparents);
long int calculateNumOffspring(int pop);
short int * checkMemoryBlocks(long int totalOffspring, long int nSitesInOffspring);
void chooseParentsAtRandom(long int *mommy, long int *daddy, long int *randomNumberLine, long int nhere);
void computeFitness(double *fitnessValues);
long int figureOutOffspringGenomeSites( unsigned long int *offsp_SiteIndexes, short int *offsp_lociStates, long int nNewMutations, unsigned long int *mutatedLoci, _Bool *copyFromParents );
void figureOutParentalSkippedSites( _Bool *skipTheseParentalSites );
void finalTasks(unsigned RNG_SEED);
void makeCumulativeFitnessNumLines(double *fitnessValues, double *fitnessNumLines, long int *individualsInDeme);
void makeDemesIndexes(long int *individualsInDeme);
//void makeOneOffspring(long int momIndex, long int dadIndex, short int *offGTpt, long int nSitesInOffspring, unsigned long int *offsp_SiteIndexes, short int *offsp_lociStates);
void makeOneOffspring(long int momIndex, long int dadIndex, short int *offGTpt, long int nSitesInOffspring, unsigned long int *offsp_SiteIndexes, _Bool *copyFromParents, _Bool *skipTheseParentalSites, unsigned long int *focalSites, int *siteLGs);
void migration(void);
void myInsertSort(unsigned long int * array, long int numElements);
void printParametersToFiles(unsigned RNG_SEED);
void putInMutations( short int *offspringGTs, short int *offsp_lociStates, unsigned long int *offsp_SiteIndexes, long int nNewMutations, long int totalOffspring, long int nSitesInOffspring );
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
            dataRecording();
    }
    
    finalTasks(RNG_SEED);
    return 0;
}


long int arraySearch(double *fnlpt, long int nparents)
{
	long int i;
	double dum;
	
	dum = gsl_rng_uniform( rngState );
	
	i = (long int) (dum * ((double) nparents)); // estimate of starting point
	// since i is truncated, it should fall between 0 and (nparents - 1), inclusive
	fnlpt = fnlpt + i;
	if ( dum > *fnlpt ) {
		// must go "up"
		while ( dum > *fnlpt && i < (nparents-1) ) {
			i++;
			fnlpt++;
		}
	}
	else if ( i > 0 ) {
		// may need to go down
		if ( dum < *(fnlpt - 1) ) {
			// must go down
			while ( dum < *(fnlpt-1) && i > 0 ) {
				i--;
				fnlpt--;
			}
		}
	}
	
#ifdef DEBUG
	if ( i < 0 || i >= nparents) {
		fprintf(stderr, "\nError in arraySearch():\n\ti = %li out of bounds!\n", i);
		exit(-1);
	}
	
#endif
	
	return i;
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
		// ensures that selfing is not allowed because later on
		// mom and dad must be different individuals (in parent choice algorithm)
        if ( FIXED_POP_SIZE )
            numOffspring = (long int) *(KvalPt + pop);
        else {
            k = *(KvalPt + pop);
            currentPop = ((double) nhere);
            expected = currentPop + (currentPop * MAX_POP_GROWTH_RATE * (k - currentPop)/k); // logistic equation
            numOffspring = gsl_ran_poisson( rngState, expected ); // potential optimization with lookup table??
        }
    }
    else
        numOffspring = 0;
    
    return numOffspring;
}



short int * checkMemoryBlocks(long int totalOffspring, long int nSitesInOffspring)
{
    int offspringBlock;
    unsigned long int neededSize, makeSize;
    
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


void chooseMutationSites( unsigned long int *mutatedLoci, long int nNewMutations )
{
	unsigned long int testSite, *ulipt = mutatedLoci;
	// want random integers between 0 and (nSITES - 1)
	long int i, counter = 1;
	_Bool looking;
	
	// pick first site; no need to worry about replacement yet because it's the first one
	testSite = gsl_rng_uniform_int(rngState, nSITES);
	*ulipt = testSite;
	ulipt++; // increment pointer to NEXT spot in destination array
	for ( counter = 1; counter < nNewMutations; counter++ ) {
		// make sure to sample without replacement:
		do {
			looking = 0;
			// get a candidate site:
			testSite = gsl_rng_uniform_int(rngState, nSITES);
			// check all previous sites to see if the candidate was used already:
			// VERY low probability, but would screw up other algorithms in
			// reproduction() if it occurred
			for ( i = 0; i < counter; i++ ) {
				if ( testSite == *(mutatedLoci + i) )
					looking = 1;
			}
		} while ( looking );
		*ulipt = testSite;
		ulipt++;
	}
	// has to be sorted for later use:
	myInsertSort( mutatedLoci, nNewMutations );
	
#ifdef DEBUG
	if ( VERBOSE ) {
		if ( t % TIME_SERIES_SAMPLE_FREQ == 0 ) {
			fprintf(stdout, "\nMutation sites in time step %lu\n", t);
			for ( i = 0; i < nNewMutations; i++ ) {
				fprintf(stdout, "\t%lu", mutatedLoci[i]);
			}
			fprintf(stdout, "\n");
		}
	}
#endif
}



void chooseParentsAtRandom(long int *mommy, long int *daddy, long int *randomNumberLine, long int nhere)
{
    long int chosenOnes[2];
    
    gsl_ran_choose( rngState, chosenOnes, 2, randomNumberLine, nhere, sizeof(long int) );
	if ( t == 1 ) {
		fprintf(stdout, "\nYo! gsl_ran_choose() is really slow for this.  Rewrite chooseParentsAtRandom()!\n");
	}
    
    *mommy = chosenOnes[0];
    *daddy = chosenOnes[1];
}


void computeFitness(double *fitnessValues)
{
    unsigned long int i, j, selectedSiteIndexesMaster[nSelectedSites];
    unsigned long int selectedSiteIndexesLocal[nSelectedSites], *ullpt;
    double cf_scoeffs[nSelectedSites], *dpt, cf_s, cf_gradient; // prefix cf_ denoting local to this function
    long int nFound;
    int cf_siteClass, cf_siteClasses[nSelectedSites], *locpt;
    short int *sipt, gtsum;
//    unsigned long int *dbpt1, *dbpt2; // pointers for debugging with lldb
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
            printf("\t%lu\t%lu\t%i\t%f\n", selectedSiteIndexesLocal[i], selectedSiteIndexesMaster[i], cf_siteClasses[i], cf_scoeffs[i]);
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
                        *dpt *= 1.0 - (2.0 * cf_s * cf_gradient); // note cf_gradient can be + or -
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
            fprintf(cf_fitvals, "%lu,%i,%f\n", i, *(locations + i), *(fitnessValues + i));
        fclose(cf_fitvals);
    }
    dpt = fitnessValues;
	// correct for negative fitness values if they occur:
    for ( i = 0; i < N; i++ ) {
        if ( *dpt < 0.0 )
            *dpt = 0.0;
        dpt++;
    }
}


long int figureOutOffspringGenomeSites( unsigned long int *offsp_SiteIndexes, short int *offsp_lociStates, long int nNewMutations, unsigned long int *mutatedLoci, _Bool *copyFromParents )
{
    long int nSitesInOffspring = 0, newMutationSiteCount = 0, parentalSiteCount = 0, nRepeatMutations = 0;
    unsigned long int locus, *newLocusIndexPt, *variableSitePt, *offsp_pt;
    long int i, totalSitesDone = 0, retainedSites = 0, nSitesToDo, nBrandNewMutations = 0;
    short int *offsp_ls;
	_Bool *cfppt;

    newLocusIndexPt = mutatedLoci;
    variableSitePt = parentalTrackedSiteIndexes;
    offsp_pt = offsp_SiteIndexes;
    offsp_ls = offsp_lociStates;
	cfppt = copyFromParents;
    
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
					*cfppt = 1; // should be copied from parents
					cfppt++;
                    
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
				*cfppt = 0; // should NOT be copied from parents
				cfppt++;
                
                newLocusIndexPt++; // advance new locus index pointer
                
                nSitesInOffspring++; // increment the total site counter
                
                newMutationSiteCount++;
                
                totalSitesDone++;
                
                nBrandNewMutations++;
            }
            else {
                // the mutation is at an already used site
                if ( *newLocusIndexPt != *variableSitePt ) {
                    fprintf(stderr, "\nError in figureOutOffspringGenomeSites():\n\t*newLocusIndexPt (%lu) != *variableSitePt (%lu)\n", *newLocusIndexPt, *variableSitePt);
                    exit(-1);
                }
                *offsp_pt = *newLocusIndexPt; // record the index
                offsp_pt++; // advance the index pointer
                
                *offsp_ls = LOCUS_STATUS_VARIABLE_PLUS_MUT; // record the status
                offsp_ls++; // advance the status pointer
				*cfppt = 1; // should be copied from parents
				cfppt++;
                
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
				*cfppt = 1; // should be copied from parents
				cfppt++;
                
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
			*cfppt = 0; // should NOT be copied from parents
			cfppt++;
            
            newLocusIndexPt++; // advance new locus index pointer
            
            nSitesInOffspring++; // increment the total site counter
            
            newMutationSiteCount++;
            
            totalSitesDone++;
            
            nBrandNewMutations++;
        }
        else {
            fprintf(stderr, "\nError in figureOutOffspringGenomeSites():\n\tcounts aren't working how you think\n\t");
            fprintf(stderr, "parental count = %li\tnewMutationSiteCount = %li\n\tnTrackedSitesInParents = %lu, nNewMutations = %li\n", parentalSiteCount, newMutationSiteCount, nTrackedSitesInParents, nNewMutations);
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
            printf("\t%li\t%lu\t%i\n", i, locus, sitesStatuses[locus]);
        }
        printf("\nOffspring loci and status codes:\n\tnumber\tindex\tcode\n");
        for ( i = 0; i < nSitesInOffspring; i++ ) {
            locus = *(offsp_SiteIndexes + i);
            printf("\t%li\t%lu\t%i\n", i, locus, offsp_lociStates[i]);
        }
        for ( i = 1; i < nSitesInOffspring; i++ ) {
            if ( *(offsp_SiteIndexes + i) <= *(offsp_SiteIndexes + i - 1) ) {
                fprintf(stderr, "\nError in figureOutOffspringGenomeSites():\n\t*(offsp_SiteIndexes + i) (%lu) <= *(offsp_SiteIndexes + i - 1) (%lu)\n", *(offsp_SiteIndexes + i), *(offsp_SiteIndexes + i - 1));
                exit(-1);
            }
        }
        
        printf("\nnSitesInOffspring = %li\nNew sites from mutations = %li, retained sites = %li, repeat mutations = %li\n", nSitesInOffspring, newMutationSiteCount, retainedSites, nRepeatMutations);
    }
    
    return nSitesInOffspring;
}



void figureOutParentalSkippedSites( _Bool *skipTheseParentalSites )
{
	// determine which parental sites (from parent perspective), should be skipped in copying because they are not variable:
	long int i;
	_Bool *bpt = skipTheseParentalSites;
	unsigned long int locus, *ulipt = parentalTrackedSiteIndexes;
	
	for ( i = 0; i < nTrackedSitesInParents; i++ ) {
		locus = *ulipt; // master index of locus
		if ( *(sitesStatuses + locus) != LOCUS_STATUS_VARIABLE_IN_PARENTS )
			*bpt = 1; // skip it
		else
			*bpt = 0; // do not skip it
		bpt++;
		ulipt++;
	}
}


void finalTasks(unsigned RNG_SEED)
{
    printParametersToFiles(RNG_SEED);
    
    fclose(dataFile_alleleFreqTS);
//    fclose(dataFile_alleleFreqTSbyPop);
    fclose(dataFile_SFS_TS);
    fclose(dataFile_segSiteTS);
    fclose(dataFile_derivedFixationTS);
	fclose(dataFile_PiAndDXY);
	fclose(dataFile_abundances);
    
    free(genotypes0);
    free(genotypes1);
    free(parentalTrackedSiteIndexes);
//    free(siteIndexes);
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
    free(expLookupTable);
	
	fprintf(stdout, "\nTotal mutations that arose during run:  %lu\n", totalMutationsForRun);
	
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
			*(dptfnl - 1) = 1.1; // safety to prevent overrun of number lines
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
	// demePosition holds the indexes of where to jump into the array for each population
    
    
    ipt = locations;  // locations is the global vector giving the deme number of each individual
    
    for ( i = 0; i < N; i++ ) {
        loc = *ipt; // deme number
        spot = demePosition[loc]; // where it belongs in the individualsInDeme array
        *(individualsInDeme + spot) = i;
        demePosition[loc] = demePosition[loc] + 1;
        
        ipt++;
    }

#ifdef DEBUG
	FILE *testMakeIndexes;
	testMakeIndexes = fopen("TestMakeDemesIndexes.txt","w");
	for ( i = 0; i < N; i++ ) {
		j = individualsInDeme[i];
		fprintf(testMakeIndexes, "%li,%i\n", j, locations[j]);
	}
	fclose(testMakeIndexes);
#endif
}


//void makeOneOffspring(long int momIndex, long int dadIndex, short int *offGTpt, long int nSitesInOffspring, unsigned long int *offsp_SiteIndexes, short int *offsp_lociStates)
void makeOneOffspring(long int momIndex, long int dadIndex, short int *offGTpt, long int nSitesInOffspring, unsigned long int *offsp_SiteIndexes, _Bool *copyFromParents, _Bool *skipTheseParentalSites, unsigned long int *focalSites, int *siteLGs)

{
    int i, currentLinkageGroup, thisSiteLG;
	unsigned long int j, focalSite, *offsp_SIpt;
	//unsigned long int *parentalLocusIndexes;
	short int *sipt, *parentPoint;
	// short int *offsp_ls;
    int chromosome;
    double meanRecombDistance = 1000.0 / RECOMBINATION_RATE_PER_KB;
    unsigned long int nextRecombinationSpot;
    long int parentalLocusCounter;
	_Bool *cfppt, *psspt;
	
#ifdef DEBUG
        if ( *(locations + momIndex) != *(locations + dadIndex) ) {
            printf("\nError in makeOneOffspring:\n\tIndexes are off giving bogus mom or dad\n");
            exit(-1);
        }
#endif
	
    for ( i = 0; i < 2; i++ ) {
 
MARK(MOO_START)
        
        sipt = offGTpt + i; // which haploid set in offspring
		
		psspt = skipTheseParentalSites; // pointer to bools for easy figuring out skips of sites
					// takes place of *(sitesStatuses + *parentalLocusIndexes) in former code versions
					// counting should be the same as old "parentalLocusIndexes"
        
        parentalLocusCounter = 0;
        //parentalLocusIndexes = parentalTrackedSiteIndexes;
        if ( nTrackedSitesInParents ) {
            //while ( *(sitesStatuses + *parentalLocusIndexes) != LOCUS_STATUS_VARIABLE_IN_PARENTS ) {
			while ( *psspt ) {
                parentalLocusCounter++; // so we can start at the first site (locus) that needs to be inherited
                //parentalLocusIndexes++;
				psspt++;
            }
        }
        
MARK(MOO_FIRSTLOCUS)
        
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

        //offsp_ls = offsp_lociStates;
		cfppt = copyFromParents;	// pointer to easy array of bools telling if a locus should be copied
									// cfppt is from offspring's perspective; replaces offsp_ls in previous versions
		
		offsp_SIpt = offsp_SiteIndexes;
        nextRecombinationSpot = (unsigned long int) randExp( meanRecombDistance );


        for ( j = 0; j < nSitesInOffspring; j++ ) {
            focalSite = focalSites[j];
MARK(MOO_STARTFORLOOP)
            // now copy alleles from parents to offspring
            //if ( *offsp_ls == LOCUS_STATUS_VARIABLE_IN_PARENTS || *offsp_ls == LOCUS_STATUS_VARIABLE_PLUS_MUT ) {
			if ( *cfppt ) {
				// yes, copy it from parents:
				thisSiteLG = siteLGs[j];
				
                // first handle recombination and independent assortment
                if ( thisSiteLG != currentLinkageGroup ) {
MARK(MOO_INSIDE_IF)
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
#ifdef DEBUG
                    if ( thisSiteLG < currentLinkageGroup ) {
                        fprintf(stderr, "\nError in makeOneOffspring():\n\t*(linkageGroupMembership + focalSite) (%i) < currentLinkageGroup (%i)\n", thisSiteLG, currentLinkageGroup );
                        exit(-1);
                    }
#endif
                    currentLinkageGroup = thisSiteLG;
                    
                    nextRecombinationSpot = focalSite + ((unsigned long int) randExp( meanRecombDistance ));
                }
                else if ( focalSite > nextRecombinationSpot ) {
MARK(MOO_INSIDE_ELSEIF)
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
                        nextRecombinationSpot += ((unsigned long int) randExp( meanRecombDistance ));
                    } while ( focalSite > nextRecombinationSpot );
                }
MARK(MOO_POINTERBUSINESS)
                *sipt = *parentPoint;
                parentalLocusCounter++;  // one more parental locus taken care of
MARK(MOO_WHILELOOP)
                if ( j < (nSitesInOffspring - 1) && parentalLocusCounter < nTrackedSitesInParents ) {
                    //parentalLocusIndexes++; // advance to index of next parental site/locus
					psspt++; // next skip site bool
                    parentPoint += PLOIDY; // advance to next site in genome
                    //while ( *(sitesStatuses + *parentalLocusIndexes) != LOCUS_STATUS_VARIABLE_IN_PARENTS && parentalLocusCounter < nTrackedSitesInParents ) {
					while ( *psspt && parentalLocusCounter < nTrackedSitesInParents ) {
						parentalLocusCounter++; // count prior locus toward those that are "done"
                        // the while part helps skip over loci that will no longer be tracked
                        if ( parentalLocusCounter < nTrackedSitesInParents ) {
                            // if statement prevents seg fault on last While check
                            //parentalLocusIndexes++; // advance to index of next parental site/locus
                            parentPoint += PLOIDY; // advance to next site in genome
							psspt++; // increment bool skip site pointer
                        }
                    }
                }
            }
            else {
MARK(MOO_ELSE)
                *sipt = ALLELE_CODE_ANCESTRAL; // just put in the locus as a placeholder; mutations added later
				// do NOT advance any parent pointers or counters
            }
            if ( *sipt )
                *(alleleCounts + focalSite) += 1; // record derived allele counts
            
            sipt += PLOIDY; // move to next locus in offspring haplotype
            offsp_SIpt++; // next site index in offspring
            //offsp_ls++; // next locus state
			cfppt++; // increment bool for copy from parents
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


void myInsertSort(unsigned long int * array, long int numElements) {
	
	// code modified from: http://www.programmingsimplified.com/c/source-code/c-program-insertion-sort
	// accessed 6/9/17
	long int c, d, d1;
	double foo;
	
	for (c = 1; c < numElements; c++) {
		
		d = c;
		d1 = d - 1;
		
		while ( d > 0 && *(array + d) < *(array + d1) ) {
			foo = *(array + d);
			*(array + d) = *(array + d1);
			*(array + d1) = foo;
			d--;
			d1--;
		}
		
	}
}



void printParametersToFiles(unsigned RNG_SEED)
{
    // make an R-sourceable script of parameter values
    int i, j;
    FILE *rfile;
    double *dpt, THETA;
    unsigned long int foo;
    
    rfile = fopen("MetadataAndParameters.R","w");
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
    fprintf(rfile, "# genome parameters\nnSITES <- %lu\n", nSITES);
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
    fprintf(rfile, "nTrackedSitesInParents <- %lu\n", nTrackedSitesInParents);
    fprintf(rfile, "abundances <- c(%li", abundances[0]);
    for ( i = 1; i < nPOPULATIONS; i++ )
        fprintf(rfile, ",%li", abundances[i]);
    fprintf(rfile, ")\n");
	fprintf(rfile, "totalMutationsForRun <- %lu\n", totalMutationsForRun);
    
    fprintf(rfile, "\n");
    
    fclose(rfile);
}


void putInMutations( short int *offspringGTs, short int *offsp_lociStates, unsigned long int *offsp_SiteIndexes, long int nNewMutations, long int totalOffspring, long int nSitesInOffspring )
{
    long int i, individual, totalMutsAdded = 0;
    short int *spot, *offsp_ls;
    unsigned long int locus;
    
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
                fprintf(stderr, "\nError in putInMutations():\n\t*(alleleCounts + locus) (%lu) out of bounds\n", *(alleleCounts + locus));
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
	// potential optimization with lookup table??
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
    unsigned long int locus, *ullipt;
    unsigned long int *offsp_SiteIndexes;
    
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
	// potential optimization with lookup table??
    nNewMutations = gsl_ran_poisson( rngState, (GENOME_MU * ((double) totalOffspring)) );
	totalMutationsForRun += nNewMutations;
	
    unsigned long int mutatedLoci[nNewMutations];
	// next was hugely expensive:
    //gsl_ran_choose( rngState, mutatedLoci, nNewMutations, siteIndexes, nSITES, sizeof(unsigned long int) );
	// trying my own choice function instead:
	chooseMutationSites( mutatedLoci, nNewMutations );
    maxSitesInOffspring = nTrackedSitesInParents + nNewMutations;
    short int offsp_lociStates[maxSitesInOffspring];
    offsp_SiteIndexes = (unsigned long int *) malloc( maxSitesInOffspring * sizeof( unsigned long int ) );
	_Bool copyFromParents[maxSitesInOffspring], skipTheseParentalSites[nTrackedSitesInParents];
    
    // determine which sites will need to be tracked among offspring
    nSitesInOffspring = figureOutOffspringGenomeSites( offsp_SiteIndexes, offsp_lociStates, nNewMutations, mutatedLoci, copyFromParents );
	// determine which parental sites (from parent perspective), should be skipped in copying because they are not variable:
	figureOutParentalSkippedSites( skipTheseParentalSites );
    
    // make some useful vectors ONCE to prevent jumping all over memory in makeOneOffspring():
    unsigned long int focalSites[nSitesInOffspring], fsfoo;
    int siteLGs[nSitesInOffspring];
    for ( j = 0; j < nSitesInOffspring; j++ ) {
        fsfoo = *(offsp_SiteIndexes + j);
        focalSites[j] = fsfoo;
        siteLGs[j] = *(linkageGroupMembership + fsfoo);
    }
    
    // debug
    //printf("\nreproduction: time, nTrackedSitesInParents, nSitesInOffspring, nNewMutations, N\n\t%li\t%lu\t%li\t%li\t%li", t, nTrackedSitesInParents, nSitesInOffspring, nNewMutations, N);
    
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
    memset( alleleCounts, 0, nSITES * sizeof(unsigned long int) );
    
    // some debug code
    /*
    if ( t >= 10 ) {
        printf("\nTracked sites in parents (%lu) @ t = %li:\n", nTrackedSitesInParents, t);
        for ( i = 0; i < nTrackedSitesInParents; i++ )
            printf("%lu\t", parentalTrackedSiteIndexes[i]);
        
        printf("\n\nTSites in offspring (%li):\n", nSitesInOffspring);
        for ( i = 0; i < nSitesInOffspring; i++ )
            printf("%lu\t", offsp_SiteIndexes[i]);
        printf("\n\n");
    }
     */
    
    
    
    for ( pop = 0; pop < nPOPULATIONS; pop++ ) { // one population at a time
        noff = noffspring[pop]; // total number of offspring to be born in this population
        nhere = abundances[pop]; // number of potential parents here
        if ( noff > 0 ) {
            for ( i = 0; i < noff; i++ ) {
				if ( INCLUDE_SELECTION ) {
					// get mom first:
					mommy = arraySearch(dpt, nhere);
					// now dad:
					do {
						daddy = arraySearch(dpt, nhere);
					} while ( mommy == daddy );
				}
                else
                    chooseParentsAtRandom(&mommy, &daddy, randomNumberLine, nhere);
                momIndex = *(lipt + mommy); // actual index of mom (again, since there is no sort)
                dadIndex = *(lipt + daddy);
                //makeOneOffspring(momIndex, dadIndex, sipt, nSitesInOffspring, offsp_SiteIndexes, offsp_lociStates);
				makeOneOffspring(momIndex, dadIndex, sipt, nSitesInOffspring, offsp_SiteIndexes, copyFromParents, skipTheseParentalSites, focalSites, siteLGs);
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
	// update abundances:
    for ( i = 0; i < nPOPULATIONS; i++ ) {
        abundances[i] = noffspring[i];
    }
	// update total across all demes:
    N = totalOffspring;
	
	// time to see what's happening at each tracked site:
    ullipt = offsp_SiteIndexes;
	// inactive by default (note that this covers ALL sites):
    memset( sitesStatuses, LOCUS_STATUS_INACTIVE, nSITES * sizeof(short int) );
    nSelectedSites = 0;
    for ( i = 0; i < nSitesInOffspring; i++ ) {
        locus = *ullipt; // site/locus index
		// offspring are now the parents for the next generation:
		*(sitesStatuses + locus) = LOCUS_STATUS_TRACKED_IN_PARENTS;
        if ( *(siteClassifications + locus) != SITE_CLASS_NEUTRAL )
            nSelectedSites++;
        ullipt++;
    }
	
	// now account for chance loss and actuall variable sites
    ullipt = alleleCounts;
    sipt = sitesStatuses;
    driftLosses = 0;
    actuallyVariable = 0;
	// we will look over ALL sites:
    for ( locus = 0; locus < nSITES; locus++ ) {
        if ( *ullipt > 0 && *ullipt < ( PLOIDY * N ) ) {
			// there are between 1 and 2N-1 copies --> must be variable:
			*sipt = LOCUS_STATUS_VARIABLE_IN_PARENTS;
            actuallyVariable++;
        }
        else if ( *ullipt == (PLOIDY * N) ) {
            fprintf( dataFile_derivedFixationTS, "%li,%lu,%i,%E\n", t, locus, *(siteClassifications + locus), *(selectionCoefficients + locus) );
            nDerivedFixations++;
        }
        else if ( *sipt == LOCUS_STATUS_TRACKED_IN_PARENTS ) {
            if ( VERBOSE )
                printf("\nLocus %lu lost in time step %li\n", locus, t);
            driftLosses++;
        }
        sipt++;
        ullipt++;
    }
    
    if ( (driftLosses + actuallyVariable + nDerivedFixations) != nSitesInOffspring ) {
        fprintf(stderr, "\nError in reproduction():\n\t(driftLosses (%li) + actuallyVariable (%li)) + nDerivedFixations (%i) != nSitesInOffspring (%li) )\n", driftLosses, actuallyVariable, nDerivedFixations, nSitesInOffspring);
        exit(-1);
    }
	
	// offspring become parents:
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
    
    fprintf(stdout, "\nNote that very limited options are available here on the command line.\nFor extensive flexible settings, put settings in parameters.ini.txt.\nOptions available on command line are as follows:\n");
    
    fprintf(stdout, "\n\t-V\tVERBOSE: print lots of human readable messages to\n\t\tstdout related to run status.\n");
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







