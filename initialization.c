// initialization functions used in fffits

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
#include <gsl_randist.h>


// read in optional command line arguments
int initializationSteps( int argc, char *argv[], char *progname )
{
	unsigned int RNG_SEED;
	int ch;
	
	while ((ch = getopt(argc, argv, "T:V?")) != -1) {
		switch (ch) {
			case 'T':
				nGENERATIONS = atoi(optarg);
				break;
			case 'V':
				VERBOSE = 1;
				break;
			case '?':
			default:
				usage(progname);
				exit(-1);
		}
	}
	
	RNG_SEED = readInParametersFromFile();
	
	// initialization steps
	initializeRNG(RNG_SEED);
	setUpGenome();
	setUpPopulations();
	setUpDataFiles();
	
	return ( RNG_SEED );
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


void setUpDataFiles(void)
{
	int i;
	
	dataFile_alleleFreqTS = fopen("AlleleFreqTS.csv", "w");
	fprintf(dataFile_alleleFreqTS, "Time,SiteIndex,LinkageGroup,DerivedAlleleCount,DerivedAlleleFreq,SelectionCoefficient,SiteClassCode\n");
	
	dataFile_alleleFreqTSbyPop = fopen("AlleleFreqTSbyPop.csv", "w");
	fprintf(dataFile_alleleFreqTSbyPop, "Time,SiteIndex");
	for ( i = 0; i < nPOPULATIONS; i++ )
		fprintf(dataFile_alleleFreqTSbyPop, ",CountInPop%i", i);
	fprintf(dataFile_alleleFreqTSbyPop, "\n");
	
	dataFile_SFS_TS = fopen("SFStimeSeries.csv", "w");
	fprintf(dataFile_SFS_TS, "Time,DerivedAlleleCopyNumber,NumberOfSites\n");
	
	dataFile_segSiteTS = fopen("SegregatingSitesTS.csv", "w");
	fprintf(dataFile_segSiteTS, "Time,nSegregatingSites,nNeutralSites,nBackgroundSelSites,nPositiveSelSites,nDivergentSelSites\n");
	
	dataFile_derivedFixationTS = fopen("DerivedFixationRecord.csv", "w");
	fprintf(dataFile_derivedFixationTS, "Time,SiteIndex,SiteClassCode,SelectionCoefficient\n");
}


void setUpGenome(void)
{
	double theta, expectedSegSites, *expectedFreq, *dpt, value;
	unsigned long long int foo, *ullpt, sitesPerLinkageGroup, lgCount, SFScounts[(PLOIDY * N)];
	long int i;
	int dumi, currentLinkageGroup;
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
	nTrackedSitesInParents = (unsigned long long int) (expectedSegSites + 0.5);
	
	
	// now to choose the sites that are variable and their frequencies
	siteIndexes = (unsigned long long int *) malloc( nSITES * sizeof(unsigned long long int) );
	parentalTrackedSiteIndexes = (unsigned long long int *) malloc( nTrackedSitesInParents * sizeof(unsigned long long int) );
	ullpt = siteIndexes;
	for ( foo = 0; foo < nSITES; foo++ ) {
		*ullpt = foo;
		ullpt++;
	}
	
	// use gsl ran choose to pick sites at random to be standing neutral variation
	gsl_ran_choose( rngState, parentalTrackedSiteIndexes, nTrackedSitesInParents, siteIndexes, nSITES, sizeof(unsigned long long int) );
	// variable site indexes now stores the ordered indexes of sites that should have segregating variants
	// siteIndexes just stores the indexes [0..nSITES-1], which is useful for random choices like this one, so we keep it
	sitesStatuses = (short int *) malloc( nSITES * sizeof(short int) );
	siteClassifications = (int *) malloc( nSITES * sizeof(int) );
	linkageGroupMembership = (int *) malloc( nSITES * sizeof(int) );
	memset( sitesStatuses, LOCUS_STATUS_INACTIVE, (nSITES * sizeof(short int)) );
	for ( foo = 0; foo < nTrackedSitesInParents; foo++ )
		*(sitesStatuses + (parentalTrackedSiteIndexes[foo])) = LOCUS_STATUS_VARIABLE_IN_PARENTS;
	// sitesStatuses stores the same information as parentalTrackedSiteIndexes, but in a different form; both are useful
	
	// now randomly assign frequencies
	alleleFrequencies = (double *) malloc( nSITES * sizeof(double)); // frequencies (0 to 1) by site
	alleleCounts = (unsigned long long int *) malloc( nSITES * sizeof(unsigned long long int));
	memset( alleleCounts, 0, nSITES * sizeof(unsigned long long int) );
	for ( foo = 0; foo < (PLOIDY * N); foo++ )
		*(SFScounts + foo) = 0;
	setUpInitialAlleleFrequencies(expectedFreq, SFScounts);
	
	// pick allele frequencies using the Ewens sampling formula from Wakeley's book.
	
	
	if ( VERBOSE )
		printf("\ntheta = %f, segSites = %llu\n", theta, nTrackedSitesInParents);
	
	
	// assign site types
	double siteProbs[nSITE_CLASSES];
	unsigned long long int siteTypeCounts[nSITE_CLASSES];
	memset( siteTypeCounts, 0, (nSITE_CLASSES * sizeof(unsigned long long int)) );
	selectionCoefficients = (double *) malloc( nSITES * sizeof(double) );
	memset( selectionCoefficients, 0, (nSITES * sizeof(double)) );
	siteProbs[0] = PROBABILITY_SITE_NEUTRAL;
	siteProbs[1] = siteProbs[0] + PROBABILITY_SITE_BGS;
	siteProbs[2] = siteProbs[1] + PROBABILITY_SITE_POS;
	siteProbs[3] = 1.0;
	currentLinkageGroup = 0;
	lgCount = 0;
	sitesPerLinkageGroup = (int) ( ((double) nSITES) / ((double) nLINKAGE_GROUPS) );
	for ( foo = 0; foo < nSITES; foo ++ ) {
		
		// assign site types and selection coefficients
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
		
		// assign linkage group membership
		*(linkageGroupMembership + foo) = currentLinkageGroup;
		lgCount++;
		if ( lgCount >= sitesPerLinkageGroup && currentLinkageGroup < (nLINKAGE_GROUPS - 1) ) {
			lgCount = 0;
			currentLinkageGroup++;
		}
	}
	
	foo = 0;
	if ( VERBOSE )
		printf("\nSite type counts:\n");
	for ( i = 0; i < nSITE_CLASSES; i++ ) {
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
	for ( i = 1; i < nSITE_CLASSES; i++ ) {
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
	fprintf(initialFreqs, "SiteIndex,LinkageGroup,DerivedAlleleCount,DerivedAlleleFreq,SelectionCoefficient,SiteClassCode,SiteClassName\n");
	for ( foo = 0; foo < nTrackedSitesInParents; foo++ ) {
		focalSiteIndex = *(parentalTrackedSiteIndexes + foo);
		
		if ( *(siteClassifications + focalSiteIndex) > SITE_CLASS_NEUTRAL )
			nSelectedSites++;
		
		
		fprintf(initialFreqs, "%llu,%i,%llu,%E,%E", focalSiteIndex, linkageGroupMembership[focalSiteIndex], alleleCounts[focalSiteIndex], alleleFrequencies[focalSiteIndex], selectionCoefficients[focalSiteIndex] );
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


void setUpInitialAlleleFrequencies(double *expectedFreq, unsigned long long int *SFScounts)
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
	
	
	for ( i = 0; i < nTrackedSitesInParents; i++ ) {
		
		// get index of focal site
		focalSiteIndex = *(parentalTrackedSiteIndexes + i);
		
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
	neededSize = PLOIDY * N * sizeof(short int) * nTrackedSitesInParents;
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
			/* the sequence of choices needs to be: 0,1, PLOIDY * nTrackedSitesInParents, PLOIDY*nTrackedSitesInParents + 1, 2*PLOIDY*nTrackedSitesInParents, 2*PLOIDY*nTrackedSitesInParents + 1, etc... */
			choiceVector[counter] = (j * PLOIDY * nTrackedSitesInParents) + i;
			counter++;
		}
	}
	
	for ( j = 0; j < nTrackedSitesInParents; j++ ) {
		focalSiteIndex = parentalTrackedSiteIndexes[j];
		siteAlleleCount = alleleCounts[focalSiteIndex];
		if ( siteAlleleCount <= 0 || siteAlleleCount >= (PLOIDY * N) ) {
			printf("\nError in setUpPopulations():\n\tsiteAlleleCount (= %llu) out of bounds.\n", siteAlleleCount);
			exit(-1);
		}
		
		// int gsl_ran_choose (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size)
		gsl_ran_choose( rngState, allelesToSwitch, siteAlleleCount, choiceVector, (PLOIDY * N), sizeof(unsigned long long int) );
		
		for ( k = 0; k < siteAlleleCount; k++ ) {
			// index into genotypes array
			// genotypes array: each individual = PLOIDY * nTrackedSitesInParents consecutive entries
			// two adjacent entries are the two entries for the two alleles an individual has at that locus
			// need to add "1" allele at jth locus for each of selected genotypes that get a derived allele copy
			// hence, index of jth locus in kth genotype copy is:
			index = (PLOIDY * j) + allelesToSwitch[k]; // kth genotype copy of jth locus to add allele to
			*(gts + index) = ALLELE_CODE_DERIVED; // derived allele
		}
		
	}
	
	
	FILE *initgts;
	if ( VERBOSE ) {
		initgts = fopen("InitialGenotypes.csv", "w");
		fprintf(initgts, "Locus0Copy0,Locus0Copy1");
		for ( k = 1; k < nTrackedSitesInParents; k++ ) {
			for ( i = 0; i < PLOIDY; i++ ) {
				fprintf(initgts, ",Locus%lluCopy%i", k, i);
			}
		}
		fprintf(initgts,"\n");
		
		sipt = gts;
		for ( k = 0; k < N; k++ ) {
			fprintf(initgts, "%i", *sipt);
			sipt++;
			for ( j = 1; j < (PLOIDY * nTrackedSitesInParents); j++ ) {
				fprintf(initgts, ",%i", *sipt);
				sipt++;
			}
			fprintf(initgts, "\n");
		}
		fclose(initgts);
	}
	
	environmentGradient = (double *) malloc( nPOPULATIONS * sizeof(double) );
	if ( ENVIRONMENT_TYPE == ENVT_TYPE_GRADIENT )
		stepSize = (ENVT_MAX - ENVT_MIN) / (((double) nPOPULATIONS) - 1.0);
	for ( i = 0; i < nPOPULATIONS; i++ ) {
		if ( ENVIRONMENT_TYPE == ENVT_TYPE_GRADIENT ) {
			environmentGradient[i] = ENVT_MIN + (((double) i) * stepSize);
			if ( i == 0 )
				fprintf(stdout, "\nNote from setUpPopulations():\n\tEnvironment gradient algorithm assumes a one-dimensional habitat\n");
		}
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


void wrongParametersIniOption(char *expected, char *previous, char *found)
{
	fprintf(stderr, "\nError in readInParametersFromFile():\n\t");
	fprintf(stderr, "Expecting '%s' after '%s'\n\t", expected, previous);
	fprintf(stderr, "in parameters.ini.txt, but instead found: '%s'.\n\tPlease fix parameters.ini.txt\n", found);
	fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
	exit(-1);
}







