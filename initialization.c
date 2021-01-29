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
#include <limits.h>
#include "fffits.h"
#include "pcg-c-master/include/pcg_variants.h"

#include <gsl_rng.h>            // gnu scientific library //
#include <gsl_randist.h>


// read in optional command line arguments
int initializationSteps( int argc, char *argv[], char *progname )
{
	unsigned int RNG_SEED;
	int ch;
	
	while ((ch = getopt(argc, argv, "V?")) != -1) {
		switch (ch) {
//			case 't':
//				nGENERATIONS = atoi(optarg);
//				break;
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

	// data type sizes allowed:
	if ( VERBOSE ) {
		fprintf(stdout, "\nULONG_MAX on this system = %lu\n", ULONG_MAX);
		fprintf(stdout, "\t--> nSITES must be <= %lu\n", ULONG_MAX);
	}
	
	// initialization steps
	initializeRNG(RNG_SEED);
    initializePCGRNG(RNG_SEED);
	setUpGenome();
	setUpPopulations();
	setUpDataFiles();
    // makeExpLookupTable(); // a method tried to replace randExp()
	
	return ( RNG_SEED );
}


void initializeRNG(unsigned int RNG_SEED)
{
	
	// initialize random number generator
	// T and r are global variables for the RNG state and calls
	gsl_rng_env_setup(); // set up the environment variables for the RNG
	rngType = gsl_rng_mt19937; // Mersenne Twister// for default you can say T = gsl_rng_default;
	rngState = gsl_rng_alloc(rngType);
	gsl_rng_set(rngState, RNG_SEED);
	
    // int i;
	// for ( i = 0; i < 10; i++ )
	//     printf("%f\n", gsl_rng_uniform(rngState) );
}


void initializePCGRNG(unsigned int RNG_SEED)
{
    pcg32_srandom(42u, RNG_SEED);
    
    printf("  32bit:");
    for (int i = 0; i < 6; ++i)
        printf("%u\t", pcg32_random());
    printf("\n");
}


unsigned readInParametersFromFile(void)
{
	char c, option[80], option2[80];
	long int INITIAL_N;
	unsigned long int foo;
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
			fscanf(pfile, "%lu", &nSITES);
			if ( nSITES < 1 ) {
				fprintf(stderr, "\nError in readInParametersFromFile():\n\tnSITES (= %lu) should be > 1\n", nSITES);
				fprintf(stderr, "Please fix parameters.ini.txt\n");
				fprintf(stderr, "\n\t\t*** Exiting *** \n\n");
				exit(-1);
			}
			// test check
			if ( VERBOSE )
				printf("Found nSITES (%s) = %lu\n", option, nSITES);
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


void setUpDataFiles(void)
{
	int i, j;
	
	dataFile_alleleFreqTS = fopen("AlleleFreqTS.csv", "w");
	fprintf(dataFile_alleleFreqTS, "Time,SiteIndex,LinkageGroup,DerivedAlleleCount,DerivedAlleleFreq,SelectionCoefficient,SiteClassCode");
	for ( i = 0; i < nPOPULATIONS; i++ )
		fprintf(dataFile_alleleFreqTS, ",CountInPop%i", i);
	fprintf(dataFile_alleleFreqTS, ",FST\n");
	
//	dataFile_alleleFreqTSbyPop = fopen("AlleleFreqTSbyPop.csv", "w");
//	fprintf(dataFile_alleleFreqTSbyPop, "Time,SiteIndex");
//	for ( i = 0; i < nPOPULATIONS; i++ )
//		fprintf(dataFile_alleleFreqTSbyPop, ",CountInPop%i", i);
//	fprintf(dataFile_alleleFreqTSbyPop, "\n");
	
	dataFile_SFS_TS = fopen("SFStimeSeries.csv", "w");
	fprintf(dataFile_SFS_TS, "Time,DerivedAlleleCopyNumber,NumberOfSites\n");
	
	dataFile_segSiteTS = fopen("SegregatingSitesTS.csv", "w");
	fprintf(dataFile_segSiteTS, "Time,nSegregatingSites,nNeutralSites,nBackgroundSelSites,nPositiveSelSites,nDivergentSelSites\n");
	
	dataFile_derivedFixationTS = fopen("DerivedFixationRecord.csv", "w");
	fprintf(dataFile_derivedFixationTS, "Time,SiteIndex,SiteClassCode,SelectionCoefficient\n");
	
	dataFile_abundances = fopen("Abundances.csv", "w");
	fprintf(dataFile_abundances, "Time,Ntotal");
	for ( i = 0; i < nPOPULATIONS; i++ )
		fprintf(dataFile_abundances, ",AbundancePop%i", i);
	fprintf(dataFile_abundances, "\n");
	
	dataFile_PiAndDXY = fopen("PiAndDXY.csv", "w");
	fprintf(dataFile_PiAndDXY, "Time");
	fprintf(dataFile_PiAndDXY, ",TotalGlobalPi,GlobalNEUTPi,GlobalBGPi,GlobalPOSPi,GlobalDIVPi");
	for ( i = 0; i < nPOPULATIONS; i++ )
		fprintf(dataFile_PiAndDXY, ",NEUTPiPop%i,BGPiPop%i,POSPiPop%i,DIVPiPop%i", i, i, i, i);
	for ( i = 0; i < (nPOPULATIONS - 1); i++ )
		for ( j = (i+1); j < nPOPULATIONS; j++ )
			fprintf(dataFile_PiAndDXY, ",TotalNEUTDxyBetween%iand%i", i, j);
	for ( i = 0; i < (nPOPULATIONS - 1); i++ )
		for ( j = (i+1); j < nPOPULATIONS; j++ )
			fprintf(dataFile_PiAndDXY, ",TotalBGDxyBetween%iand%i", i, j);
	for ( i = 0; i < (nPOPULATIONS - 1); i++ )
		for ( j = (i+1); j < nPOPULATIONS; j++ )
			fprintf(dataFile_PiAndDXY, ",TotalPOSDxyBetween%iand%i", i, j);
	for ( i = 0; i < (nPOPULATIONS - 1); i++ )
		for ( j = (i+1); j < nPOPULATIONS; j++ )
			fprintf(dataFile_PiAndDXY, ",TotalDIVDxyBetween%iand%i", i, j);
	fprintf(dataFile_PiAndDXY, "\n");

}


void setUpGenome(void)
{
	double theta, expectedSegSites, *expectedFreq, *dpt, value;
	unsigned long int foo, *ullpt, sitesPerLinkageGroup, lgCount, SFScounts[(PLOIDY * N)];
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
	nTrackedSitesInParents = (unsigned long int) (expectedSegSites + 0.5);
	
	
	// now to choose the sites that are variable and their frequencies
	// siteIndexes = (unsigned long int *) malloc( nSITES * sizeof(unsigned long int) );
	parentalTrackedSiteIndexes = (unsigned long int *) malloc( nTrackedSitesInParents * sizeof(unsigned long int) );
//	ullpt = siteIndexes;
//	for ( foo = 0; foo < nSITES; foo++ ) {
//		*ullpt = foo;
//		ullpt++;
//	}
	
	// use gsl ran choose to pick sites at random to be standing neutral variation
	// gsl_ran_choose( rngState, parentalTrackedSiteIndexes, nTrackedSitesInParents, siteIndexes, nSITES, sizeof(unsigned long int) );
	chooseMutationSites( parentalTrackedSiteIndexes, nTrackedSitesInParents );
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
	alleleCounts = (unsigned long int *) malloc( nSITES * sizeof(unsigned long int));
	memset( alleleCounts, 0, nSITES * sizeof(unsigned long int) );
	for ( foo = 0; foo < (PLOIDY * N); foo++ )
		*(SFScounts + foo) = 0;
	setUpInitialAlleleFrequencies(expectedFreq, SFScounts);
	
	// pick allele frequencies using the Ewens sampling formula from Wakeley's book.
	
	
	if ( VERBOSE )
		printf("\ntheta = %f, segSites = %lu\n", theta, nTrackedSitesInParents);
	
	
	// assign site types
	double siteProbs[nSITE_CLASSES];
	unsigned long int siteTypeCounts[nSITE_CLASSES];
	memset( siteTypeCounts, 0, (nSITE_CLASSES * sizeof(unsigned long int)) );
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
			printf("\t%lu", siteTypeCounts[i]);
		foo += siteTypeCounts[i];
	}
	if ( VERBOSE )
		printf("\n");
	if ( foo != nSITES ) {
		fprintf(stderr, "\nError in setUpGenome():\n\tsite total (%lu) != nSITES (%lu)\n\t*** Exiting ***\n\n", foo, nSITES);
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
					fprintf( siteDesignations, "%lu", foo );
					fprintf( sd2, "%lu", foo );
					firstOne = 0;
				}
				else {
					fprintf( siteDesignations, ",%lu", foo );
					fprintf( sd2, " %lu", foo );
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
	unsigned long int focalSiteIndex;
	int SiteClassCode;
	initialFreqs = fopen("InitialAlleleFreqs.csv", "w");
	fprintf(initialFreqs, "SiteIndex,LinkageGroup,DerivedAlleleCount,DerivedAlleleFreq,SelectionCoefficient,SiteClassCode,SiteClassName\n");
	for ( foo = 0; foo < nTrackedSitesInParents; foo++ ) {
		focalSiteIndex = *(parentalTrackedSiteIndexes + foo);
		
		if ( *(siteClassifications + focalSiteIndex) > SITE_CLASS_NEUTRAL )
			nSelectedSites++;
		
		
		fprintf(initialFreqs, "%lu,%i,%lu,%E,%E", focalSiteIndex, linkageGroupMembership[focalSiteIndex], alleleCounts[focalSiteIndex], alleleFrequencies[focalSiteIndex], selectionCoefficients[focalSiteIndex] );
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
			fprintf(initialFreqs, "%lu,%lu\n", foo, SFScounts[foo]);
		}
	}
	fclose(initialFreqs);
	
	
	
	
	free(expectedFreq);
	
	//printf("\nWarning: setUpGenome() not done yet.  Still need to set selection coefficients.\n");
	//exit(0);
}


void setUpInitialAlleleFrequencies(double *expectedFreq, unsigned long int *SFScounts)
{
	unsigned long int i, copies, focalSiteIndex;
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
	unsigned long int neededSize, focalSiteIndex;
	int i;
	unsigned long int j, k, index, counter, siteAlleleCount;
	unsigned long int choiceVector[(PLOIDY * N)], allelesToSwitch[(PLOIDY * N)];
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
			printf("\nError in setUpPopulations():\n\tsiteAlleleCount (= %lu) out of bounds.\n", siteAlleleCount);
			exit(-1);
		}
		
		// int gsl_ran_choose (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size)
		gsl_ran_choose( rngState, allelesToSwitch, siteAlleleCount, choiceVector, (PLOIDY * N), sizeof(unsigned long int) );
		
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
				fprintf(initgts, ",Locus%luCopy%i", k, i);
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
			if ( i == 0 && nPOPULATIONS > 2 )
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







