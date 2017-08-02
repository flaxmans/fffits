#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fffits.h"

// "main" function that calls others here:
void dataRecording(void)
{
	writeAbundances();
	calculatePopGenMetrics();
}


void calculateAlleleCountsByPop(long int focalSite, long int *alleleCountsByPopulation)
{
	// allele counts by population/deme/patch
	int *locpt, pop, hap, j;
	short int *sipt;
	long int dumIndex;
	
	locpt = locations;
	for ( j = 0; j < N; j++ ) {
		sipt = gts + (PLOIDY * focalSite) + (PLOIDY * j * nTrackedSitesInParents);
		pop = *locpt; // location of individual
		dumIndex = (pop * nTrackedSitesInParents) + focalSite; // where in alleleCountsByPopulation to store
		for ( hap = 0; hap < PLOIDY; hap++ ) {
			if ( *sipt == ALLELE_CODE_DERIVED )
				*(alleleCountsByPopulation + dumIndex) += 1;
			sipt++;
		}
		// alleleCountsByPopulation array: sites in consecutive order ("rows")
		// and populations across "columns"
		
		locpt++; // increment location pointer
	}
}


void calculateAlleleFreqsByPop(long int *alleleCountsByPopulation, double *alleleFreqsByPop)
{
	long int locus, *lipt;
	int pop;
	double *dpt, popDoub, countDoub;
	
	dpt = alleleFreqsByPop; // pointer to frequency array
	lipt = alleleCountsByPopulation; // pointer to count array
	
	// loop over populations, site by site:
	for ( pop = 0; pop < nPOPULATIONS; pop++ ) {
		popDoub = 1.0 / ( 2.0 * ((double) abundances[pop]) );
		for ( locus = 0; locus < nTrackedSitesInParents; locus++ ) {
			countDoub = (double) *lipt;
			*dpt = countDoub * popDoub; // allele frequency calculated
			dpt++; // increment pointers
			lipt++;
		}
	}
}



// basic stats first: allele frequencies/counts, used for SFS and other stats:
void calculatePopGenMetrics(void)
{
	long int i, j, alleleCountsByPopulation[(nPOPULATIONS * nTrackedSitesInParents)];
	unsigned long long int SFScountsByPopulation[(nPOPULATIONS * 2 * N)], nDivSites = 0, nPosSites = 0;
	unsigned long long int siteMasterIndex, nSegSites = 0, nNeutralSites = 0, nBGsites = 0, SFScounts[(PLOIDY * N)];
	short int *sipt, siteClass;
	long int nhere, count, dumIndex, globalCounts[nTrackedSitesInParents];
	double oneOver2N, FSTarray[nTrackedSitesInParents], globalFreqs[nTrackedSitesInParents];
	double abundanceDoubs[nPOPULATIONS], alleleFreqsByPop[(nPOPULATIONS * nTrackedSitesInParents)];
	
	oneOver2N = 1.0 / (((double) PLOIDY) * ((double) N));
	memset( &alleleCountsByPopulation[0], 0, (sizeof(long int) * nPOPULATIONS * nTrackedSitesInParents) );
	memset( &alleleFreqsByPop[0], 0, (sizeof(double) * nPOPULATIONS * nTrackedSitesInParents) );
	memset( &globalCounts[0], 0, (sizeof(long int) * nTrackedSitesInParents) );
	memset( &globalFreqs[0], 0, (sizeof(double) * nTrackedSitesInParents) );
	
	for ( count = 0; count < (PLOIDY * N); count++ )
		SFScounts[count] = 0;
	
	for ( i = 0; i < nTrackedSitesInParents; i++ ) {
		
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
			
			// global allele frequencies
			count = *(alleleCounts + siteMasterIndex);
			if ( TEST_MODE ) {
				if ( count <= 0 ) {
					fprintf(stderr, "\nError in calculatePopGenMetrics():\n\tcount (%li) <= 0 for siteMasterIndex = %llu\n", count, siteMasterIndex);
					exit(-1);
				}
			}
			// storage in arrays for later calcs:
			globalCounts[i] = count;
			globalFreqs[i] = (((double) count) * oneOver2N);
			// counts by population:
			calculateAlleleCountsByPop(i, alleleCountsByPopulation);
			// site frequency spectrum:
			*(SFScounts + count) += 1;
			
		}
	}
	
	// get allele frequencies from counts:
	calculateAlleleFreqsByPop(alleleCountsByPopulation, alleleFreqsByPop);
	
	// calculate FST:
	calculateFST( FSTarray, alleleFreqsByPop, globalFreqs );
	
	// print allele count data to file:
	// pointers to big arrays:
	unsigned long long int *ullipt = parentalTrackedSiteIndexes; // master site index
	
	for ( i = 0; i < nTrackedSitesInParents; i++ ) {
		siteMasterIndex = *ullipt;
		if ( globalCounts[i] > 0 ) {
			fprintf(dataFile_alleleFreqTS, "%li,%llu,%i,%li,%E,%E,%i", t, siteMasterIndex, *(linkageGroupMembership + siteMasterIndex), globalCounts[i], globalFreqs[i], *(selectionCoefficients + siteMasterIndex), *(siteClassifications + siteMasterIndex) );
			for ( j = 0; j < nPOPULATIONS; j++ ) {
				dumIndex = (j * nTrackedSitesInParents) + i;
				fprintf(dataFile_alleleFreqTS, ",%li", ( *(alleleCountsByPopulation + dumIndex) ) );
			}
			fprintf(dataFile_alleleFreqTS, ",%E\n", FSTarray[i]);
		}
		ullipt++;
	}
	
	// error checking:
	if ( TEST_MODE ) {
		// check allele counts
		long int alleleSums[nTrackedSitesInParents];
		
		for ( i = 0; i < nTrackedSitesInParents; i++ ) {
			count = 0;
			for ( j = 0; j < nPOPULATIONS; j++ ) {
				count += *(alleleCountsByPopulation + (nTrackedSitesInParents * j) + i);
			}
			siteMasterIndex = *(parentalTrackedSiteIndexes + i);
			if ( count != ( *(alleleCounts + siteMasterIndex) ) ) {
				fprintf(stderr, "\nError in calculatePopGenMetrics():\n\tsum count (%li) by popn != value from alleleCounts (%llu)\n", count, *(alleleCounts + siteMasterIndex));
				exit(-1);
			}
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
	
	/* Note for Joint SFS: the raw data needed for that are in AlleleFreqTS */
	
	// calculate other common summary stats: dxy, Tajima's D
	
	//calculateDXY();
	//calculateTajimasD();
	//calculateLD();
	// sliding or non-overlapping windows?
	
}


void calculateFST( double *FSTarray, double *alleleFreqsByPop, double *globalFreqs )
{
	if ( PLOIDY != 2 ) {
		fprintf(stderr, "\nError: calculateFST() only works for diploids\n\tExiting ...\n");
		exit(-1);
	}
	
	// calculateFST( FSTarray, alleleFreqsByPop, globalCounts, globalFreqs )
	double HT, localFreq, globalFreq, weights[nPOPULATIONS];
	double HS, ndoub, localHet, bigNdoub;
	double *afbppt;
	int pop;
	long int i, j;
	
	bigNdoub = (double) N;
	for ( pop = 0; pop < nPOPULATIONS; pop++ ) {
		ndoub = ((double) abundances[pop]);
		weights[pop] = ndoub / bigNdoub;
	}
	
	
	for ( i = 0; i < nTrackedSitesInParents; i++ ) {
		globalFreq = *(globalFreqs + i);
		HT = 2.0 * (1.0 - globalFreq) * globalFreq; // total expected heterozygosity
		HS = 0.0; // within heterozygosity; will be built as a sum
		
		afbppt = alleleFreqsByPop + i; // pointer to allele count in first populations at focal site
		for ( pop = 0; pop < nPOPULATIONS; pop++ ) {
			localFreq = *afbppt; // local "p"
			localHet = 2.0 * localFreq * (1.0 - localFreq);
			HS += localHet * weights[pop]; // contribution of deme weighted by population size
			afbppt += nTrackedSitesInParents; // advance pointer
			
			// error checking:
			if ( TEST_MODE ) {
				if ( localFreq < 0.0 || localFreq > 1.0 || localHet < 0.0 || localHet > 0.5 ) {
					fprintf(stderr, "\nError in calculateFST(): Values off\n\tlocalFreq = %f, localHet = %f\n", localFreq, localHet);
					exit(-1);
				}
			}
			
		}
		
		FSTarray[i] = (HT - HS) / HT;
		
		// error checking:
		if ( TEST_MODE ) {
			if ( FSTarray[i] < 0.0 || FSTarray[i] > 1.0 || HT < 0.0 || HT > 0.5 || HS > HT || HS < 0.0 || HS > 0.5 ) {
				fprintf(stderr, "\nError in calculateFST(): Values off\n\tFSTthisLocus = %f, HT = %f, HS = %f\n", FSTarray[i], HT, HS);
				exit(-1);
			}
			
		}
	}
}


void writeAbundances(void)
{
	int i;
	
	fprintf(dataFile_abundances, "%li", t);
	for ( i = 0; i < nPOPULATIONS; i++ ) {
		fprintf(dataFile_abundances, ",%li", abundances[i]);
	}
	fprintf(dataFile_abundances, "\n");
}
