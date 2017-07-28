#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fffits.h"

void calculatePopGenMetrics(void)
{
	long int i, j, alleleCountsByPopulation[(nPOPULATIONS * nTrackedSitesInParents)];
	unsigned long long int SFScountsByPopulation[(nPOPULATIONS * 2 * N)], nDivSites = 0, nPosSites = 0;
	unsigned long long int siteMasterIndex, nSegSites = 0, nNeutralSites = 0, nBGsites = 0, SFScounts[(PLOIDY * N)];
	int pop, *locpt, hap;
	short int *sipt, siteClass;
	long int nhere, count, dumIndex;
	double oneOver2N, FSTarray[nTrackedSitesInParents], globalFreq;
	
	oneOver2N = 1.0 / (((double) PLOIDY) * ((double) N));
	memset( &alleleCountsByPopulation[0], 0, (sizeof(long int) * nPOPULATIONS * nTrackedSitesInParents) );
	
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
			
			// global allele frequencies and SFS
			count = *(alleleCounts + siteMasterIndex);
			if ( count <= 0 ) {
				fprintf(stderr, "\nError in calculatePopGenMetrics():\n\tcount (%li) <= 0 for siteMasterIndex = %llu\n", count, siteMasterIndex);
				exit(-1);
			}
			*(SFScounts + count) += 1;
			
			// allele counts by population/deme/patch
			locpt = locations;
			for ( j = 0; j < N; j++ ) {
				sipt = gts + (PLOIDY * i) + (PLOIDY * j * nTrackedSitesInParents);
				pop = *locpt; // location of individual
				dumIndex = (pop * nTrackedSitesInParents) + i; // where in alleleCountsByPopulation to store
				for ( hap = 0; hap < PLOIDY; hap++ ) {
					if ( *sipt == ALLELE_CODE_DERIVED )
						*(alleleCountsByPopulation + dumIndex) += 1;
					sipt++;
				}
				// alleleCountsByPopulation array: sites in consecutive order ("rows")
				// and populations across "columns"
				
				locpt++; // increment location pointer
			}
			
			globalFreq = (((double) count) * oneOver2N);
			// FST:
			FSTarray[i] = calculateFST( i, alleleCountsByPopulation, count, globalFreq );
			
			// print allele count data to file:
			fprintf(dataFile_alleleFreqTS, "%li,%llu,%i,%li,%E,%E,%i", t, siteMasterIndex, *(linkageGroupMembership + siteMasterIndex), count, globalFreq, *(selectionCoefficients + siteMasterIndex), *(siteClassifications + siteMasterIndex) );
			for ( j = 0; j < nPOPULATIONS; j++ ) {
				dumIndex = (j * nTrackedSitesInParents) + i;
				fprintf(dataFile_alleleFreqTS, ",%li", ( *(alleleCountsByPopulation + dumIndex) ) );
			}
			fprintf(dataFile_alleleFreqTS, ",%E\n", FSTarray[i]);
			
		}
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
	
	// calculate other common summary stats: FST, dxy, Tajima's D
	// need allele counts by population
	
	//calculateDXY();
	//calculateTajimasD();
	//calculateLD();
	
}


double calculateFST( long int i, long int *alleleCountsByPopulation, long int totalAlleleCount, double globalFreq )
{
	double FSTthisLocus = 0.0, HT, localFreq;
	double HS, ndoub, localHet, bigNdoub;
	int pop;
	long int j, *acbppt;
	
	if ( PLOIDY != 2 ) {
		fprintf(stderr, "\nError: calculateFST() only works for diploids\n\tExiting ...\n");
		exit(-1);
	}
	
	HT = 2.0 * (1.0 - globalFreq) * globalFreq; // total heterozygosity
	HS = 0.0; // within heterozygosity; will be built as a sum
	bigNdoub = (double) N;
	
	acbppt = alleleCountsByPopulation + i; // pointer to allele count in first populations at focal site
	for ( pop = 0; pop < nPOPULATIONS; pop++ ) {
		ndoub = (double) abundances[pop]; // local abundance as double
		localFreq = ((double) *acbppt) / (2.0 * ndoub); // local "p"
		localHet = 2.0 * localFreq * (1.0 - localFreq);
		HS += localHet * ( ndoub / bigNdoub ); // contribution of deme weighted by population size
		acbppt += nTrackedSitesInParents; // advance pointer
		
		if ( TEST_MODE ) {
			if ( localFreq < 0.0 || localFreq > 1.0 || localHet < 0.0 || localHet > 0.5 ) {
				fprintf(stderr, "\nError in calculateFST(): Values off\n\tlocalFreq = %f, localHet = %f\n", localFreq, localHet);
				exit(-1);
			}
		}
	}
	
	FSTthisLocus = (HT - HS) / HT;
	if ( TEST_MODE ) {
		if ( FSTthisLocus < 0.0 || FSTthisLocus > 1.0 || HT < 0.0 || HT > 0.5 || HS > HT || HS < 0.0 || HS > 0.5 ) {
			fprintf(stderr, "\nError in calculateFST(): Values off\n\tFSTthisLocus = %f, HT = %f, HS = %f\n", FSTthisLocus, HT, HS);
			exit(-1);
		}
		
	}
	
	return FSTthisLocus;
}
