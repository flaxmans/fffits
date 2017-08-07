#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fffits.h"

// "main" function that calls others here:
void dataRecording(void)
{
	double alleleFreqsByPop[(nPOPULATIONS * nTrackedSitesInParents)], globalFreqs[nTrackedSitesInParents];
	
	writeAbundances();
	
	// fitness time series?
	
	calcAlleleCountsFreqsFSTandSFS( alleleFreqsByPop, globalFreqs );
	/* Note for Joint SFS: the raw data needed for that are in AlleleFreqTS,
	 so it isn't calculated and printed; it would be a ton of data/files due to combinatorics (like LD) */
	
	calculateAndPrintPi( alleleFreqsByPop, globalFreqs );
	calculateAndPrintDXY( alleleFreqsByPop );
	fprintf(dataFile_PiAndDXY, "\n");
}



void calcAlleleCountsFreqsFSTandSFS(double *alleleFreqsByPop, double *globalFreqs)
{
	long int i, j, alleleCountsByPopulation[(nPOPULATIONS * nTrackedSitesInParents)];
	unsigned long int SFScountsByPopulation[(nPOPULATIONS * 2 * N)], nDivSites = 0, nPosSites = 0;
	unsigned long int siteMasterIndex, nSegSites = 0, nNeutralSites = 0, nBGsites = 0, SFScounts[(PLOIDY * N)];
	short int *sipt, siteClass;
	long int nhere, count, dumIndex, globalCounts[nTrackedSitesInParents], twoN = 2 * N;
	double oneOver2N, FSTarray[nTrackedSitesInParents];
	double abundanceDoubs[nPOPULATIONS];
	
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
		count = *(alleleCounts + siteMasterIndex);
		// storage in arrays for later calcs:
		globalCounts[i] = count;
		globalFreqs[i] = (((double) count) * oneOver2N);
		// counts by population:
		calculateAlleleCountsByPop(i, alleleCountsByPopulation);
		
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
#ifdef DEBUG
			if ( count < 0 ) {
				fprintf(stderr, "\nError in calculatePopGenMetrics():\n\tcount (%li) < 0 for siteMasterIndex = %lu, siteStatus = %i\n", count, siteMasterIndex, *(sitesStatuses + siteMasterIndex));
				exit(-1);
			}
#endif
			// site frequency spectrum:
			*(SFScounts + count) += 1;
		}
	}
	
	// get allele frequencies from counts:
	calculateAlleleFreqsByPop(alleleCountsByPopulation, alleleFreqsByPop);
	
	// calculate FST:
	calculateFST( FSTarray, alleleFreqsByPop, globalFreqs );

#ifdef DEBUG
	// error checking:
	// check allele counts
	long int alleleSums[nTrackedSitesInParents];
	
	for ( i = 0; i < nTrackedSitesInParents; i++ ) {
		count = 0;
		for ( j = 0; j < nPOPULATIONS; j++ ) {
			count += *(alleleCountsByPopulation + (nTrackedSitesInParents * j) + i);
		}
		siteMasterIndex = *(parentalTrackedSiteIndexes + i);
		if ( count != ( *(alleleCounts + siteMasterIndex) ) ) {
			fprintf(stderr, "\nError in calculatePopGenMetrics():\n\tsum count (%li) by popn != value from alleleCounts (%lu); globalFreqs[i] = %E\n", count, *(alleleCounts + siteMasterIndex), globalFreqs[i]);
			fprintf(stderr, "\n\tN = %li, t = %li, i = %li, j = %li, siteMasterIndex = %lu\n", N, t, i, j, siteMasterIndex);
			fprintf(stderr, "\n\tsiteStatus = %i, nTrackedSitesInParents = %lu\n", *(sitesStatuses + siteMasterIndex), nTrackedSitesInParents);
			for ( j = 0; j < nPOPULATIONS; j++ ) {
				fprintf(stderr, "\tpop%li count = %li\n", j, *(alleleCountsByPopulation + (nTrackedSitesInParents * j) + i) );
			}
			exit(-1);
		}
	}
#endif
	
	
	// print calculations to files:
	// pointers to big arrays:
	unsigned long int *ullipt = parentalTrackedSiteIndexes; // master site index
	
	// print allele counts and FST:
	for ( i = 0; i < nTrackedSitesInParents; i++ ) {
		siteMasterIndex = *ullipt;
		if ( globalCounts[i] > 0 && globalCounts[i] < twoN ) {
			fprintf(dataFile_alleleFreqTS, "%li,%lu,%i,%li,%E,%E,%i", t, siteMasterIndex, *(linkageGroupMembership + siteMasterIndex), globalCounts[i], globalFreqs[i], *(selectionCoefficients + siteMasterIndex), *(siteClassifications + siteMasterIndex) );
			for ( j = 0; j < nPOPULATIONS; j++ ) {
				dumIndex = (j * nTrackedSitesInParents) + i;
				fprintf(dataFile_alleleFreqTS, ",%li", ( *(alleleCountsByPopulation + dumIndex) ) );
			}
			fprintf(dataFile_alleleFreqTS, ",%E\n", FSTarray[i]);
		}
		ullipt++;
	}
	
	// print segregating site counts
	fprintf(dataFile_segSiteTS, "%li,%lu,%lu,%lu,%lu,%lu\n", t, nSegSites, nNeutralSites, nBGsites, nPosSites, nDivSites);
	
	// print the SFS
	for ( i = 1; i < (PLOIDY * N); i++ ) {
		if ( *(SFScounts + i) > 0 && *(SFScounts + i) < twoN ) {
			fprintf(dataFile_SFS_TS, "%li,%li,%lu\n", t, i, *(SFScounts + i));
		}
	}
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


void calculateAndPrintPi(double *alleleFreqsByPop, double *globalFreqs)
{
	long int i, j;
	double mySums[nSITE_CLASSES], mySum, p, *afbppt;
	unsigned long int *ullipt;
	short int siteClass;
	double piWithin[nSITE_CLASSES][nPOPULATIONS];
	
	// make sure we are at 0.0 since these are built as sums:
	for ( i = 0; i < nSITE_CLASSES; i++ ) {
		mySums[i] = 0.0;
		for ( j = 0; j < nPOPULATIONS; j++ )
			piWithin[i][j] = 0.0;
	}
	
	// going to calculate global averages as sums which can
	// then later be divided by kb or divided by number of segregating sites
	
	// global from overall frequencies, which could be called "between" heterozygosity
	mySum = 0.0;
	ullipt = parentalTrackedSiteIndexes; // pointer to master site indexes
	for ( i = 0; i < nTrackedSitesInParents; i++ ) {
		
		p = globalFreqs[i];
		p = 2.0 * p * (1.0 - p);
		mySum += p;
		
		// see fffits.h for site class codes
		siteClass = *(siteClassifications + (*ullipt)); // site class; note this is an int
		mySums[siteClass] += p;
		
		// now population specific
		afbppt = alleleFreqsByPop + i;
		for ( j = 0; j < nPOPULATIONS; j++ ) {
			p = *afbppt; // frequency at this site in this population
			piWithin[siteClass][j] += (2.0 * p * (1.0 - p));
			
			afbppt += nTrackedSitesInParents;
		}
		
		ullipt++;
	}

#ifdef DEBUG
	// error checking:
	double dumsum = 0.0, diff;
	for ( i = 0; i < nSITE_CLASSES; i++ )
		dumsum += mySums[i];
	diff = fabs(dumsum - mySum);
	if ( diff > 0.0001 )
		fprintf(stderr, "\nWarning from calculateAndPrintPi():\n\tmySum = %E but global dumsum = %E\n", mySum, dumsum);
#endif
	
	// print pi data:
	// order of prints: PiGlobalAllSites, PiGlobalNeutral, PiGlobalBG, PiGlobalPOS, PiGlobalDIV
	fprintf(dataFile_PiAndDXY, "%li,%E", t, mySum);
	for ( i = 0; i < nSITE_CLASSES; i++ )
		fprintf(dataFile_PiAndDXY, ",%E", mySums[i]);
	// continued prints: all populations, and within each, the summed pi for each site class in the same order as global pi values, but for each population
	for ( j = 0; j < nPOPULATIONS; j++ )
		for ( i = 0; i < nSITE_CLASSES; i++ )
			fprintf(dataFile_PiAndDXY, ",%E", piWithin[i][j]);
	// fprintf(dataFile_PiAndDXY, "\n"); // newline print is in dataRecording() since multiple functions write to this file
}


void calculateAndPrintDXY(double *alleleFreqsByPop)
{
	long int i, j, pop1, pop2;
	double p1, p2, *afbppt1, *afbppt2;
	unsigned long int *ullipt;
	short int siteClass;
	int numCombos;
	numCombos = (nPOPULATIONS * (nPOPULATIONS - 1)) / 2; // # of population pairs to compare
	double dxySums[nSITE_CLASSES][numCombos];
	
	// make sure we are at 0.0 since these are built as sums:
	for ( i = 0; i < nSITE_CLASSES; i++ ) {
		for ( j = 0; j < numCombos; j++ )
			dxySums[i][j] = 0.0;
	}
	
	// going to calculate values as sums which can
	// then later be divided by kb or divided by number of segregating sites
	
	// global from overall frequencies, which could be called "between" heterozygosity
	ullipt = parentalTrackedSiteIndexes; // pointer to master site indexes
	for ( i = 0; i < nTrackedSitesInParents; i++ ) {
		
		siteClass = *(siteClassifications + (*ullipt)); // site class; note this is an int
		
		// now population specific
		afbppt1 = alleleFreqsByPop + i;
		// j will be a population combo counter
		j = 0;
		for ( pop1 = 0; pop1 < (nPOPULATIONS-1); pop1++ ) {
			p1 = *afbppt1; // frequency at this site in this population 1
			afbppt2 = afbppt1 + nTrackedSitesInParents; // pointer to pop2's frequency at this site
			for ( pop2 = (pop1 + 1); pop2 < nPOPULATIONS; pop2++ ) {
				p2 = *afbppt2;
				// now calculate expected differences:
				dxySums[siteClass][j] += (p1 * (1.0 - p2)) + ((1.0 - p1) * p2);
				afbppt2++; // increment pop2's pointer
				j++; // increment counter of which pair we are on
			}
			afbppt1 += nTrackedSitesInParents; // increment pop1's pointer
		}
		
		ullipt++;
	}
	
#ifdef DEBUG
	// error checking:
	if ( j != numCombos ) {
		fprintf(stderr, "\nError from calculateAndPrintDXY():\n\tj = %li but numCombos = %i\n", j, numCombos);
		exit(-1);
	}
#endif
	
	// print DXY data:
	// order of prints: cover each site class; see initialization.c for header order prints
	for ( i = 0; i < nSITE_CLASSES; i++ ) {
		for ( j = 0; j < numCombos; j++ ) {
			fprintf(dataFile_PiAndDXY, ",%E", dxySums[i][j]);
		}
	}
	
	
	
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
		if ( globalFreq <= 0.0 || globalFreq >= 1.0 )
			FSTarray[i] = NAN;
		else {
			HT = 2.0 * (1.0 - globalFreq) * globalFreq; // total expected heterozygosity
			HS = 0.0; // within heterozygosity; will be built as a sum
			
			afbppt = alleleFreqsByPop + i; // pointer to allele count in first populations at focal site
			for ( pop = 0; pop < nPOPULATIONS; pop++ ) {
				localFreq = *afbppt; // local "p"
				localHet = 2.0 * localFreq * (1.0 - localFreq);
				HS += localHet * weights[pop]; // contribution of deme weighted by population size
				afbppt += nTrackedSitesInParents; // advance pointer
				
#ifdef DEBUG
				// error checking:
				if ( localFreq < 0.0 || localFreq > 1.0 || localHet < 0.0 || localHet > 0.5 ) {
					fprintf(stderr, "\nError in calculateFST(): Values off\n\tlocalFreq = %f, localHet = %f\n", localFreq, localHet);
					exit(-1);
				}
#endif
				
			}
			
			FSTarray[i] = (HT - HS) / HT;
			
#ifdef DEBUG
			// error checking:
			if ( FSTarray[i] < 0.0 || FSTarray[i] > 1.0 || HT < 0.0 || HT > 0.5 || HS > HT || HS < 0.0 || HS > 0.5 ) {
				fprintf(stderr, "\nError in calculateFST(): Values off\n\tFSTthisLocus = %f, HT = %f, HS = %f\n", FSTarray[i], HT, HS);
				exit(-1);
			}
#endif
		}
	}
}



void writeAbundances(void)
{
	int i;
	
	fprintf(dataFile_abundances, "%li,%li", t, N);
	for ( i = 0; i < nPOPULATIONS; i++ ) {
		fprintf(dataFile_abundances, ",%li", abundances[i]);
	}
	fprintf(dataFile_abundances, "\n");
	
#ifdef DEBUG
	long int totalCount = 0;
	for ( i = 0; i < nPOPULATIONS; i++ )
		totalCount += abundances[i];
	if ( N != totalCount ) {
		fprintf(stderr, "\nError from writeAbundances():\n\t N (%li) != totalCount (%li)\n", N, totalCount);
		exit(-1);
	}
#endif
	
}
