#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl_rng.h>            // gnu scientific library //
#include <gsl_randist.h>        // gnu scientific library //
#include "fffits.h"

// functions to create and use an exp() lookup table for (hopefully) optimization purposes
// rather than using randExp which is costly.


// first function called once to make the lookup table
void makeExpLookupTable(void)
{
    // idea of this function is to create a table
    // that will serve as an inverse exp CDF
    double meanRecombDistance, recombRatePerBP, dist;
    unsigned long int maxDistance, nTableElements, i;
    unsigned long int increment = 1;
    
    meanRecombDistance = 1000.0 / RECOMBINATION_RATE_PER_KB;
    maxDistance = ((unsigned long int) meanRecombDistance * 100.0); // upper threshold far past
    if ( maxDistance > nSITES )
        maxDistance = nSITES;
    nTableElements = maxDistance;
    
    while ( nTableElements > MAX_EXP_TABLE_SIZE ) {
        increment = increment * 2;
        nTableElements = nTableElements / 2;
    }
    if ( increment > 1000 )
        fprintf(stdout, "\nWarning from makeExpLookupTable():\n\tincrement (= %lu) larger than 1000 nucleotides\n", increment);
    // note recomb distances are truly discrete because nucleotides are discrete
    // however, if the possibilities are huge, we might choose an
    // increment of distances between them that is greater than 1
    // so that our table isn't too huge

    expLookupTable = (double *) malloc( (nTableElements + 1) * sizeof(double) );
    // +1 is for indexing including zero
    recombRatePerBP = RECOMBINATION_RATE_PER_KB / 1000.0;
    for ( i = 0; i < nTableElements; i++ ) {
        dist = ((double) (i * increment)) + 0.5; // the 0.5 is to make the divisions halfway between integer values
        expLookupTable[i] = 1.0 - exp( -recombRatePerBP * dist );
    }
    expLookupTable[nTableElements] = 1.1; // ensure no running off the end
    
    // global assignments:
    // median of exponential is ln(2)/lambda or ln(2) * mean
    expLookupMedian = (unsigned long int) ( (log(2.0) * meanRecombDistance) / ((double) increment) );
    if (expLookupMedian >= nTableElements)
        expLookupMedian = nTableElements - 1;
    expMedianVal = expLookupTable[expLookupMedian];
    expIncrement = increment;
    halfIncr = expIncrement / 2;
    
#ifdef DEBUG
    fprintf(stdout, "\nexpLookupTable status:\n\tmeanRecombDistance = %E, maxDistance = %lu\n\tnTableElements = %lu, increment = %lu\n", meanRecombDistance, maxDistance, (nTableElements + 1), increment);
    fprintf(stdout, "\texpLookupMedian = %lu, expMedianVal = %E\n", expLookupMedian, expMedianVal);
    fprintf(stdout, "\nSome entries:\n\ti\ti*incr\texpLookupTable[i]\n");
    fprintf(stdout, "At the beginning:\n");
    for ( i = 0; i <= 10; i++ )
        fprintf(stdout, "\t%lu\t%lu\t%E\n", i, (increment * i), expLookupTable[i]);
    fprintf(stdout, "At the end:\n");
    for ( i = nTableElements-10; i <= nTableElements; i++ )
        fprintf(stdout, "\t%lu\t%lu\t%E\n", i, (increment * i), expLookupTable[i]);
    unsigned long int expected;
    expected = (unsigned long int) (meanRecombDistance / ((double) increment));
    if ( expected < (nTableElements - 5) ) {
        fprintf(stdout, "Around the mean:\n");
        for ( i = (expected - 5); i <= (expected + 5); i++ )
            fprintf(stdout, "\t%lu\t%lu\t%E\n", i, (increment * i), expLookupTable[i]);
    }
    else
        fprintf(stdout, "Expected table index (= %lu) > nTableElements (= %lu) - 5\n", expected, nTableElements);
    if ( expLookupMedian < (nTableElements - 5) ) {
    fprintf(stdout, "Around the median:\n");
    for ( i = (expLookupMedian - 5); i <= (expLookupMedian + 5); i++ )
        fprintf(stdout, "\t%lu\t%lu\t%E\n", i, (increment * i), expLookupTable[i]);
    }
    else
        fprintf(stdout, "Median table index (= %lu) > nTableElements (= %lu) - 5\n", expected, nTableElements);
#endif
    
}


unsigned long int randExpLookup(void)
{
    // function to generate exponentially distributed random variates specifically for
    // crossover locations
    unsigned long int i = expLookupMedian;
    double dum, *testVal = (expLookupTable + expLookupMedian);
    
    dum = gsl_rng_uniform( rngState ); // test value for doing reverse CDF
    
    // this algorithim starts at the median and looks out
    if ( dum < expMedianVal ) {
        while ( dum < *(--testVal) )
            i--;
    }
    else {
        do {
            i++;
            testVal++;
        } while ( dum > *testVal );
    }

#ifdef DEBUG
    if ( i < 1 || i > nSITES ) {
        fprintf(stderr, "\nError in randExpLookup(): \n\ti = %lu\n", i);
        exit(-1);
    }
    
//    static int maxPrintTimes = 10;
//    static int last_t = 1;
//    long int thisMethod, randExpMethod, mydiff, myexpincr = (long int) expIncrement;
//    static double myCounter = 0.0, mySum = 0.0;
//    thisMethod = ( (i * expIncrement) - halfIncr );
//    randExpMethod = (unsigned long int) ((log(1.0 - dum)) * (-1000.0 / RECOMBINATION_RATE_PER_KB));
//    mydiff = (thisMethod - randExpMethod);
//    if ( thisMethod < nSITES ) {
//        if ( ((thisMethod - randExpMethod) > myexpincr) || ((thisMethod - randExpMethod) < -myexpincr ) ) {
//            fprintf(stderr, "\nError in randExpLookup(): \n\tthisMethod = %li, randExpMethod = %li, diff = %li, expIncrement = %lu\n", thisMethod, randExpMethod, mydiff, myexpincr);
//            exit(-1);
//        }
//        myCounter++;
//        mySum += (double) mydiff;
//    }
//    if ( t % 10 == 1 && t <= 101 ) {
//        if ( t != last_t )  {
//            last_t = t;
//            maxPrintTimes = 10;
//            fprintf(stdout, "\nt = %lu\n", t);
//        }
//        if ( maxPrintTimes ) {
//            maxPrintTimes--;
//            fprintf(stdout, "\tdum = %f, i = %lu, *testVal = %f, this = %li, randExp() = %li, diff = %li, avgDev = %f\n", dum, i, *testVal, thisMethod, randExpMethod, mydiff, mySum/myCounter );
//        }
//        // with algorithms as written above, this check confirms that the average deviation between
//        // this method and the randExp() method is 0.5, which is exactly what we expect if this
//        // algorithm is working properly, because the randExp() results were simply truncated (floored)
//        // to the integer component
//    }
#endif
    
    
    return ( (i * expIncrement) - halfIncr );
}








