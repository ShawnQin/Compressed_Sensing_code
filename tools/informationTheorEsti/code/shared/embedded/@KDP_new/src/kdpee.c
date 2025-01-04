/**
    This file is part of "kdpee", kd-partitioning entropy estimator.
    (c) 2009 Dan Stowell and Queen Mary University of London
    All rights reserved.

    kdpee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    kdpee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with kdpee.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "../kdpee/kdpee.h"
//#include <math.h>
//#include <stdlib.h>
//#include <stdbool.h>
//typedef int bool;

#ifdef _MSC_VER
const double loge_2 =  0.693147181; //natural logarithm of 2
// Calculates log2 of a number. You need this function in case of using Microsoft (Microsoft does not provide log2...).
    double log2( double n ){
    return log(n) / loge_2; // log(n)/log(2) is log2(n)
}
#endif

#ifndef false
	#define false 0
	#define true 1
#endif

#define LN 100000.0  //define a large integer used in initialization
#define PRECISION 1e-8 //minimum bin intercal, small than which is not valid
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

floatval kdpee(const floatval **dimrefs, const int n, const int d, floatval *mins, 
								floatval *maxs, const floatval zcut, int *keys){
	int minlev = ceil(0.5 * log2(n)); // Min partitioning level
	int i;
	floatval result;

	for(i = 0; i < n; ++i){
		keys[i] = i; // initialise keys
	}
	result = kdpee_recurse(dimrefs, n, d, mins, maxs, zcut, keys,
						false, 0, 1./n, 0, n-1, minlev);
	return result;
}

floatval kdpee_recurse(const floatval **dimrefs, const int n, const int d, floatval *mins, 
									floatval *maxs, const floatval zcut, int *keys,
									bool mayTerminate, int curlev, const floatval n_rec,
									int minindex, int maxindex, int minlev
									){
	
	int dimno = curlev % d; // The dimension along which we're intending to split
	int thesize = 1+maxindex-minindex; // No of points in this subset
	
	// As well as returning the median, this PARTITIONS the data (keys) in-place
	floatval median = kdpee_hoareMedian(dimrefs[dimno], keys, minindex, maxindex);
	floatval zscore = LN;
    

	if(curlev == minlev){
		mayTerminate = true; // We have passed the lower termination depth
	}
	/*
	if(mayTerminate){
		zscore = (sqrt(thesize) * (median+median-mins[dimno]-maxs[dimno]) 
						/ 
				(maxs[dimno]-mins[dimno]));
		if(zscore<0.){
			zscore = 0. - zscore;
		}
	}
*/
	
	// revised on 12/23/2017, based on another criterior for uniformity
	
    floatval minVal = mins[dimno];
    floatval maxVal = maxs[dimno];
    
	if(mayTerminate){
        //return deviation(thesize,dimrefs[dimno],dimno,maxs,mins);
        zscore = deviation(thesize,dimrefs[dimno],minindex,minVal,maxVal,keys);
        if(zscore<0.){
            zscore = 0. - zscore;
        }
        
	}

	if(thesize<=1 || (mayTerminate && (zscore < zcut))){
		// allowed to terminate, and z-score doesn't disprove uniformity, 
		// so let's calc the negsurprisal!
		floatval frac = thesize * n_rec;
		floatval volume = fmax(maxs[0] - mins[0],PRECISION); //revised on 01/03/2018
		int i;
		for(i=1; i<d; ++i){
			volume *= fmax(maxs[i] - mins[i],PRECISION); //revised 01/03/2018
		}
		
		if(volume==0.){
			return 0.;
		}else{
			return log(volume / frac) * frac;
		}
	}else{
		// We need to partition and recurse
		floatval oldextremum;
		floatval left = LN;
		floatval right = LN;

		int newmaxindex, newminindex;
		if((thesize & 1) == 0){ // even # points
			newmaxindex = minindex + thesize/2 - 1;
			newminindex = minindex + thesize/2;
		}else{ // odd # points
			newmaxindex = maxindex - (thesize+1)/2;
			newminindex = minindex + (thesize-1)/2;
		}
		
		// Remember the outer extremum, replace with median, then recurse
		oldextremum = maxs[dimno];
		maxs[dimno] = median;
		left = kdpee_recurse(dimrefs, n, d, mins, maxs, zcut, keys,
								mayTerminate, curlev+1, n_rec,
								minindex, newmaxindex, minlev
								);
		// put the extremum back in place
		maxs[dimno] = oldextremum;
		
		// Remember the outer extremum, replace with median, then recurse
		oldextremum = mins[dimno];
		mins[dimno] = median;
		right = kdpee_recurse(dimrefs, n, d, mins, maxs, zcut, keys,
								mayTerminate, curlev+1, n_rec,
								newminindex, maxindex, minlev
								);
		// put the extremum back in place
		mins[dimno] = oldextremum;
		
		return left + right;
	}
}

// by rearranging the keys between keys[minindex] & keys[maxindex] (inclusive),
//   find the median value in that section of oneRow
floatval kdpee_hoareMedian(const floatval *oneRow, int *keys, int minindex, int maxindex){
	int num = 1 + maxindex - minindex;
	
	if((num & 1) == 1){ // odd
		return kdpee_hoareFind(oneRow, keys, minindex, maxindex, (num - 1) / 2);
	}else{ // even
		return (
			kdpee_hoareFind(oneRow, keys, minindex, maxindex, (num / 2))
				+
			kdpee_hoareFind(oneRow, keys, minindex, maxindex, (num / 2) - 1)
				) * 0.5;
	}
}

floatval kdpee_hoareFind(const floatval *oneRow, int *keys, int minindex, int maxindex, int findThis){
	int l, r, i, j;
	i = minindex;
	j = maxindex;
	findThis += minindex; // offset so that we're actually in the desired section
	while(i < j){
		kdpee_hoarePartition(oneRow, keys, minindex, maxindex, i, j, oneRow[keys[findThis]], &l, &r);
		if(r < findThis){
			// kth smallest is in right split
			i = l;
		}
		if(findThis < l){
			// kth smallest is in left split
			j = r;
		}
	}
	
	// The desired element is now in desired position, so return its value
	return oneRow[keys[findThis]];
}

void kdpee_hoarePartition(const floatval *oneRow, int *keys, int minindex, int maxindex, 
											int l0, int r0, floatval fulcrum, int* l, int* r){
	int tmp;
	*l = l0;
	*r = r0;
	while(*l <= *r){
		// left_scan
		while(*l <= maxindex && oneRow[keys[*l]] < fulcrum){
			++*l;
		}
		// right_scan
		while(*r >= minindex && oneRow[keys[*r]] > fulcrum){
			--*r;
		}
		// check and exchange (keys)
		if(*l <= *r){
			tmp = keys[*l];
			keys[*l] = keys[*r];
			keys[*r] = tmp;
			// then
			++*l;
			--*r;
		}
	}
}


//definition of function "deviation" for test of uniformity

floatval deviation(int thesize,const floatval *oneRow,int minindex,floatval minVal, floatval maxVal,int *keys)
{
    // if interval between minimum and maximum is too small
	if(maxVal-minVal <PRECISION){
		return 0.0;
	}
	else{
	int nbin = ceil(sqrt(thesize));  //number of bins
	floatval binWidth = (maxVal-minVal)/nbin; //bin width
    
    //floatval rowData[thesize] = dimrefs[dimno];
    int bin_idx;
	//int *counts = NULL;
	//counts = (int *)malloc(nbin * sizeof(int));
	//if(counts==NULL)
	//	printf("can't get memory!\n");
	//memset(counts,0,nbin*sizeof(int));
    //int *counts = (int *)malloc(nbin * sizeof(int)); // revised on 1/1/2018
	//printf("binss:%d\n",nbin);
	//printf("the value in counts: %d\n",counts[0]);
	int *counts = (int *) calloc(nbin,sizeof(int));
    
    //int counts[nbin] = { 0 };
    for (int i = 0; i < thesize; ++i){
		//printf("the subgroup size: %d\n",thesize);
        bin_idx = floor((oneRow[keys[i+minindex]] - minVal) / binWidth);
		//printf("the bin index: %d\n",bin_idx);
        if (bin_idx == nbin) {bin_idx = nbin - 1;}
        counts[bin_idx] += 1;
    }
    
    //probability
    //floatval *prob = (floatval *)malloc(nbin * sizeof(floatval));
	floatval *prob = (floatval *)calloc(nbin,sizeof(floatval));
	//printf("the first value of initial prob: %f\n",prob[0]);
    //floatval prob[nbin] = {0};
    floatval normVal = 0.;
    for (int i = 0; i<nbin; ++i){
        //floatval douCounts = counts[i];
        //prob[i] = douCounts/thesize;
		prob[i] = (double) counts[i]/thesize;
        normVal += prob[i]*prob[i];
    }
    normVal = sqrt(normVal);
    
    floatval zscore = (normVal*sqrt(nbin) - 1)/(sqrt(nbin) - 1);
    if(zscore<0.){
        zscore = 0. - zscore;
    }
    free(counts);
    free(prob);
    counts = NULL;
    prob = NULL;
    return zscore;
}
    
}
