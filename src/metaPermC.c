/*
// Fast permutations for rCI using a naive matrix based approach.
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

// #include "xoroshiro128+.h"


/* This is xoroshiro128+ 1.0, our best and fastest small-state generator
   for floating-point numbers. We suggest to use its upper bits for
   floating-point generation, as it is slightly faster than
   xoroshiro128**. It passes all tests we are aware of except for the four
   lower bits, which might fail linearity tests (and just those), so if
   low linear complexity is not considered an issue (as it is usually the
   case) it can be used to generate 64-bit outputs, too; moreover, this
   generator has a very mild Hamming-weight dependency making our test
   (http://prng.di.unimi.it/hwd.php) fail after 5 TB of output; we believe
   this slight bias cannot affect any application. If you are concerned,
   use xoroshiro128** or xoshiro256+.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s.

   NOTE: the parameters (a=24, b=16, b=37) of this version give slightly
   better results in our test than the 2016 version (a=55, b=14, c=36).
*/

static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}


// static uint64_t s[2];

uint64_t next(uint64_t *s) {
  const uint64_t s0 = s[0];
  uint64_t s1 = s[1];
  const uint64_t result = s0 + s1;

  s1 ^= s0;
  s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
  s[1] = rotl(s1, 37); // c

  return result;
}

// Code to trim random numbers from :https://stackoverflow.com/questions/822323/how-to-generate-a-random-int-in-c
// maxindex here is exclusive of the right edge
uint64_t generate_random_index(uint64_t *state, uint64_t maxindex){
   // uint64_t maxindex = llround(maxindexd);

   if ((maxindex-1) == UINT64_MAX){
      return next(state);
   } else {
    // Supporting larger values for n would requires an even more
    // elaborate implementation that combines multiple calls to rand()
    // assert (maxindex <= UINT64_MAX)

    // Chop off all of the values that would cause skew...
    uint64_t end = UINT64_MAX / (maxindex); // truncate skew
    // assert (end > 0);
    end *= maxindex;

    // ... and ignore results from rand() that fall above that limit.
    // (Worst case the loop condition should succeed 50% of the time,
    // so we can expect to bail out of this loop pretty quickly.)
    uint64_t r;
    while ((r = next(state)) >= end);

    return r % maxindex;
  }
}



// void printVec(uint64_t *list, uint64_t N){
//   for(uint64_t i = 0; i < N; i ++){
//     printf("%lld, ", list[i]);
//   }
// }

// Using "inside out" Fisher Yates https://en.wikipedia.org/wiki/Fisher–Yates_shuffle
void sampleIdx(uint64_t N, uint64_t *permPointer, uint64_t *state){
  uint64_t j;
  permPointer[0] = 0;
  // printf("%lld", N);

  for(uint64_t i = 0; i <= N-1; i++){
      j = generate_random_index(state, i+1);
      if(j != i){
        permPointer[i] = permPointer[j];
      }
      permPointer[j] = i;
  }

}


// This code is adapted from the quickstop reference implementation here: https://github.com/julianhecker/QUICK-STOP

double log_denom(uint64_t suc,uint64_t n, double p)
{
   double tmp=0;
   if(fabs(p)>pow(10,-16))
   {
    tmp+=(double)(suc)*log(p);
   }
   if(fabs(1.0-p)>pow(10,-16))
   {
    tmp+=((double)(n)-(double)(suc))*log(1.0-p);
   }
   return tmp;
}



void runPerm(double *out,
               double *xvec,
               double *yvec,
               double obsCor,
               int *GroupFactor,
               uint64_t N,
               int *GroupSize,
               int numGroup,
               double req_alpha,
               double tolerance_par,
               int log_decision_boundary,
               uint64_t max_iter,
               uint64_t *state){



  uint64_t num_larger = 0;
  uint64_t cur_success = 0;
  uint64_t cur_iter = 1;

  double log_cur_PiN = log(1); // Initialization so the logic can stay the same through all loops
  double log_cur_suph1;
  double log_cur_suph2;

  double p1 = req_alpha;
  double p2 = p1 + tolerance_par;

  double pr_min_1 = (double)1.0/2;
  // int success;


  double currCor;

  // uint64_t j = 0;
  uint64_t *permIdxX = malloc(N * sizeof(uint64_t));
  uint64_t *permIdxY = malloc(N * sizeof(uint64_t));

  double dSDx;
  double dSDy;

  double *dGroupMeanY = malloc(numGroup*sizeof(double));
  double *dGroupMeanX = malloc(numGroup*sizeof(double));

  double *yShuffled = malloc(N * sizeof(double));
  double *xShuffled = malloc(N * sizeof(double));
  double significant = NA_REAL;



          // sample_function <- function(){

          //   partial.dp <- sample(dd[,1], nrow(dd))
          //   partial.x <- sample(dd[,2], nrow(dd))

          //   for(gp in unique(dd[,3])){
          //     partial.x[dd[,3]==gp] <- partial.x[dd[,3]==gp]-mean(partial.x[dd[,3]==gp])
          //     partial.dp[dd[,3]==gp] <- partial.dp[dd[,3]==gp]-mean(partial.dp[dd[,3]==gp])
          //   }

          //   perm.cor <- coop::pcor(partial.dp, partial.x, use="complete.obs")
          //   return(abs(obs.cor) < abs(perm.cor))

  while(cur_iter <= max_iter){

    sampleIdx(N, permIdxX, state);
    sampleIdx(N, permIdxY, state);
    dSDx = (double)0;
    dSDy = (double)0;

    for(int j = 0; j < numGroup; j++){
      dGroupMeanX[j] = 0;
      dGroupMeanY[j] = 0;
    }

    for(int j = 0; j < N; j++){
      yShuffled[j] = yvec[permIdxY[j]];
      dGroupMeanY[GroupFactor[j]] += yvec[permIdxY[j]];

      xShuffled[j] = xvec[permIdxX[j]];
      dGroupMeanX[GroupFactor[j]] += xvec[permIdxX[j]];

    }

    // for(int j = 0; j < N; j++){
    //   yShuffled[j] = yvec[j];
    //   dGroupMeanY[GroupFactor[j]] += yvec[j];

    //   xShuffled[j] = xvec[j];
    //   dGroupMeanX[GroupFactor[j]] += xvec[j];

    // }

    for(int j = 0; j < numGroup; j++){
      dGroupMeanX[j] /= (double)GroupSize[j];
      dGroupMeanY[j] /= (double)GroupSize[j];
    }

    for(int j = 0; j < N; j++){
      yShuffled[j] = yShuffled[j] - dGroupMeanY[GroupFactor[j]];
      xShuffled[j] = xShuffled[j] - dGroupMeanX[GroupFactor[j]];

      dSDy += pow(yShuffled[j],2);

      dSDx += pow(xShuffled[j],2);

    }


    dSDy = sqrt(dSDy);
    dSDx = sqrt(dSDx);

    currCor = 0;

    for(int j = 0; j < N; j++){
      currCor += (yShuffled[j])/dSDy * (xShuffled[j])/dSDx;
    }
    // printf("%f\n", currCor);



    if(fabs(currCor)>=fabs(obsCor)){
      cur_success = cur_success + 1;
      log_cur_PiN = log_cur_PiN + log_denom(1,1,pr_min_1);
    } else {
      log_cur_PiN = log_cur_PiN + log_denom(0,1,pr_min_1);
    }

    if(pr_min_1<p1) {
      log_cur_suph1 = log_denom(cur_success, cur_iter, pr_min_1);
    } else {
      log_cur_suph1 = log_denom(cur_success, cur_iter, p1);
    }

    if(pr_min_1>p2) {
      log_cur_suph2 = log_denom(cur_success, cur_iter, pr_min_1);
    } else {
      log_cur_suph2 = log_denom(cur_success, cur_iter, p2);
    }

    cur_iter = cur_iter + 1;
    pr_min_1 = ((double)cur_success + 1.0/2.0)/cur_iter;
    // printf("%f\n", pr_min_1);


    // TODO: make a list to return everything to R
    if(log_cur_PiN - log_cur_suph2 > log_decision_boundary){
      significant = (double) 1;
      break;
    }
    if(log_cur_PiN - log_cur_suph1 > log_decision_boundary){
      significant = (double) 0;
      break;
    }

  }



  free(permIdxX);
  free(permIdxY);

  free(dGroupMeanY);
  free(dGroupMeanX);
  free(yShuffled);
  free(xShuffled);


  out[0] = significant;
  out[1] = pr_min_1;
  out[2] = (double) cur_iter;
  out[3] = (double) cur_success;

  return;
  // return(currCor);

}

// Tested the computation of correlations by comparing against R code, perm distributions look
// very similar
SEXP partialCorQUICKSTOP(SEXP pin_x,
               SEXP pin_y,
               SEXP pobsCor,
               SEXP pGroupFactor,
               SEXP pGroupSize,
               SEXP pnumGroup,
               SEXP pMaxIter,
               SEXP pn,
               SEXP preq_alpha,
               SEXP ptolerance_par,
               SEXP plog_decision_boundary,
               SEXP pseed){

  double Ndouble = *REAL(pn);

  double MaxIterdouble = *REAL(pMaxIter);
  double obsCor = *REAL(pobsCor);
  double req_alpha = *REAL(preq_alpha);
  double tolerance_par = *REAL(ptolerance_par);
  int log_decision_boundary = *INTEGER(plog_decision_boundary);


  double temp;

  uint64_t N = (uint64_t) Ndouble;
  uint64_t max_iter = (uint64_t) MaxIterdouble;

  SEXP pout = PROTECT(allocVector(REALSXP,4));

  double *out = REAL(pout);

  double *seed = REAL(pseed);
  uint64_t *state = (uint64_t*) seed;

  int numGroup = *INTEGER(pnumGroup);

  runPerm(out, REAL(pin_x), REAL(pin_y),
                  obsCor,
                  INTEGER(pGroupFactor),
                  N,
                  INTEGER(pGroupSize),
                  numGroup,
                  req_alpha,
                  tolerance_par,
                  log_decision_boundary,
                  max_iter,
                  state);

  UNPROTECT(1);

  return pout;

}
