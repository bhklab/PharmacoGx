#include <Rcpp.h>


//' QUICKSTOP significance testing for partial correlation
//'
//' This function will test whether the observed partial correlation is significant
//' at a level of req_alpha, doing up to MaxIter permutations. Currently, it
//' supports only grouping by discrete categories when calculating a partial correlation.
//' Currenlty, only does two sided tests.
//'
//' @param pin_x one of the two vectors to correlate.
//' @param pin_y the other vector to calculate
//' @param pobsCor the observed (partial) correlation between these varaiables
//' @param pGroupFactor an integer vector labeling group membership, to correct
//' for in the partial correlation. NEEDS TO BE ZERO BASED!
//' @param pGroupSize an integer vector of size length(unique(pGroupFactor)), counting
//' the number of members of each group (basically table(pGroupFactor)) as integer vector
//' @param pnumGroup how many groups are there (len(pGroupSize))
//' @param pMaxIter maximum number of iterations to do, as a REAL NUMBER
//' @param pn length of x and y, as a REAL NUMBER
//' @param preq_alpha the required alpha for significance
//' @param ptolerance_par the tolerance region for quickstop. Suggested to be 1/100th of req_alpha'
//' @param plog_decision_boundary log (base e) of 1/probability of incorrectly calling significance, as
//' per quickstop paper (used to determine the log-odds)
//' @param pseed A numeric vector of length 2, used to seed the internal xoroshiro128+ 1.0
//' random number generator. Note that currently, these values get modified per call, so pass in a copy
//' if you wish to keep a seed for running same simulation twice
//'
//' @return a double vector of length 4, entry 1 is either 0, 1 (for TRUE/FALSE) or NA_REAL_ for significance determination
//' NA_REAL_ is returned when the MaxIter were reached before a decision is made. Usually, this occurs when the real p value is close to, or
//' falls within the tolerance region of (req_alpha, req_alpha+tolerance_par). Entry 2 is the current p value estimate. entry 3 is the total
//' number of iterations performed. Entry 4 is the number of time a permuted value was larger in absolute value than the observed cor.
//'
//' @useDynLib PharmacoGx _PharmacoGx_partialCorQUICKSTOP
//'
//'
// [[Rcpp::export]]
extern "C" SEXP partialCorQUICKSTOP(SEXP pin_x,
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
               SEXP pseed);
