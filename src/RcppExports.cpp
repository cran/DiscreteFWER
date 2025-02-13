// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// kernel_DFWER_singlestep_fast
NumericVector kernel_DFWER_singlestep_fast(const List& pCDFlist, const NumericVector& pvalues, const bool independence, const Nullable<IntegerVector>& pCDFcounts);
RcppExport SEXP _DiscreteFWER_kernel_DFWER_singlestep_fast(SEXP pCDFlistSEXP, SEXP pvaluesSEXP, SEXP independenceSEXP, SEXP pCDFcountsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type pCDFlist(pCDFlistSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pvalues(pvaluesSEXP);
    Rcpp::traits::input_parameter< const bool >::type independence(independenceSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector>& >::type pCDFcounts(pCDFcountsSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_DFWER_singlestep_fast(pCDFlist, pvalues, independence, pCDFcounts));
    return rcpp_result_gen;
END_RCPP
}
// kernel_DFWER_singlestep_crit
List kernel_DFWER_singlestep_crit(const List& pCDFlist, const NumericVector& support, const NumericVector& sorted_pv, const double alpha, const bool independence, const Nullable<IntegerVector>& pCDFcounts);
RcppExport SEXP _DiscreteFWER_kernel_DFWER_singlestep_crit(SEXP pCDFlistSEXP, SEXP supportSEXP, SEXP sorted_pvSEXP, SEXP alphaSEXP, SEXP independenceSEXP, SEXP pCDFcountsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type pCDFlist(pCDFlistSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type support(supportSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sorted_pv(sorted_pvSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const bool >::type independence(independenceSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector>& >::type pCDFcounts(pCDFcountsSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_DFWER_singlestep_crit(pCDFlist, support, sorted_pv, alpha, independence, pCDFcounts));
    return rcpp_result_gen;
END_RCPP
}
// kernel_DFWER_stepwise_fast
NumericVector kernel_DFWER_stepwise_fast(const List& pCDFlist, const NumericVector& sorted_pv, const bool independence, const Nullable<List>& pCDFindices);
RcppExport SEXP _DiscreteFWER_kernel_DFWER_stepwise_fast(SEXP pCDFlistSEXP, SEXP sorted_pvSEXP, SEXP independenceSEXP, SEXP pCDFindicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type pCDFlist(pCDFlistSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sorted_pv(sorted_pvSEXP);
    Rcpp::traits::input_parameter< const bool >::type independence(independenceSEXP);
    Rcpp::traits::input_parameter< const Nullable<List>& >::type pCDFindices(pCDFindicesSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_DFWER_stepwise_fast(pCDFlist, sorted_pv, independence, pCDFindices));
    return rcpp_result_gen;
END_RCPP
}
// kernel_DFWER_stepwise_crit
List kernel_DFWER_stepwise_crit(const List& pCDFlist, const NumericVector& support, const NumericVector& sorted_pv, const double alpha, const bool independence, const Nullable<List>& pCDFindices);
RcppExport SEXP _DiscreteFWER_kernel_DFWER_stepwise_crit(SEXP pCDFlistSEXP, SEXP supportSEXP, SEXP sorted_pvSEXP, SEXP alphaSEXP, SEXP independenceSEXP, SEXP pCDFindicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type pCDFlist(pCDFlistSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type support(supportSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sorted_pv(sorted_pvSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const bool >::type independence(independenceSEXP);
    Rcpp::traits::input_parameter< const Nullable<List>& >::type pCDFindices(pCDFindicesSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_DFWER_stepwise_crit(pCDFlist, support, sorted_pv, alpha, independence, pCDFindices));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DiscreteFWER_kernel_DFWER_singlestep_fast", (DL_FUNC) &_DiscreteFWER_kernel_DFWER_singlestep_fast, 4},
    {"_DiscreteFWER_kernel_DFWER_singlestep_crit", (DL_FUNC) &_DiscreteFWER_kernel_DFWER_singlestep_crit, 6},
    {"_DiscreteFWER_kernel_DFWER_stepwise_fast", (DL_FUNC) &_DiscreteFWER_kernel_DFWER_stepwise_fast, 4},
    {"_DiscreteFWER_kernel_DFWER_stepwise_crit", (DL_FUNC) &_DiscreteFWER_kernel_DFWER_stepwise_crit, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_DiscreteFWER(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
