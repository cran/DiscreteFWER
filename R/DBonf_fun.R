#' @name DBonferroni
#' 
#' @title
#' Discrete Bonferroni Procedure
#' 
#' @description 
#' `DBonferroni()` is a wrapper function of [`discrete_FWER()`] for computing
#' the discrete Bonferroni procedure for tests with an arbitrary dependency
#' structure. It simply passes its arguments to [`discrete_FWER()`] with fixed
#' `independence = FALSE` and `single_step = TRUE`.
#' 
#' @templateVar test_results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar critical_values TRUE
#' @templateVar select_threshold TRUE
#' @templateVar pCDFlist_indices TRUE
#' @templateVar triple_dots TRUE
#' @template param
#' 
#' @template details_crit
#' 
#' @template return
#' 
#' @seealso
#' [`discrete_FWER()`], [`DHolm()`], [`DSidak()`], [`DHochberg()`]
#' 
#' @references
#' Döhler, S. (2010). Validation of credit default probabilities using
#'   multiple-testing procedures. *Journal of Risk Model Validation*, *4*(4),
#'   59-92. \doi{10.21314/JRMV.2010.062}
#' 
#' Zhu, Y., & Guo, W. (2019). Family-Wise Error Rate Controlling Procedures for
#'   Discrete Data. *Statistics in Biopharmaceutical Research*, *12*(1), 
#'   117-128. \doi{10.1080/19466315.2019.1654912}
#'  
#' @template example
#' @examples
#' # d-Bonferroni without critical values; using extracted p-values and supports
#' DBonferroni_fast <- DBonferroni(raw_pvalues, pCDFlist)
#' summary(DBonferroni_fast)
#' 
#' # d-Bonferroni with critical values; using test results object
#' DBonferroni_crit <- DBonferroni(test_results, critical_values = TRUE)
#' summary(DBonferroni_crit)
#' 
#' @export
DBonferroni <- function(test_results, ...) UseMethod("DBonferroni")

#' @rdname DBonferroni
#' @export
DBonferroni.default <- function(
    test_results,
    pCDFlist,
    alpha            = 0.05,
    critical_values  = FALSE,
    select_threshold = 1,
    pCDFlist_indices = NULL,
    ...
) {
  out <- discrete_FWER.default(
    test_results     = test_results,
    pCDFlist         = pCDFlist,
    alpha            = alpha,
    independence     = FALSE,
    single_step      = TRUE,
    critical_values  = critical_values,
    select_threshold = select_threshold,
    pCDFlist_indices = pCDFlist_indices,
    ...
  )
  
  out$Data$Data_name <- paste(
    deparse(substitute(test_results)),
    "and",
    deparse(substitute(pCDFlist))
  )
  
  return(out)
}

#' @rdname DBonferroni
#' @export
DBonferroni.DiscreteTestResults <- function(
    test_results,
    alpha            = 0.05,
    critical_values  = FALSE,
    select_threshold = 1,
    ...
) {
  out <- discrete_FWER.DiscreteTestResults(
    test_results     = test_results,
    alpha            = alpha,
    independence     = FALSE,
    single_step      = TRUE,
    critical_values  = critical_values,
    select_threshold = select_threshold,
    ...
  )
  
  out$Data$Data_name <- deparse(substitute(test_results))
  
  return(out)
}
