#' @name DSidak
#' 
#' @title
#' Discrete Šidák Procedure for Independent Tests
#' 
#' @description 
#' `DSidak()` is a wrapper function of [`discrete_FWER()`] for computing the 
#' discrete Šidák procedure for independent discrete tests. It simply passes its
#' arguments to [`discrete_FWER()`] with fixed `independence = TRUE` and
#' `single_step = TRUE`.
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
#' [`discrete_FWER()`], [`DHochberg()`], [`DBonferroni()`], [`DHolm()`]
#' 
#' @references
#' Döhler, S. (2010). Validation of credit default probabilities using
#'   multiple-testing procedures. *Journal of Risk Model Validation*, *4*(4),
#'   59-92. \doi{10.21314/JRMV.2010.062}
#' 
#' @template example
#' @examples
#' # d-Šidák without critical values; using extracted p-values and supports
#' DSidak_fast <- DSidak(raw_pvalues, pCDFlist)
#' summary(DSidak_fast)
#' 
#' # d-Šidák with critical values; using test results object
#' DSidak_crit <- DSidak(test_results, critical_values = TRUE)
#' summary(DSidak_crit)
#' 
#' @export
DSidak <- function(test_results, ...) UseMethod("DSidak")

#' @rdname DSidak
#' @export
DSidak.default <- function(
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
    independence     = TRUE,
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

#' @rdname DSidak
#' @export
DSidak.DiscreteTestResults <- function(
    test_results,
    alpha            = 0.05,
    critical_values  = FALSE,
    select_threshold = 1,
    ...
) {
  out <- discrete_FWER.DiscreteTestResults(
    test_results     = test_results,
    alpha            = alpha,
    independence     = TRUE,
    single_step      = TRUE,
    critical_values  = critical_values,
    select_threshold = select_threshold,
    ...
  )
  
  out$Data$Data_name <- deparse(substitute(test_results))
  
  return(out)
}
