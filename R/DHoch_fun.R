#' @name DHochberg
#' 
#' @title
#' Discrete Hochberg Procedure
#' 
#' @description 
#' `DHochberg()` is a wrapper function of [`discrete_FWER()`] for computing the
#' discrete Hochberg step-up procedure for independent or positively correlated
#' discrete tests. It simply passes its arguments to [`discrete_FWER()`] with
#' fixed `independence = TRUE` and `single_step = FALSE`.
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
#' [`discrete_FWER()`], [`DSidak()`], [`DBonferroni()`], [`DHolm`]
#' 
#' @references
#' Zhu, Y., & Guo, W. (2019). Family-Wise Error Rate Controlling Procedures for
#'   Discrete Data. *Statistics in Biopharmaceutical Research*, *12*(1), 
#'   117-128. \doi{10.1080/19466315.2019.1654912}
#'  
#' @template example
#' @examples
#' # d-Hochberg without critical values; using test results object
#' DHoch_fast <- DHochberg(test_results)
#' summary(DHoch_fast)
#' 
#' # d-Hochberg with critical values; using extracted p-values and supports
#' DHoch_crit <- DHochberg(raw_pvalues, pCDFlist, critical_values = TRUE)
#' summary(DHoch_crit)
#' 
#' @export
DHochberg <- function(test_results, ...) UseMethod("DHochberg")

#' @rdname DHochberg
#' @export
DHochberg.default <- function(
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
    single_step      = FALSE,
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

#' @rdname DHochberg
#' @export
DHochberg.DiscreteTestResults <- function(
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
    single_step      = FALSE,
    critical_values  = critical_values,
    select_threshold = select_threshold,
    ...
  )
  
  out$Data$Data_name <- deparse(substitute(test_results))
  
  return(out)
}
