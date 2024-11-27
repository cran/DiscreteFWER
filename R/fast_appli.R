#' @title 
#' Direct Application of Multiple Testing Procedures to Dataset
#' 
#' @description
#' Apply one of the various FWER adaptation procedures, with or without
#' computing the critical constants, to a data set of 2x2 contingency tables
#' using statistical test functions from package
#' [`DiscreteTests`][DiscreteTests::DiscreteTests-package]. If necessary,
#' functions for pre-processing can be passed as well.
#' 
#' @templateVar dat TRUE
#' @templateVar test_fun TRUE
#' @templateVar test_args TRUE
#' @templateVar alpha TRUE
#' @templateVar independence TRUE
#' @templateVar single_step TRUE
#' @templateVar critical_values TRUE
#' @templateVar select_threshold TRUE
#' @templateVar preprocess_fun TRUE
#' @templateVar preprocess_args TRUE
#' @template param
#' 
#' @template return
#' 
#' @template example
#' @examples
#' DBonf <- direct_discrete_FWER(df, "fisher")
#' summary(DBonf)
#' 
#' DHolm <- direct_discrete_FWER(df, "fisher_test_pv", single_step = FALSE)
#' summary(DHolm)
#' 
#' DBonf_bin <- direct_discrete_FWER(X1 + X2, "binom_test_pv", 
#'                                   list(n = N1 + N2, p = 0.05))
#' summary(DBonf_bin)
#' 
#' DHolm_bin <- direct_discrete_FWER(X1 + X2, "binom", 
#'                                   list(n = N1 + N2, p = 0.05),
#'                                   single_step = TRUE)
#' summary(DHolm_bin)
#' 
#' @export
#' @importFrom DiscreteFDR generate.pvalues
direct_discrete_FWER <- function(
  dat,
  test_fun, 
  test_args        = NULL,
  alpha            = 0.05, 
  independence     = FALSE,
  single_step      = TRUE,
  critical_values  = FALSE,
  select_threshold = 1,
  preprocess_fun   = NULL, 
  preprocess_args  = NULL
) {
  out <- discrete_FWER.DiscreteTestResults(
    test_results = generate.pvalues(
      dat             = dat,
      test.fun        = test_fun,
      test.args       = test_args,
      preprocess.fun  = preprocess_fun,
      preprocess.args = preprocess_args
    ),
    alpha            = alpha,
    independence     = independence,
    single_step      = single_step,
    critical_values  = critical_values,
    select_threshold = select_threshold
  )
  
  out$Data$Data.name <- deparse(substitute(dat))
  
  return(out)
}
