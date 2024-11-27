#' @name discrete_FWER
#' 
#' @title
#' Discrete Family-wise Error Rate (FWER) Adaptation Procedures
#' 
#' @description
#' Apply a discrete FWER adaptation procedure, with or without computing the
#' critical values, to a set of p-values and their discrete support.
#' 
#' @templateVar test_results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar independence TRUE
#' @templateVar single_step TRUE
#' @templateVar critical_values TRUE
#' @templateVar select_threshold TRUE
#' @templateVar pCDFlist_indices TRUE
#' @templateVar triple_dots TRUE
#' @template param
#'  
#' @template details_crit
#' @details
#'  
#' Depending on the choices of `independence` and `single_step`, one of the
#' following procedures, is applied:
#' 
#' |                 | single-step |      stepwise      |
#' |:----------------|:-----------:|:------------------:|
#' | independent     |    Šidák    | Hochberg (step-up) |
#' | not independent |  Bonferroni |  Holm (step-down)  |
#' 
#' Each procedure is available by its own shortcut function:
#' 
#' |                 |    single-step   |    stepwise   |
#' |:----------------|:----------------:|:-------------:|
#' | independent     |     `DSidak()`   | `DHochberg()` |
#' | not independent |  `DBonferroni()` |   `DHolm()`   |
#' 
#' @template return
#' 
#' @seealso
#' [`DiscreteFWER`][DiscreteFWER-package], [`DBonferroni()`], [`DHolm()`],
#' [`DSidak()`], [`DHochberg()`]
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
#' # d-Holm without critical values; using test results object
#' DFWER_dep_sd_fast <- discrete_FWER(test_results)
#' summary(DFWER_dep_sd_fast)
#' 
#' # d-Holm with critical values; using extracted p-values and supports
#' DFWER_dep_sd_crit <- discrete_FWER(raw_pvalues, pCDFlist, 
#'                                    critical_values = TRUE)
#' summary(DFWER_dep_sd_crit)
#' 
#' # d-Bonferroni without critical values; using test results object
#' DFWER_dep_fast <- discrete_FWER(test_results, single_step = TRUE)
#' summary(DFWER_dep_fast)
#' 
#' # d-Bonferroni with critical values; using extracted p-values and supports
#' DFWER_dep_crit <- discrete_FWER(raw_pvalues, pCDFlist, single_step = TRUE,
#'                                 critical_values = TRUE)
#' summary(DFWER_dep_crit)
#' 
#' # d-Hochberg without critical values; using test results object
#' DFWER_ind_su_fast <- discrete_FWER(test_results, independence = TRUE)
#' summary(DFWER_ind_su_fast)
#' 
#' # d-Hochberg with critical values; using extracted p-values and supports
#' DFWER_ind_su_crit <- discrete_FWER(raw_pvalues, pCDFlist, 
#'                                    independence = TRUE,
#'                                    critical_values = TRUE)
#' summary(DFWER_ind_su_crit)
#' 
#' # d-Šidák without critical values; using extracted p-values and supports
#' DFWER_ind_fast <- discrete_FWER(raw_pvalues, pCDFlist,
#'                                 independence = TRUE,
#'                                 single_step = TRUE)
#' summary(DFWER_ind_fast)
#' 
#' # d-Šidák with critical values; using test results object
#' DFWER_ind_crit <- discrete_FWER(test_results, independence = TRUE,
#'                                 single_step = TRUE, 
#'                                 critical_values = TRUE)
#' summary(DFWER_ind_crit)
#' 
#' @export
discrete_FWER <- function(test_results, ...) UseMethod("discrete_FWER")

#' @rdname discrete_FWER
#' @importFrom checkmate assert_integerish assert_list assert_numeric qassert
#' @export
discrete_FWER.default <- function(
    test_results,
    pCDFlist,
    alpha            = 0.05,
    independence     = FALSE,
    single_step      = FALSE,
    critical_values  = FALSE,
    select_threshold = 1,
    pCDFlist_indices = NULL,
    ...
) {
  #----------------------------------------------------
  #       check arguments
  #----------------------------------------------------
  # test results (p-values)
  qassert(x = test_results, rules = "N+[0, 1]")
  n <- length(test_results)
  
  # list structure of p-value distributions
  assert_list(
    x = pCDFlist,
    types = "numeric",
    any.missing = FALSE,
    min.len = 1,
    max.len = n
  )
  # individual p-value distributions
  for(i in seq_along(pCDFlist)) {
    assert_numeric(
      x = pCDFlist[[i]],
      lower = 0,
      upper = 1,
      any.missing = FALSE,
      min.len = 1,
      sorted = TRUE
    )
    if(max(pCDFlist[[i]]) != 1)
      stop("Last value of each vector in 'pCDFlist' must be 1!")
  }
  m <- length(pCDFlist)
  
  # FWERlevel
  qassert(x = alpha, rules = "N1[0, 1]")
  
  # independence
  qassert(independence, "B1")
  
  # step-down or single-step
  qassert(single_step, "B1")
  
  # compute and return critical values?
  qassert(critical_values, "B1")
  
  # selection threshold
  qassert(x = select_threshold, rules = "N1(0, 1]")
  
  # list structure of indices
  assert_list(
    x = pCDFlist_indices,
    types = "numeric",
    any.missing = FALSE,
    len = m,
    unique = TRUE,
    null.ok = TRUE
  )
  # individual index vectors (if not NULL)
  if(is.null(pCDFlist_indices)) {
    if(n != m) {
      stop(
        paste(
          "If no indices for the p-value CDFs are provided, the lengths of",
          "'test_results' and 'pCDFlist' must be equal!"
        )
      )
    }
    pCDFlist_indices <- as.list(1:n)
    #pCDFlist_counts <- rep(1, n)
  } else {
    set <- 1L:n
    for(i in seq_along(pCDFlist_indices)) {
      pCDFlist_indices[[i]] <- assert_integerish(
        x = pCDFlist_indices[[i]],
        lower = 1,
        upper = n,
        any.missing = FALSE,
        min.len = 1,
        max.len = n,
        unique = TRUE,
        sorted = TRUE,
        coerce = TRUE
      )
      set <- setdiff(set, pCDFlist_indices[[i]])
    }
    if(length(set))
      stop("'pCDFlist_indices' must contain each p-value index exactly once!")
    #pCDFlist_counts <- sapply(pCDFlist_indices, length)
  }
  
  #----------------------------------------------------
  #       check and prepare p-values for processing
  #----------------------------------------------------
  pvec <- match_pvals(test_results, pCDFlist, pCDFlist_indices)
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete_fwer_int(
    pvec             = test_results,
    pCDFlist         = pCDFlist,
    pCDFlist_indices = pCDFlist_indices,
    alpha            = alpha,
    independence     = independence,
    single_step      = single_step,
    crit_consts      = critical_values,
    threshold        = select_threshold,
    data_name        = paste(
                         deparse(substitute(test_results)),
                         "and",
                         deparse(substitute(pCDFlist))
                       )
  )
  
  return(output)
}

#' @rdname discrete_FWER
#' @importFrom checkmate assert_r6 qassert
#' @export
discrete_FWER.DiscreteTestResults <- function(
    test_results,
    alpha            = 0.05,
    independence     = FALSE,
    single_step      = FALSE,
    critical_values  = FALSE,
    select_threshold = 1,
    ...
) {
  #----------------------------------------------------
  #       check arguments
  #----------------------------------------------------
  # discrete test results object
  assert_r6(
    x = test_results,
    classes = "DiscreteTestResults",
    public = c("get_pvalues", "get_pvalue_supports", "get_support_indices")
  )
  
  # FWERlevel
  qassert(x = alpha, rules = "N1[0, 1]")
  
  # independence
  qassert(independence, "B1")
  
  # step-down or single-step
  qassert(single_step, "B1")
  
  # compute and return critical values?
  qassert(critical_values, "B1")
  
  # selection threshold
  qassert(x = select_threshold, rules = "N1(0, 1]")
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete_fwer_int(
    pvec             = test_results$get_pvalues(),
    pCDFlist         = test_results$get_pvalue_supports(unique = TRUE),
    pCDFlist_indices = test_results$get_support_indices(),
    alpha            = alpha,
    independence     = independence,
    single_step      = single_step,
    crit_consts      = critical_values,
    threshold        = select_threshold,
    data_name        = deparse(substitute(test_results))
  )
  
  return(output)
}
