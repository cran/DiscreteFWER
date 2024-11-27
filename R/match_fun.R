#' @name match_pvals
#' 
#' @title
#' Matching Raw P-Values with Supports
#' 
#' @keywords internal
#' 
#' @description 
#' Constructs the observed p-values from the raw observed p-values, by rounding
#' them to their nearest neighbour, matching with the supports of their
#' respective CDFs (as in function `p.discrete.adjust()` of package
#' `discreteMTP`, which is no longer available on CRAN).
#' 
#' **Note**: This is an internal function and has to be called directly via
#' `:::`, i.e. `DiscreteFWER:::match_pvals()`.
#' 
#' @details
#' Well computed raw p-values should already belong to their respective CDF
#' support. So this function is called at the beginning of
#' [`discrete_FWER.default()`] and its wrappers, just in case raw p-values may
#' be biased.
#'
#' For each raw p-value that needs to be rounded, a warning is issued.
#'
#' @seealso
#' [`discrete_FWER()`]
#'
#' @templateVar test_results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar pCDFlist_indices TRUE
#' @template param
#' 
#' @return
#' A vector where each raw p-value has been replaced by its nearest neighbour,
#' if necessary.
#'
match_pvals <- function(test_results, pCDFlist, pCDFlist_indices = NULL) {
  m <- length(test_results)
  if(!is.null(pCDFlist_indices)) {
    idx <- unlist(pCDFlist_indices)
    counts <- sapply(pCDFlist_indices, length)
    pCDFlist <- rep(pCDFlist, counts)[order(idx)]
  }
  n <- length(pCDFlist)
  if(m > 0 && m == n) {
    pvec <- test_results
    in.CDF <- numeric(m)
    for(k in seq_len(m)) {
      in.CDF[k] <- match(pvec[k], pCDFlist[[k]])
      if(is.na(in.CDF[k])) {
        in.CDF[k] <- which.min(abs(pCDFlist[[k]] - pvec[k]))
        pvec[k] <- pCDFlist[[k]][in.CDF[k]]
        ordinal <- "th"
        if(k %% 10 == 1) ordinal <- "st"
        if(k %% 10 == 2) ordinal <- "nd"
        if(k %% 10 == 3) ordinal <- "rd"
        if(k %% 100 - k %% 10 == 10) ordinal <- "th"
        warning("Since ", test_results[k], 
                " is not a value of the CDF of the ", k, ordinal, " p-value,\n",
                "  the p-value is rounded to be ", pCDFlist[[k]][in.CDF[k]],
                call. = F)
      }
    }
    return(pvec)
  } else {
    stop("'pCDFlist' and 'test_results' do not match")
  }
}