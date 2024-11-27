#' @name print.DiscreteFWER
#' 
#' @title
#' Printing discrete FWER results
#' 
#' @description
#' Prints the results of discrete FWER analysis, stored in a `DiscreteFWER`
#' S3 class object.
#' 
#' @param x     object of class `DiscreteFWER`.
#' @param ...   further arguments to be passed to or from other methods. They
#'              are ignored in this function.
#' 
#' @return
#' The respective input object is invisibly returned via `invisible(x)`. 
#' 
#' @template example
#' @examples
#' # d-Holm with critical values; using test results object
#' DHolm_crit <- DHolm(test_results, critical.values = TRUE)
#' # print results
#' print(DHolm_crit)
#' 
#' @method print DiscreteFWER
#' @importFrom stats p.adjust
#' @export
## S3 method for class 'DiscreteFWER'
print.DiscreteFWER <- function(x, ...){
  if(!any(c("DiscreteFWER", "summary.DiscreteFWER") %in% class(x)))
    return(print(x))
  
  # determine if selection was performed
  select <- exists('Select', x)
  if(select) m <- x$Select$Number
  
  # number of tests
  n <- length(x$Data$Raw_pvalues)
  # number of rejected null hypotheses
  k <- x$Num_rejected
  
  if(x$Data$Independence) {
    if(x$Data$Single_step) {
      k_orig <- sum(x$Data$Raw_pvalues <= 1 - exp(log(1 - x$Data$FWER_level)/n))
      orig <- "Sidak"
    } else {
      k_orig <- sum(p.adjust(x$Data$Raw_pvalues, "hochberg") <= x$Data$FWER_level)
      orig <- "Hochberg"
    }
  } else if(x$Data$Single_step) {
    k_orig <- sum(p.adjust(x$Data$Raw_pvalues, "bonferroni") <= x$Data$FWER_level)
    orig <- "Bonferroni"
  } else {
    k_orig <- sum(p.adjust(x$Data$Raw_pvalues, "holm") <= x$Data$FWER_level)
    orig <- "Holm"
  }
  
  # print title (i.e. algorithm)
  cat("\n")
  cat("\t", x$Data$Method, "\n")
  
  # print dataset name(s)
  cat("\n")
  cat("Data: ", x$Data$Data_name, "\n")
  
  # print short results overview
  if(!select) {
    cat("Number of tests =", n, "\n")
  } else {
    cat("Number of selected tests =", m, "out of", n, "\n")
    cat("Selection threshold =", x$Select$Threshold, "\n")
  }
    
  cat("Number of rejections =", k, "at global FWER level", x$Data$FWER_level, "\n")
  cat("(Original ", orig, " rejections = ", k_orig, ")\n", sep = "")
  
  if(k && !select) {
    cat("Largest rejected p value: ", max(x$Rejected), "\n")
  }
  
  cat("\n")
  invisible(x)
}
