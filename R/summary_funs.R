#' @title
#' Summarizing Discrete FWER Results
#' 
#' @description
#' `summary` method for class `DiscreteFWER`.
#'
#' @param object       an object of class `DiscreteFWER`.
#' @param x            an object of class `summary.DiscreteFWER`.
#' @param max          numeric or `NULL`, specifying the maximal number of
#'                     *rows* of the p-value table to be printed. By default,
#'                     when `NULL`, `getOption("max.print")` is used.
#' @param ...          further arguments passed to or from other methods.
#'
#' @details
#' `summary.DiscreteFWER` objects contain all data of an `DiscreteFWER` object,
#' but also include an additional table which includes the raw p-values,
#' their indices, the respective critical values (if present), the adjusted
#' p-values (if present) and a logical column to indicate rejection. The table
#' is sorted in ascending order by the raw p-values.
#' 
#' `print.summary.DiscreteFWER` simply prints the same output as
#' `print.DiscreteFWER`, but also prints the p-value table.
#' 
#' @return
#' `summary.DiscreteFWER` computes and returns a list that includes all the
#' data of an input `DiscreteFWER` object, plus
#' \item{Table}{`data.frame`, sorted by the raw p-values, that contains the
#'              indices, the raw p-values themselves, their respective critical
#'              values (if present), their adjusted p-values (if present) and a
#'              logical column to indicate rejection.}
#' 
#' @template example
#' @examples
#' # d-Holm procedure without critical values; using test results object
#' DFWER_dep_sd_fast <- discrete_FWER(test_results)
#' summary(DFWER_dep_sd_fast)
#' 
#' # d-Bonferroni procedure with critical values; using test results object
#' DFWER_dep_crit <- discrete_FWER(test_results, single_step = TRUE,
#' critical_values = TRUE)
#' summary(DFWER_dep_crit)
#' 
#' @rdname summary.DiscreteFWER
#' @export
## S3 method for class 'DiscreteFWER'
summary.DiscreteFWER <- function(object, ...){
  if(!("DiscreteFWER" %in% class(object)))
    return(summary(object))
    
  # determine if selection as performed
  select <- exists('Select', object)
  if(select) m <- object$Select$Number
  
  # number of tests
  n <- length(object$Data$Raw_pvalues)
  # determine order of raw p-values
  o <- order(object$Data$Raw_pvalues)
  # ordered indices
  i <- seq_len(n)
  # determine for each p-value if its corresponding null hypothesis is rejected
  r <- i %in% object$Indices #if(!select) o %in% object$Indices else o %in% object$Select.Indices[object$Indices]
  
  # create summary table
  tab <- data.frame('Index' = i, 'P.value' = object$Data$Raw_pvalues)
  if(select) {
    tab$Selected <- i %in% object$Select$Indices #rep(c(TRUE, FALSE), c(m, n - m))
    tab$Scaled <- NA
    tab$Scaled[tab$Selected] <- object$Select$Scaled
  }
  tab <- tab[o, ]
  if(exists('Critical_values', object)) {
    if(select) {
      tab$Critical.value <- NA
      tab$Critical.value[seq_len(m)] <- object$Critical_values[seq_len(m)][order(order(tab$Scaled[seq_len(m)]))]
    } else tab$Critical.value <- object$Critical_values
  }
  tab$Adjusted <- object$Adjusted[o]
  tab <- data.frame(tab, 'Rejected' = r[o])
  
  # if row names are numbers, rearrange them to represent sort order
  if(all(rownames(tab) == tab$Index)) rownames(tab) <- i
  
  # return output object
  out <- c(object, list(Table = tab))
  class(out) <- "summary.DiscreteFWER"
  return(out)
}

#'@rdname summary.DiscreteFWER
#'@export
## S3 method for class 'summary.DiscreteFWER'
print.summary.DiscreteFWER <- function(x, max = NULL, ...){
  if(!("summary.DiscreteFWER" %in% class(x)))
    return(print(x))
  
  # print 'DiscreteFWER' part of the object
  print.DiscreteFWER(x)
  
  # rows to print: number of rejections + 5 (if not requested otherwise)
  max <- if(!is.null(max)) ncol(x$Table) * max else getOption("max.print")
  
  # print additional summary table
  print(x$Table, max = max, ...)
  
  cat("\n")
  invisible(x)
}
