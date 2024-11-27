discrete_fwer_int <- function(
  pvec,
  pCDFlist,
  pCDFlist_indices,
  alpha        = 0.05,
  independence = FALSE,
  single_step  = TRUE,
  crit_consts  = FALSE,
  threshold    = 1,
  data_name    = NULL
) {
  # original number of hypotheses
  n <- length(pvec)
  
  #--------------------------------------------
  #       prepare output object
  #--------------------------------------------
  input_data <- list()
  input_data$Method <- if(independence) {
    # independence, i.e. Sidak/Hochberg
    paste("Discrete", ifelse(single_step, "Sidak", "Hochberg"), "procedure")
  } else {
    # not independence, i.e. Bonferroni/Holm
    paste("Discrete", ifelse(single_step, "Bonferroni", "Holm"), "procedure")
  }
  input_data$Raw_pvalues <- pvec
  if(length(pCDFlist) == n) {
    input_data$pCDFlist <- pCDFlist
  } else {
    idx <- unlist(pCDFlist_indices)
    pCDFlist_counts <- sapply(pCDFlist_indices, length)
    input_data$pCDFlist <- rep(pCDFlist, pCDFlist_counts)[order(idx)]
  }
  input_data$FWER_level   <- alpha
  input_data$Independence <- independence
  input_data$Single_step  <- single_step
  input_data$Data_name    <- ifelse(
    !is.null(data_name),
    data_name,
    paste(deparse(substitute(pvec)), "and", deparse(substitute(pCDFlist)))
  )
  
  #--------------------------------------------
  #       apply p-value selection
  #--------------------------------------------
  if(threshold < 1) {
    # which p-values are not above threshold?
    select <- which(pvec <= threshold)
    # number of selected p-values
    m <- length(select)
    # filter pCDFs, indices and counts of selected p-values
    if(is.null(pCDFlist_indices)) {
      # keep CDFs of selected p-values
      pCDFlist         <- pCDFlist[select]
      # create indices according to selction
      pCDFlist_indices <- as.list(seq_len(m))
      # set counts (each CDF exists exactly once)
      pCDFlist_counts  <- rep(1, m)
    } else{
      # remove indices that were not selected
      pCDFlist_indices <- lapply(pCDFlist_indices, setdiff, seq_len(n)[-select])
      # determine counts
      pCDFlist_counts  <- sapply(pCDFlist_indices, length)
      # determine CDFs with non-zero counts
      idx_nonempty     <- which(pCDFlist_counts > 0)
      # keep CDFs with non-zero counts
      pCDFlist         <- pCDFlist[idx_nonempty]
      pCDFlist_indices <- pCDFlist_indices[idx_nonempty]
      pCDFlist_counts  <- pCDFlist_counts[idx_nonempty]
      # determine new indices (selection changes numbering!)
      new_idx          <- rep(NA, n)
      new_idx[select]  <- seq_len(m)
      # change indices according to selection
      pCDFlist_indices <- lapply(pCDFlist_indices, function(l) new_idx[l])
    }
    pCDFlist_idx <- order(unlist(pCDFlist_indices))
    # rescale pCDFs
    F_thresh <- sapply(pCDFlist, function(X) {t <- which(X <= threshold); if(length(t)) X[max(t)] else 0})
    pCDFlist <- sapply(seq_along(pCDFlist), function(k) pCDFlist[[k]] / F_thresh[k])
    pCDFlist <- lapply(pCDFlist, function(f) f[f <= 1])
    # rescale selected p-values
    pvec <- pvec[select] / rep(F_thresh, pCDFlist_counts)[pCDFlist_idx]
  } else {
    # all p-values were selected
    select <- seq_len(n)
    m <- n
    # use original counts (or 1 for all, if all pCDFs are unique)
    pCDFlist_counts <- if(is.null(pCDFlist_indices))
      rep(1, m) else
        sapply(pCDFlist_indices, length)
    # F_i(1) = 1 for all i = 1, ..., n
    F_thresh <- rep(1.0, n)
  }
  #--------------------------------------------
  #       determine sort order and do sorting
  #--------------------------------------------
  ord <- order(pvec)
  org_ord <- order(ord)
  sorted_pvals <- pvec[ord]
  sorted_pCDFlist_indices <- if(!is.null(pCDFlist_indices))
    lapply(pCDFlist_indices, function(l) sort(org_ord[l])) else
      as.list(org_ord)
  
  #--------------------------------------------
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  support <- unique(sort(pmin(as.numeric(unlist(pCDFlist)), 1.0)))
  
  #--------------------------------------------
  #        compute significant p-values, their
  #        indices and the number of rejections
  #--------------------------------------------
  if(crit_consts) {
    if(single_step) {
      res <- kernel_DFWER_singlestep_crit(
        pCDFlist, support, sorted_pvals, alpha, independence, pCDFlist_counts
      )
      crit_constants <- res$crit_consts
      idx_rej <- which(sorted_pvals <= crit_constants)
    } else {
      res <- kernel_DFWER_stepwise_crit(
        pCDFlist, support, sorted_pvals, alpha, independence, sorted_pCDFlist_indices
      )
      crit_constants <- res$crit_consts
      idx_rej <- if(independence) 
        which(sorted_pvals <= crit_constants) else
          which(sorted_pvals > crit_constants)
    }
  } else {
    if(single_step) {
      res <- kernel_DFWER_singlestep_fast(
        pCDFlist, sorted_pvals, independence, pCDFlist_counts
      )
      idx_rej <- which(res <= alpha)
    } else {
      res <- kernel_DFWER_stepwise_fast(
        pCDFlist, sorted_pvals, independence, sorted_pCDFlist_indices
      )
      idx_rej <- if(independence) 
        which(res <= alpha) else
          which(res > alpha)
    }
  }
  
  k <- length(idx_rej)
  if(single_step || (!single_step && independence)) {
    if(k > 0) {
      m_rej <- max(idx_rej)
      # determine significant (observed) p-values in sorted_pvals
      idx_rej <- which(pvec <= sorted_pvals[m_rej]) 
      pvec_rej <- input_data$Raw_pvalues[select][idx_rej]
    } else {
      m_rej <- 0
      idx_rej <- integer(0)
      pvec_rej <- numeric(0)
    }
  } else {
    if(k > 0) {
      m_rej <- min(idx_rej) - 1
      if(m_rej) {
        # determine significant (observed) p-values in sorted_pvals
        idx_rej <- which(pvec <= sorted_pvals[m_rej])
        pvec_rej <- input_data$Raw_pvalues[select][idx_rej]
      } else {
        idx_rej <- numeric(0)
        pvec_rej <- numeric(0)
      }
    } else {
      m_rej <- m
      idx_rej <- seq_len(m)
      pvec_rej <- input_data$Raw_pvalues[select]
    }
  }
  
  #--------------------------------------------
  #       create output object
  #--------------------------------------------
  # rejections
  output <- list(
    Rejected = pvec_rej,
    Indices = select[idx_rej],
    Num_rejected = m_rej
  )
  
  # adjusted p-values
  pv_adj <- if(crit_consts) res$pval_transf else res
  # add adjusted p-values to output list
  output$Adjusted          <- numeric(n)
  output$Adjusted[select]  <- pv_adj[org_ord]
  output$Adjusted[-select] <- NA
    
  # add critical values to output list
  if(crit_consts) {
    output$Critical_values          <- numeric(n)
    output$Critical_values[select]  <- crit_constants
    output$Critical_values[-select] <- NA
  }
  
  # original test data
  output$Data <- input_data
  
  # include selection data, if selection was applied
  if(threshold < 1) {
    output$Select <- list()
    output$Select$Threshold <- threshold
    output$Select$Effective_Thresholds <- F_thresh
    output$Select$Pvalues <- input_data$Raw_pvalues[select]
    output$Select$Indices <- select
    output$Select$Scaled <- pvec
    output$Select$Number <- m
  }
  
  class(output) <- "DiscreteFWER"
  return(output)
}
