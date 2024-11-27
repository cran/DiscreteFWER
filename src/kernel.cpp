#include "kernel.h"

NumericVector kernel_DFWER_singlestep_fast(
  const List& pCDFlist,
  const NumericVector& pvalues,
  const bool independence,
  const Nullable<IntegerVector>& pCDFcounts
) {
  // Number of p-values
  int numValues = pvalues.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // counts of the CDFs
  IntegerVector CDFcounts;
  if(numCDF == numValues || pCDFcounts.isNull() || as<IntegerVector>(pCDFcounts).length() == 0) 
    CDFcounts = IntegerVector(numCDF, 1.0);
  else 
    CDFcounts = pCDFcounts;
  
  // extract p-value CDF vectors
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  // vector to store transformed p-values
  NumericVector pval_transf(numValues);
  // evaluation of current p-value CDF
  NumericVector f_eval(numValues);
  for(int i = 0; i < numCDF; i++) {
    checkUserInterrupt();
    
    int pos = 0;
    int len = sfuns[i].length();
    for(int j = 0; j < numValues; j++) {
      f_eval[j] = eval_pv(pvalues[j], sfuns[i], len, pos);
    }
    
    if(independence)
      pval_transf += (double)CDFcounts[i] * log(1 - f_eval);
    else 
      pval_transf += (double)CDFcounts[i] * f_eval;
  }
  
  if(independence)
    pval_transf = 1 - exp(pval_transf);
  
  // garbage collection
  delete[] sfuns;
  
  // compute adjustments
  for(int i = 0; i < numValues; i++)
    if(pval_transf[i] > 1.0) pval_transf[i] = 1.0;
  
  return pval_transf;
}

List kernel_DFWER_singlestep_crit(
  const List& pCDFlist,
  const NumericVector& support,
  const NumericVector& sorted_pv,
  const double alpha,
  const bool independence,
  const Nullable<IntegerVector>& pCDFcounts
) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // number of all attainable p-values in the support
  int numValues = support.length();
  
  // extract p-value CDF vectors
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  // get count of each unique p-value distribution
  IntegerVector CDFcounts;
  if(pCDFcounts.isNull() || as<IntegerVector>(pCDFcounts).length() == 0)
    CDFcounts = IntegerVector(numCDF, 1.0);
  else
    CDFcounts = pCDFcounts;
  
  // restrict support to values <= alpha (critical value cannot exceed alpha)
  int idx_max = binary_search(support, alpha, numValues);
  //NumericVector pv_list = support[Range(0, index_max)];
  
  // transform support with fast kernel
  NumericVector support_transf = kernel_DFWER_singlestep_fast(
    pCDFlist, support, independence, CDFcounts
  );
  
  // get index of critical value
  int idx_pval = binary_search(support_transf, alpha, idx_max + 1);
  // vector to store critical value
  NumericVector crit(1, support[idx_pval]);
  
  // store transformed sorted pvalues
  NumericVector pval_transf(numTests);
  // search for sorted p-values in 'pv_list' and save their adjustments
  idx_pval = 0;
  for(int i = 0; i < numTests; i++) {
    checkUserInterrupt();
    while(idx_pval < numValues && support[idx_pval] < sorted_pv[i]) idx_pval++;
    pval_transf[i] = std::min<double>(1.0, support_transf[idx_pval]);
  }
  
  // garbage collection
  delete[] sfuns;
  
  // return critical values and adjusted sorted p-values
  return List::create(Named("crit_consts") = crit, Named("pval_transf") = pval_transf);
}

NumericVector kernel_DFWER_stepwise_fast(
  const List& pCDFlist,
  const NumericVector& sorted_pv,
  const bool independence,
  const Nullable<List>& pCDFindices
) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // indices of the CDFs and their counts
  IntegerVector* CDFindices = new IntegerVector[numCDF];
  int* CDFcounts = new int[numCDF];
  if(pCDFindices.isNull() || as<List>(pCDFindices).length() == 0) {
    for(int i = 0; i < numCDF; i++) {
      CDFindices[i] = IntegerVector(1, i + 1);
      CDFcounts[i] = 1;
    }
  } else {
    for(int i = 0; i < numCDF; i++) {
      CDFindices[i] = as<IntegerVector>(as<List>(pCDFindices)[i]);
      CDFcounts[i] = CDFindices[i].length();
    }
  }
  // extract p-value CDF vectors and their lengths
  NumericVector* sfuns = new NumericVector[numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // evaluation of current p-value CDF
  NumericVector f_eval(numTests);
  for(int i = 0; i < numCDF; i++) {
    checkUserInterrupt();
    
    // current position in i-th CDF
    int pos = 0;
    // length of i-th CDF
    int len = sfuns[i].length();
    // current sorted p-value to which i-th CDF belongs
    int k = 0;
     
    for(int j = 0; j < CDFindices[i][CDFcounts[i] - 1]; j++) {
      // evaluate i-th CDF for all RELEVANT p-values and multiply with count 
      f_eval[j] = eval_pv(sorted_pv[j], sfuns[i], len, pos);
      //if(independence) f_eval[j] = std::log(1 - f_eval[j]);
      f_eval[j] *= (CDFcounts[i] - k);
      if(CDFindices[i][k] == j + 1) k++;
    }
    for(int j = CDFindices[i][CDFcounts[i] - 1]; j < numTests; j++)
      f_eval[j] = 0;
    
    // add evaluations to overall sums
    pval_transf += f_eval;
  }
  //if(independence)
    // revert log
  //  pval_transf = 1 - exp(pval_transf);
  
  // compute adjustments
  pval_transf[numTests - 1] = std::min<double>(1.0, pval_transf[numTests - 1]);
  if(independence)
    for(int i = numTests - 2; i >= 0; i--)
      pval_transf[i] = std::min<double>(pval_transf[i], pval_transf[i + 1]);
  else
    for(int i = 1; i < numTests; i++)
      pval_transf[i] = std::max<double>(pval_transf[i - 1], std::min<double>(1.0, pval_transf[i]));
  
  // garbage collection
  delete[] sfuns;
  delete[] CDFindices;
  delete[] CDFcounts;
  
  return pval_transf;
}

List kernel_DFWER_stepwise_crit(
    const List& pCDFlist,
    const NumericVector& support,
    const NumericVector& sorted_pv,
    const double alpha,
    const bool independence,
    const Nullable<List>& pCDFindices
) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // support size
  int numValues = support.length();
  
  // extract p-value CDF vectors and their lengths
  NumericVector* sfuns = new NumericVector[numCDF];
  int* lens = new int[numCDF];
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  // indices of the CDFs and their counts
  int* CDFcounts = new int[numCDF];
  int* pv2CDFindices = new int[numTests];
  if(pCDFindices.isNull() || as<List>(pCDFindices).length() == 0) {
    for(int i = 0; i < numCDF; i++) {
      CDFcounts[i] = 1;
      pv2CDFindices[i] = i;
    }
  } else {
    for(int i = 0; i < numCDF; i++) {
      IntegerVector CDFindices = as<IntegerVector>(as<List>(pCDFindices)[i]);
      CDFcounts[i] = CDFindices.length();
      for(int j = 0; j < CDFcounts[i]; j++)
        pv2CDFindices[CDFindices[j] - 1] = i;
    }
  }
  // threshold
  //double beta = independence ? -std::log(1 - alpha) : alpha;
  
  // vector to store transformed p-values
  NumericVector pval_transf;
  // evaluation of current p-value CDF
  NumericVector f_eval(numValues);
  
  // finding critical value of [d-Bonf]; reduce support first
  int limit = binary_search(support, alpha / numTests, numValues);
  NumericVector pv_list = support[Range(limit, numValues - 1)];
  numValues -= limit;
  limit = binary_search(pv_list, alpha, numValues);
  pv_list = pv_list[Range(0, limit)];
  numValues = limit + 1;
  pval_transf = NumericVector(limit + 1);
  for(int i = 0; i < numCDF; i++) {
    checkUserInterrupt();
    int pos = 0;
    for(int j = 0; j <= limit; j++) 
      f_eval[j] = eval_pv(pv_list[j], sfuns[i], lens[i], pos);
    //if(independence)
    //  pval_transf += -CDFcounts[i] * log(1 - f_eval);
    //else 
      pval_transf += CDFcounts[i] * f_eval;
  }
  
  int idx_pval = binary_search(pval_transf, alpha, limit + 1);
  pv_list = pv_list[Range(idx_pval, limit)];
  double crit_1 = pv_list[0];
  pv_list = sort_combine(pv_list, sorted_pv);
  numValues = pv_list.length();
  
  // critical values indices
  NumericVector crit(numTests, crit_1);
  // index of current critical value to be computed
  int idx_crit = numTests - 1;
  // vector to store transformed p-values
  pval_transf = NumericVector(numTests);
  // current position in transformed support for transforming observed p-values
  int idx_transf = numValues - 1;
  // vector to store CDF sums
  NumericVector pval_sums(numValues);
  // array to store if a p-value is in the current combined support
  LogicalVector supported(numValues);
  // number of observed p-values in i,...,m equal to the current one
  int count_pv = 0;
  // array for storing counts of unique CDFs of a p-value "block"
  int* CDFcounts_running = new int[numCDF];
  
  // search for critical values and transform observed p-values
  while(idx_crit >= 0) {
    checkUserInterrupt();
    
    // number of observed p-values equal to current one ("block" size)
    count_pv = 1;
    while(
      count_pv <= idx_crit && 
        sorted_pv[idx_crit - count_pv] == sorted_pv[idx_crit]
    )
      count_pv++;
    
    // find current p-value in support for adjustment
    while(idx_transf > 0 && pv_list[idx_transf] > sorted_pv[idx_crit]) 
      idx_transf--;
    
    // determine critical value and transformation
    if(count_pv == 1) {  // current p-value is unique
      // index of CDF belonging to current p-value
      int idx_CDF = pv2CDFindices[idx_crit];
      // vector for evaluating current CDF
      NumericVector f_eval(numValues);
      
      // evaluate CDF and add its attainable values to support
      for(int i = 0; i < numValues; i++) {
        int pos = 0;
        f_eval[i] = eval_pv(pv_list[i], sfuns[idx_CDF], lens[idx_CDF], pos);
        supported[i] = supported[i] || (f_eval[i] == pv_list[i]);
      }
      
      // add evaluations to overall sums
      //if(independence)
      //  pval_sums += -log(1 - f_eval);
      //else
        pval_sums += f_eval;
      
      // find critical value
      idx_pval = numValues - 1;
      while(
        idx_pval > 0 && 
          (pval_sums[idx_pval] > alpha || !supported[idx_pval])
      )
        idx_pval--;
      
      // save critical value
      crit[idx_crit] = pv_list[idx_pval];
      
      // compute transformed p-value
      if(pv_list[idx_transf] == sorted_pv[idx_crit]) 
      //  pval_transf[idx_crit] = independence 
      //    ? 1 - std::exp(-pval_sums[idx_transf])
      //    : pval_sums[idx_transf];
        pval_transf[idx_crit] = std::min<double>(1.0, pval_sums[idx_transf]);
      
      // go to next critical value
      idx_crit--;
    } else {  // current p-value is not unique (i.e. in a "block")
      // determine CDF counts for current "block"
      for(int i = 0; i < numCDF; i++) CDFcounts_running[i] = 0;
      // largest index of CDFs in "block"
      int max_CDF = 0;
      // index of current CDF
      int idx_CDF = 0;
      // determine counts and largest index
      for(int i = idx_crit - count_pv + 1; i <= idx_crit; i++) {
        idx_CDF = pv2CDFindices[i];
        max_CDF = std::max<int>(max_CDF, idx_CDF);
        CDFcounts_running[idx_CDF]++;
      }
      // index of last CDF for current p-value "block"
      int idx_last = idx_CDF;
      // last sum of current p-value
      double pval_sum_last = pval_sums[idx_transf];
      
      for(idx_CDF = 0; idx_CDF <= max_CDF; idx_CDF++) {
        if(CDFcounts_running[idx_CDF] > 0) {
          // vector for evaluating current CDF
          NumericVector f_eval(numValues);
          // evaluate CDF and add its attainable values to support
          for(int i = 0; i < numValues; i++) {
            int pos = 0;
            f_eval[i] = eval_pv(pv_list[i], sfuns[idx_CDF], lens[idx_CDF], pos);
            supported[i] = supported[i] || (f_eval[i] == pv_list[i]);
          }
          // add evaluations to overall sums
          //if(independence)
          //  pval_sums += -log(1 - f_eval) * CDFcounts_running[idx_CDF];
          //else
            pval_sums += f_eval * CDFcounts_running[idx_CDF];
          // compute adjustment for Hochberg procedure
          if(independence && idx_CDF == idx_last)
            pval_sum_last += f_eval[idx_transf];
        }
      }
      
      // find critical value
      idx_pval = numValues - 1;
      while(
        idx_pval > 0 && 
          (pval_sums[idx_pval] > alpha || !supported[idx_pval])
      )
        idx_pval--;
      
      // save critical values and transformed p-values
      for(int i = idx_crit - count_pv + 1; i <= idx_crit; i++) {
        // critical values
        crit[i] = pv_list[idx_pval];
        // transform p-value
        if(pv_list[idx_transf] == sorted_pv[idx_crit]) 
          //pval_transf[i] = independence
          //  ? 1 - std::exp(-pval_sums[idx_transf])
          //  : pval_sums[idx_transf];
          pval_transf[i] = independence
            ? std::min<double>(1.0, pval_sum_last)
            : std::min<double>(1.0, pval_sums[idx_transf]);
      }
      
      // go to next critical values
      idx_crit -= count_pv;
    }
  }
  
  /*if(independence) {
    pval_transf[numTests - 1] = std::min<double>(1.0, pval_transf[numTests - 1]);
    for(int i = numTests - 2; i >= 0; i--)
      pval_transf[i] = std::min<double>(pval_transf[i], pval_transf[i + 1]);
  } else {
    pval_transf[0] = std::min<double>(1.0, pval_transf[0]);
    for(int i = 1; i < numTests; i++)
      pval_transf[i] = std::max<double>(pval_transf[i - 1], std::min<double>(1.0, pval_transf[i]));
  }*/
  
  // garbage collection
  delete[] CDFcounts_running;
  delete[] pv2CDFindices;
  delete[] CDFcounts;
  delete[] lens;
  delete[] sfuns;
  
  // output results
  return List::create(Named("crit_consts") = crit, Named("pval_transf") = pval_transf);
}
/*
// [[Rcpp::export]]
List kernel_DFWER_stepwise_crit2(
    const List &pCDFlist,
    const NumericVector &support,
    const NumericVector &sorted_pv,
    const double alpha,
    const bool independence,
    const Nullable<List> &pCDFindices
) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // support size
  int numValues = support.length();
  
  // extract p-value CDF vectors and their lengths
  NumericVector* sfuns = new NumericVector[numCDF];
  int* lens = new int[numCDF];
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  // indices of the CDFs and their counts
  int* CDFcounts = new int[numCDF];
  int* pv2CDFindices = new int[numTests];
  if(pCDFindices.isNull() || as<List>(pCDFindices).length() == 0) {
    for(int i = 0; i < numCDF; i++) {
      CDFcounts[i] = 1;
      pv2CDFindices[i] = i;
    }
  } else {
    for(int i = 0; i < numCDF; i++) {
      IntegerVector CDFindices = as<IntegerVector>(as<List>(pCDFindices)[i]);
      CDFcounts[i] = CDFindices.length();
      for(int j = 0; j < CDFcounts[i]; j++)
        pv2CDFindices[CDFindices[j] - 1] = i;
    }
  }
  // threshold
  double beta = independence ? -std::log(1 - alpha) : alpha;
  
  // vector to store transformed p-values
  NumericVector pval_transf;
  // evaluation of current p-value CDF
  NumericVector f_eval(numValues);
  
  // finding critical value of [d-Bonf]; reduce support first
  int limit = binary_search(support, alpha / numTests, numValues);
  NumericVector pv_list = support[Range(limit, numValues - 1)];
  numValues -= limit;
  limit = binary_search(pv_list, alpha, numValues);
  pv_list = pv_list[Range(0, limit)];
  numValues = limit + 1;
  pval_transf = NumericVector(limit + 1);
  for(int i = 0; i < numCDF; i++) {
    checkUserInterrupt();
    int pos = 0;
    for(int j = 0; j <= limit; j++) {
      f_eval[j] = eval_pv(pv_list[j], sfuns[i], lens[i], pos);
    }
    if(independence)
      pval_transf += -CDFcounts[i] * log(1 - f_eval);
    else 
      pval_transf += CDFcounts[i] * f_eval;
  }
  
  int idx_pval = binary_search(pval_transf, beta, limit + 1);
  pv_list = pv_list[Range(idx_pval, limit)];
  double crit_1 = pv_list[0];
  pv_list = sort_combine(pv_list, sorted_pv);
  numValues = pv_list.length();
  
  // critical values indices
  NumericVector crit(numTests, crit_1);
  // index of current critical value to be computed
  int idx_crit = numTests - 1;
  // vector to store transformed p-values
  pval_transf = NumericVector(numTests);
  // current position in transformed support for transforming observed p-values
  int idx_transf = numValues - 1;
  // vector to store CDF sums
  NumericVector pval_sums(numValues);
  // array to store if a p-value is in the current combined support
  LogicalVector supported(numValues);
  // number of observed p-values in i,...,m equal to the current one
  int count_pv = 0;
  // array for storing counts of unique CDFs of a p-value "block"
  int* CDFcounts_running = new int[numCDF];
  
  // search for critical values and transform observed p-values
  while(idx_crit >= 0) {
    checkUserInterrupt();
    // number of observed p-values equal to current one ("block" size)
    count_pv = 1;
    while(
      count_pv <= idx_crit && 
        sorted_pv[idx_crit - count_pv] == sorted_pv[idx_crit]
    )
      count_pv++;
    // determine critical value and transformation
    if(count_pv == 1) {  // current p-value is unique
      // index of CDF belonging to current p-value
      int idx_CDF = pv2CDFindices[idx_crit];
      // vector for evaluating current CDF
      NumericVector f_eval(numValues);
      // evaluate CDF and add its attainable values to support
      for(int i = 0; i < numValues; i++) {
        int pos = 0;
        f_eval[i] = eval_pv(pv_list[i], sfuns[idx_CDF], lens[idx_CDF], pos);
        supported[i] = supported[i] || (f_eval[i] == pv_list[i]);
      }
      // add evaluations to overall sums
      if(independence)
        pval_sums += -log(1 - f_eval);
      else
        pval_sums += f_eval;
      // find critical value
      idx_pval = numValues - 1;
      while(
        idx_pval > 0 && 
          (pval_sums[idx_pval] > beta || !supported[idx_pval])
      )
        idx_pval--;
      // save critical value
      crit[idx_crit] = pv_list[idx_pval];
      Rcout << idx_crit << ": " << crit[idx_crit] << "\n";
      
      // compute adjusted p-value
      while(idx_transf > 0 && pv_list[idx_transf] > sorted_pv[idx_crit]) 
        idx_transf--;
      if(pv_list[idx_transf] == sorted_pv[idx_crit]) 
        pval_transf[idx_crit] = independence 
      ? 1 - std::exp(-pval_sums[idx_transf])
        : pval_sums[idx_transf];
      
      // go to next critical value
      idx_crit--;
    } else {  // current p-value is not unique (i.e. in a "block")
      // reset running CDF counts
      for(int i = 0; i < numCDF; i++) CDFcounts_running[i] = 0;
      // largest index of CDFs in "block"
      int max_CDF = 0;
      // determine counts and largest index
      for(int i = idx_crit - count_pv + 1; i <= idx_crit; i++) {
        int idx_CDF = pv2CDFindices[i];
        max_CDF = std::max<int>(max_CDF, idx_CDF);
        CDFcounts_running[idx_CDF]++;
      }
      
      // determine all critical and transformations in block
      while(count_pv > 0) {
        checkUserInterrupt();
        // temporary index of CDF that allows largest transformation <= alpha
        int idx_crit_max = 0;
        // its newly added p-values for running support
        LogicalVector supported_max(numValues);
        // its new evaluation vector
        NumericVector f_eval_max(numValues);
        // cycle through CDFs
        for(int idx_CDF = 0; idx_CDF <= max_CDF; idx_CDF++) {
          if(CDFcounts_running[idx_CDF] > 0) {
            // vector for evaluating current CDF
            NumericVector f_eval(numValues);
            // temporary array to store if p-values are in current CDF support
            LogicalVector supported_temp(numValues);
            // evaluate CDF and add its attainable values to support
            for(int i = 0; i < numValues; i++) {
              int pos = 0;
              f_eval[i] = eval_pv(pv_list[i], sfuns[idx_CDF], lens[idx_CDF], pos);
              supported_temp[i] = (f_eval[i] == pv_list[i]);
            }
            if(independence) f_eval = -log(1 - f_eval);
            // find critical value
            idx_pval = numValues - 1;
            while(
              idx_pval >= 0 && (
                pval_sums[idx_pval] + f_eval[idx_pval] > beta || 
                  (!supported[idx_pval] && !supported_temp[idx_pval])
              )
            )
              idx_pval--;
            
            // store best results
            if(idx_pval > 0 && crit[idx_crit] < pv_list[idx_pval]) {
              crit[idx_crit] = pv_list[idx_pval];
              Rcout << idx_crit << ": " << crit[idx_crit] << " (" << pval_sums[idx_pval] + f_eval[idx_pval] << ", " << idx_CDF << ", " << CDFcounts_running[idx_CDF] << ")\n";
              idx_crit_max = idx_CDF;
              f_eval_max = clone(f_eval);
              supported_max = clone(supported_temp);
            }
          }
        }
        
        // update sums
        pval_sums += f_eval_max;
        // update support
        for(idx_pval = 0; idx_pval < numValues; idx_pval++)
          supported[idx_pval] = supported[idx_pval] || supported_max[idx_pval];
        
        // compute adjusted p-value
        while(idx_transf > 0 && pv_list[idx_transf] > sorted_pv[idx_crit]) 
          idx_transf--;
        if(pv_list[idx_transf] == sorted_pv[idx_crit]) 
          pval_transf[idx_crit] = independence
        ? 1 - std::exp(-pval_sums[idx_transf])
          : pval_sums[idx_transf];
        
        // reduce count of selected CDF
        CDFcounts_running[idx_crit_max]--;
        // reduce p-value count
        count_pv--;
        // go to next critical value
        idx_crit--;
      }
    }
  }
  
  // garbage collection
  delete[] CDFcounts_running;
  delete[] pv2CDFindices;
  delete[] CDFcounts;
  delete[] lens;
  delete[] sfuns;
  
  // output results
  return List::create(Named("crit_consts") = crit, Named("pval_transf") = pval_transf);
}
*/
