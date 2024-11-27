#' @title
#' FWER-Based Multiple Testing Procedures with Adaptation for Discrete Tests
#'
#' @docType package
#' @importFrom Rcpp evalCpp
#' @useDynLib DiscreteFWER
#' 
#' @description
#' This package implements adaptions for discrete tests of the Bonferroni, Holm,
#' Hochberg and Šidák procedures for control of the family-wise error rate
#' (FWER). 
#' 
#' @details
#' The main function [`discrete_FWER()`] makes all four procedures available to
#' the user. [`DBonferroni()`], [`DHolm()`], [`DHochberg()`] and [`DSidak()`]
#' are wrapper functions that enable the user to access them directly. Their
#' main parameters are either a 
#' [`DiscreteTestResults`][DiscreteTests::DiscreteTestResults] object from
#' package [DiscreteTests][DiscreteTests::DiscreteTests-package] or a vector of
#' raw observed p-values and a list whose elements are the discrete supports
#' of the CDFs of the \eqn{p}-values.
#' 
#' The function [`direct_discrete_FWER()`]is a wrapper for
#' [`DiscreteFDR::generate.pvalues()`] and [`discrete_FWER()`], which applies
#' discrete procedures directly to data.
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
"_PACKAGE"
