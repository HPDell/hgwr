#ifndef HGWRR_RCPP
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

void save_mongo(const size_t iter, const arma::mat& gamma, const arma::vec& beta, const arma::mat& mu, const double rss, const double diff, const double loglik);