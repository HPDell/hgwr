#ifndef HLMGWR_H
#define HLMGWR_H

#include <armadillo>

struct HLMGWRArgs {
    arma::mat G;
    arma::mat X;
    arma::mat Z;
    arma::vec y;
    arma::mat u;
    arma::uvec group;
    double bw;
};

struct HLMGWRParams {
    arma::mat gamma;
    arma::mat beta;
    arma::mat mu;
    arma::mat D;
};

struct HLMGWROptions {
    double alpha = 0.01;
    double eps_iter = 1e-6;
    double eps_gradient = 1e-6;
    size_t max_iters = (size_t)1e6;
    size_t max_retries = 10;
    size_t verbose = 0;
    size_t ml_type = 0;
};

HLMGWRParams backfitting_maximum_likelihood(const HLMGWRArgs& args, const HLMGWROptions& options);

#endif  // HLMGWR_H
