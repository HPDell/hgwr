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
    double alpha;
    double eps_iter;
    double eps_gradient;
    size_t max_iters;
    size_t max_retries;
    size_t verbose;
    size_t ml_type;
};

HLMGWRParams backfitting_maximum_likelihood(const HLMGWRArgs& args, const HLMGWROptions& options);

#endif  // HLMGWR_H
