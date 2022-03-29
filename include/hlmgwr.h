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

    HLMGWRArgs() : G(), X(), Z(), y(), u(), group(), bw(0.0)
    {
    }

    HLMGWRArgs(arma::mat in_G, arma::mat in_X, arma::mat in_Z, arma::vec in_y, arma::mat in_u, arma::uvec in_group, double in_bw) :
        G(in_G),
        X(in_X),
        Z(in_Z),
        y(in_y),
        u(in_u),
        group(in_group),
        bw(in_bw)
    {
    }
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

    HLMGWROptions() 
    {
        alpha = 0.01;
        eps_iter = 1e-6;
        eps_gradient = 1e-6;
        max_iters = (size_t)1e6;
        max_retries = (size_t)10;
        verbose = (size_t)0;
        ml_type = (size_t)0;
    }
};

HLMGWRParams backfitting_maximum_likelihood(const HLMGWRArgs& args, const HLMGWROptions& options);

#endif  // HLMGWR_H
