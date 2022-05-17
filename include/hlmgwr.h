#ifndef HLMGWR_H
#define HLMGWR_H

#include <string>

#ifndef HGWRR_RCPP
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

enum GWRKernelType {
    GAUSSIAN,
    BISQUARED
};

struct HLMGWRArgs {
    arma::mat G;
    arma::mat X;
    arma::mat Z;
    arma::vec y;
    arma::mat u;
    arma::uvec group;
    double bw;
    GWRKernelType kernel;

    HLMGWRArgs() : G(), X(), Z(), y(), u(), group(), bw(0.0), kernel(GWRKernelType::GAUSSIAN)
    {
    }

    HLMGWRArgs(arma::mat in_G, arma::mat in_X, arma::mat in_Z, arma::vec in_y, arma::mat in_u, arma::uvec in_group, double in_bw, GWRKernelType in_kernel) :
        G(in_G),
        X(in_X),
        Z(in_Z),
        y(in_y),
        u(in_u),
        group(in_group),
        bw(in_bw),
        kernel(in_kernel)
    {
    }
};

struct HLMGWRParams {
    arma::mat gamma;
    arma::mat beta;
    arma::mat mu;
    arma::mat D;
    double sigma;
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

    HLMGWROptions(
        double in_alpha,
        double in_eps_iter,
        double in_eps_gradient,
        size_t in_max_iters,
        size_t in_max_retries,
        size_t in_verbose,
        size_t in_ml_type
    )
    {
        alpha = in_alpha;
        eps_iter = in_eps_iter;
        eps_gradient = in_eps_gradient;
        max_iters = in_max_iters;
        max_retries = in_max_retries;
        verbose = in_verbose;
        ml_type = in_ml_type;
    }
};

typedef void (*PrintFunction)(std::string);

HLMGWRParams backfitting_maximum_likelihood(const HLMGWRArgs& args, const HLMGWROptions& options, const PrintFunction pcout);

#endif  // HLMGWR_H
