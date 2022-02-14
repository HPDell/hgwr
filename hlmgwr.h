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

#endif  // HLMGWR_H