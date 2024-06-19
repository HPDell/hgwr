#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <armadillo>
#include "hlmgwr.h"

using namespace std;
using namespace arma;

void pcout(string message)
{
    cout << message;
}

TEST_CASE("HGWR(BFML)")
{
    mat G,X,Z,u;
    vec y;
    uvec group;
    X.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_x.csv"));
    G.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_g.csv"));
    Z.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_z.csv"));
    u.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_u.csv"));
    y.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_y.csv"));
    group.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_group.csv"));

    double bw = 10.0;
    GWRKernelType kernel = GWRKernelType::GAUSSIAN;
    HLMGWROptions options { 0.01, 1e-6, 1e-6, 100000, 10, 0, 0 };

    HLMGWRArgs alg_args { G, X, Z, y, u, group, bw, kernel };
    HLMGWRParams alg_params;
    REQUIRE_NOTHROW([&](){
        alg_params = backfitting_maximum_likelihood(alg_args, options, pcout);
        return 0;
    }());
    INFO("Results:");
    CAPTURE(alg_params.bw, alg_params.sigma, alg_params.gamma, alg_params.beta, alg_params.mu, alg_params.D);
    CHECK_THAT(alg_params.bw, Catch::Matchers::WithinAbs(10.0, 1e-6));
    CHECK_THAT(alg_params.sigma, Catch::Matchers::WithinAbs(1.9212546213601509, 1e-6));
}
