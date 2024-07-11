#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <armadillo>
#include "hlmgwr.h"
#include "helper.h"

using namespace std;
using namespace arma;
using namespace hgwr;

void pcout(const string& message)
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

    // SECTION("Specified bandwidth 10") {
    //     auto kernel = HGWR::KernelType::GAUSSIAN;
    //     double bw = 10.0;
    //     HGWR::Options options { 0.1, 1e-6, 1e-6, 100000, 10, 0, 0 };
    //     HGWR algorithm { G, X, Z, y, u, group, kernel, bw, options, pcout };
    //     REQUIRE_NOTHROW(algorithm.fit());
    //     INFO("Results:");
    //     CAPTURE(algorithm.get_bw(), algorithm.get_sigma(), algorithm.get_gamma(), algorithm.get_beta(), algorithm.get_mu(), algorithm.get_D());
    //     CHECK_THAT(algorithm.get_bw(), Catch::Matchers::WithinAbs(10.0, 1e-6));
    //     CHECK_THAT(algorithm.get_sigma(), Catch::Matchers::WithinAbs(1.95, 1e-2));
    // }

    SECTION("Optimise bandwidth") {
        auto kernel = HGWR::KernelType::BISQUARED;
        HGWR::Options options { 1e-16, 1e-6, 1e-6, 100000, 10, 0, 0 };
        HGWR algorithm { G, X, Z, y, u, group, kernel, options };
        algorithm.set_printer(pcout);
        REQUIRE_NOTHROW(algorithm.fit());
        INFO("Results:");
        CAPTURE(algorithm.get_bw(), algorithm.get_sigma(), algorithm.get_gamma(), algorithm.get_beta(), algorithm.get_mu(), algorithm.get_D());
        CHECK_THAT(algorithm.get_bw(), Catch::Matchers::WithinAbs(8, 1));
        CHECK_THAT(algorithm.get_sigma(), Catch::Matchers::WithinAbs(1.90, 1e-2));
    }
}
