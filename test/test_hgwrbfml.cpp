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

    SECTION("Specified bandwidth 10") {
        auto kernel = HGWR::KernelType::GAUSSIAN;
        double bw = 10.0;
        HGWR::Options options { 0.1, 1e-6, 1e-6, 100000, 10, 0, 0 };
        HGWR algorithm { G, X, Z, y, u, group, kernel, bw, options, pcout };
        REQUIRE_NOTHROW(algorithm.fit());
        INFO("Results:");
        CAPTURE(algorithm.get_bw(), algorithm.get_sigma(), algorithm.get_gamma(), algorithm.get_beta(), algorithm.get_mu(), algorithm.get_D());
        CHECK_THAT(algorithm.get_bw(), Catch::Matchers::WithinAbs(10.0, 1e-6));
        CHECK_THAT(algorithm.get_sigma(), Catch::Matchers::WithinAbs(1.966, 1e-2));
        CHECK(eig_sym(algorithm.get_D()).min() > 0.0);
        CHECK_FALSE(approx_equal(algorithm.get_D(), eye<mat>(algorithm.get_D().n_rows, algorithm.get_D().n_cols), "absdiff", 1e-6));
    }

    SECTION("Optimise bandwidth") {
        auto kernel = HGWR::KernelType::BISQUARED;
        HGWR::Options options { 0.01, 1e-6, 1e-6, 100000, 10, 0, 0 };
        HGWR algorithm { G, X, Z, y, u, group, kernel, options };
        algorithm.set_printer(pcout);
        REQUIRE_NOTHROW(algorithm.fit());
        INFO("Results:");
        CAPTURE(algorithm.get_bw(), algorithm.get_sigma(), algorithm.get_gamma(), algorithm.get_beta(), algorithm.get_mu(), algorithm.get_D());
        CHECK_THAT(algorithm.get_bw(), Catch::Matchers::WithinAbs(8, 1));
        CHECK_THAT(algorithm.get_sigma(), Catch::Matchers::WithinAbs(1.966, 1e-2));
        CHECK(eig_sym(algorithm.get_D()).min() > 0.0);
    }

    SECTION("F test") {
        auto kernel = HGWR::KernelType::GAUSSIAN;
        double bw = 10.0;
        HGWR::Options options { 0.1, 1e-6, 1e-6, 100000, 10, 0, 0 };
        HGWR algorithm { G, X, Z, y, u, group, kernel, bw, options, pcout };
        REQUIRE_NOTHROW(algorithm.fit());
        vector<vec4> fResults = algorithm.test_glsw();
        INFO("Results:");
        for (auto &&i : fResults)
        {
            CAPTURE(i(0), i(1), i(2), i(3));
        }
    }
}

TEST_CASE("HGWR covariance gradients match finite differences")
{
    mat X,Z;
    vec y;
    uvec group;
    X.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_x.csv"));
    Z.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_z.csv"));
    y.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_y.csv"));
    group.load(arma::csv_name(string(TEST_DATA_DIR) + "/hlmgwr_group.csv"));
    vec beta = solve(X.t() * X, X.t() * y);
    mat D = {{0.5, 0.1}, {0.1, 0.3}};
    CHECK(ml_gradient_relative_error(X, Z, y, group, beta, D, false) < 1e-5);
    CHECK(ml_gradient_relative_error(X, Z, y, group, beta, D, true) < 1e-5);
}
