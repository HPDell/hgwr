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

TEST_CASE("HGWR Multiscale (BFML)")
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

    SECTION("Optimise bandwidth per column") {
        auto kernel = HGWR::KernelType::BISQUARED;
        HGWR::Options options { 1e-16, 1e-6, 1e-6, 100000, 10, 0, 0, true };
        HGWR algorithm { G, X, Z, y, u, group, kernel, options };
        algorithm.set_printer(pcout);
        REQUIRE_NOTHROW(algorithm.fit());
        INFO("Results:");
        CAPTURE(algorithm.get_bw(), algorithm.get_sigma(), algorithm.get_gamma(), algorithm.get_beta(), algorithm.get_mu(), algorithm.get_D());
        const arma::vec& bws = algorithm.get_bws();
        REQUIRE(bws.n_elem == G.n_cols);
        for (arma::uword k = 0; k < bws.n_elem; k++)
        {
            CHECK_THAT(bws(k), Catch::Matchers::WithinAbs(8, 6));
        }
        CHECK_THAT(algorithm.get_sigma(), Catch::Matchers::WithinAbs(1.90, 5e-2));
    }

    SECTION("Specified per-column bandwidths") {
        auto kernel = HGWR::KernelType::GAUSSIAN;
        // Assign different bandwidths per column
        arma::vec bws(G.n_cols);
        bws(0) = 10.0;
        bws(1) = 12.0;
        bws(2) = 8.0;
        HGWR::Options options { 0.1, 1e-6, 1e-6, 100000, 10, 0, 0, true };
        HGWR algorithm { G, X, Z, y, u, group, kernel, bws, options };
        algorithm.set_printer(pcout);
        REQUIRE(algorithm.get_multiscale());
        REQUIRE_NOTHROW(algorithm.fit());
        INFO("Results:");
        CAPTURE(algorithm.get_bw(), algorithm.get_sigma(), algorithm.get_gamma());
        const arma::vec& result_bws = algorithm.get_bws();
        CHECK_THAT(result_bws(0), Catch::Matchers::WithinAbs(10.0, 1e-6));
        CHECK_THAT(result_bws(1), Catch::Matchers::WithinAbs(12.0, 1e-6));
        CHECK_THAT(result_bws(2), Catch::Matchers::WithinAbs(8.0, 1e-6));
        CHECK_THAT(algorithm.get_sigma(), Catch::Matchers::WithinAbs(1.95, 3e-2));
    }

    SECTION("F test") {
        auto kernel = HGWR::KernelType::GAUSSIAN;
        arma::vec bws(G.n_cols);
        bws.fill(10.0);
        HGWR::Options options { 0.1, 1e-6, 1e-6, 100000, 10, 0, 0, true };
        HGWR algorithm { G, X, Z, y, u, group, kernel, bws, options, pcout };
        REQUIRE_NOTHROW(algorithm.fit());
        vector<vec4> fResults = algorithm.test_glsw();
        INFO("Results:");
        REQUIRE(fResults.size() == G.n_cols);
        for (auto &&i : fResults)
        {
            CAPTURE(i(0), i(1), i(2), i(3));
        }
    }

    SECTION("Switch mode at runtime") {
        auto kernel = HGWR::KernelType::GAUSSIAN;
        double bw = 10.0;
        HGWR::Options options { 0.1, 1e-6, 1e-6, 100000, 10, 0, 0 };
        HGWR algorithm { G, X, Z, y, u, group, kernel, bw, options, pcout };
        REQUIRE_FALSE(algorithm.get_multiscale());
        algorithm.set_multiscale(true);
        REQUIRE(algorithm.get_multiscale());
        REQUIRE_NOTHROW(algorithm.fit());
        CHECK_THAT(algorithm.get_sigma(), Catch::Matchers::WithinAbs(1.95, 3e-2));
        // All bandwidths should equal 10.0 since we specified bw=10
        const arma::vec& result_bws = algorithm.get_bws();
        for (arma::uword k = 0; k < result_bws.n_elem; k++)
        {
            CHECK_THAT(result_bws(k), Catch::Matchers::WithinAbs(10.0, 1e-6));
        }
    }
}
