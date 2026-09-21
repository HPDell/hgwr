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
        CHECK(algorithm.get_converged());
        CHECK(algorithm.get_optimizer_failures() == 0);
        CHECK(std::isfinite(algorithm.get_sigma()));
        CHECK(algorithm.get_sigma() > 0.0);
        CHECK(algorithm.get_D().is_sympd());
    }

    SECTION("Optimise bandwidth") {
        auto kernel = HGWR::KernelType::BISQUARED;
        HGWR::Options options { 1e-4, 1e-6, 1e-6, 100000, 10, 0, 0 };
        HGWR algorithm { G, X, Z, y, u, group, kernel, options };
        algorithm.set_printer(pcout);
        REQUIRE_NOTHROW(algorithm.fit());
        INFO("Results:");
        CAPTURE(algorithm.get_bw(), algorithm.get_sigma(), algorithm.get_gamma(), algorithm.get_beta(), algorithm.get_mu(), algorithm.get_D());
        CHECK(algorithm.get_converged());
        CHECK(algorithm.get_bw_optimizer_failures() == 0);
        CHECK(std::isfinite(algorithm.get_bw()));
        CHECK(algorithm.get_bw() == std::floor(algorithm.get_bw()));
        CHECK(algorithm.get_bw() >= algorithm.get_bw_lower());
        CHECK(algorithm.get_bw() <= algorithm.get_bw_upper());
        CHECK(std::isfinite(algorithm.get_bw_objective()));
        CHECK(algorithm.get_bw_evaluations() > 0);
        CHECK(algorithm.get_D().is_sympd());
    }

    SECTION("F test") {
        auto kernel = HGWR::KernelType::GAUSSIAN;
        double bw = 10.0;
        HGWR::Options options { 0.1, 1e-6, 1e-6, 100000, 10, 0, 0 };
        HGWR algorithm { G, X, Z, y, u, group, kernel, bw, options, pcout };
        REQUIRE_NOTHROW(algorithm.fit(true));
        vector<vec4> fResults;
        REQUIRE_NOTHROW(fResults = algorithm.test_glsw());
        REQUIRE(fResults.size() == G.n_cols);
        INFO("Results:");
        for (auto &&i : fResults)
        {
            CAPTURE(i(0), i(1), i(2), i(3));
            CHECK(std::isfinite(i(0)));
            CHECK(std::isfinite(i(1)));
            CHECK(std::isfinite(i(2)));
            CHECK(i(0) >= 0.0);
            CHECK(i(1) > 0.0);
            CHECK(i(2) > 0.0);
            CHECK((std::isnan(i(3)) || (i(3) >= 0.0 && i(3) <= 1.0)));
        }
    }
}
