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
        REQUIRE_NOTHROW(algorithm.fit(true));
        vector<vec4> fResults = algorithm.test_glsw();
        REQUIRE(fResults.size() == G.n_cols);

        const uword n = y.n_elem;
        const uword m = G.n_rows;
        const uword p = G.n_cols;
        vector<uvec> indices(m);
        vector<mat> V(m), V_inv(m), GV(m), GVG(m);
        mat V_full(n, n, fill::zeros);
        vector<mat> coefficient_map(p);
        for (uword k = 0; k < p; ++k) coefficient_map[k].zeros(m, n);
        for (uword j = 0; j < m; ++j)
        {
            indices[j] = find(group == j);
            mat Zj = Z.rows(indices[j]);
            V[j] = Zj * algorithm.get_D() * Zj.t() + eye(indices[j].n_elem, indices[j].n_elem);
            V_inv[j] = inv_sympd(V[j]);
            V_full.submat(indices[j], indices[j]) = V[j];
            GV[j] = G.row(j).t() * ones<rowvec>(indices[j].n_elem) * V_inv[j];
            GVG[j] = GV[j] * ones<vec>(indices[j].n_elem) * G.row(j);
        }

        mat smoother_by_group(m, n, fill::zeros);
        for (uword i = 0; i < m; ++i)
        {
            mat delta_u = u.each_row() - u.row(i);
            vec distance = sqrt(sum(delta_u % delta_u, 1));
            double local_bw = HGWR::actual_bw(distance, bw);
            vec weight = HGWR::gwr_kernel_gaussian2(distance % distance, local_bw * local_bw);
            mat lhs(p, p, fill::zeros), rhs(p, n, fill::zeros);
            for (uword j = 0; j < m; ++j)
            {
                lhs += weight(j) * GVG[j];
                rhs.cols(indices[j]) = weight(j) * GV[j];
            }
            mat Cit = rhs.t() * inv(lhs).t();
            for (uword k = 0; k < p; ++k) coefficient_map[k].row(i) = Cit.col(k).t();
            smoother_by_group.row(i) = (Cit * G.row(i).t()).t();
        }

        mat S = smoother_by_group.rows(group);
        mat residual_maker = eye(n, n) - S;
        mat QV = residual_maker.t() * residual_maker * V_full;
        double delta1 = trace(QV);
        double delta2 = trace(QV * QV);
        double expected_df2 = delta1 * delta1 / delta2;

        vec fitted_glsw = sum(G.rows(group) % algorithm.get_gamma().rows(group), 1);
        vec residual = y - fitted_glsw - X * algorithm.get_beta();
        double expected_sigma2 = 0.0;
        for (uword j = 0; j < m; ++j)
        {
            vec rj = residual.rows(indices[j]);
            expected_sigma2 += as_scalar(rj.t() * V_inv[j] * rj);
        }
        expected_sigma2 /= double(n);
        CHECK_THAT(algorithm.get_sigma() * algorithm.get_sigma(),
            Catch::Matchers::WithinRel(expected_sigma2, 1e-10));

        INFO("Results:");
        for (uword k = 0; k < p; ++k)
        {
            vec gamma_k = algorithm.get_gamma().col(k);
            vec expanded_gamma = gamma_k.rows(group);
            double vk2 = accu(square(expanded_gamma - mean(expanded_gamma))) / double(n);
            mat expanded_map = coefficient_map[k].rows(group);
            mat centring = eye(n, n) - ones<mat>(n, n) / double(n);
            mat BV = expanded_map.t() * centring * expanded_map * V_full / double(n);
            double xi1 = trace(BV);
            double xi2 = trace(BV * BV);
            double expected_f = vk2 / (xi1 * expected_sigma2);
            double expected_df1 = xi1 * xi1 / xi2;

            const vec4& result = fResults[k];
            CAPTURE(result(0), result(1), result(2), result(3));
            CHECK(result.is_finite());
            CHECK_THAT(result(0), Catch::Matchers::WithinRel(expected_f, 1e-10));
            CHECK_THAT(result(1), Catch::Matchers::WithinRel(expected_df1, 1e-10));
            CHECK_THAT(result(2), Catch::Matchers::WithinRel(expected_df2, 1e-10));
            CHECK(result(3) >= 0.0);
            CHECK(result(3) <= 1.0);
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
