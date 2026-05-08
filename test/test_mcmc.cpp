#include <iostream>
#include <armadillo>
#include "hlmgwr.h"

using namespace std;
using namespace arma;
using namespace hgwr;

void printer(const string& msg)
{
    cout << msg;
}

int main()
{
    cout << "========================================" << endl;
    cout << "  HGWR MCMC Test" << endl;
    cout << "========================================" << endl;

    mat G, X, Z, u;
    vec y;
    uvec group;

    string data_dir = string(TEST_DATA_DIR);

    X.load(csv_name(data_dir + "/hlmgwr_x.csv"));
    G.load(csv_name(data_dir + "/hlmgwr_g.csv"));
    Z.load(csv_name(data_dir + "/hlmgwr_z.csv"));
    u.load(csv_name(data_dir + "/hlmgwr_u.csv"));
    y.load(csv_name(data_dir + "/hlmgwr_y.csv"));
    group.load(csv_name(data_dir + "/hlmgwr_group.csv"));

    cout << "Data loaded: " << X.n_rows << " observations, "
         << G.n_rows << " groups" << endl;
    cout << "G cols: " << G.n_cols << ", X cols: " << X.n_cols
         << ", Z cols: " << Z.n_cols << endl;
    cout << endl;

    auto kernel = HGWR::KernelType::GAUSSIAN;
    double bw = 10.0;

    {
        cout << "========================================" << endl;
        cout << "  Test 1: MLE (ml_type=0, original method)" << endl;
        cout << "========================================" << endl;
        HGWR::Options options_mle;
        options_mle.alpha = 0.1;
        options_mle.eps_iter = 1e-6;
        options_mle.eps_gradient = 1e-6;
        options_mle.max_iters = 100000;
        options_mle.max_retries = 10;
        options_mle.verbose = 0;
        options_mle.ml_type = 0;

        HGWR alg_mle(G, X, Z, y, u, group, kernel, bw, options_mle);
        auto res_mle = alg_mle.fit();

        cout << "  MLE Results:" << endl;
        cout << "  bw    = " << res_mle.bw << endl;
        cout << "  sigma = " << res_mle.sigma << endl;
        cout << "  beta  = " << res_mle.beta.t() << endl;
        cout << "  D     = " << endl << res_mle.D << endl;
        cout << endl;
    }

    {
        cout << "========================================" << endl;
        cout << "  Test 2: MCMC (ml_type=2, Bayesian)" << endl;
        cout << "========================================" << endl;
        HGWR::Options options_mcmc;
        options_mcmc.alpha = 0.1;
        options_mcmc.eps_iter = 1e-6;
        options_mcmc.eps_gradient = 1e-6;
        options_mcmc.max_iters = 10000;
        options_mcmc.max_retries = 10;
        options_mcmc.verbose = 1;
        options_mcmc.ml_type = 2;

        HGWR alg_mcmc(G, X, Z, y, u, group, kernel, bw, options_mcmc);
        alg_mcmc.set_printer(printer);
        auto res_mcmc = alg_mcmc.fit_mcmc_backfitting();

        cout << "  MCMC Results:" << endl;
        cout << "  bw    = " << res_mcmc.bw << endl;
        cout << "  sigma = " << res_mcmc.sigma << endl;
        cout << "  beta  = " << res_mcmc.beta.t() << endl;
        cout << "  D     = " << endl << res_mcmc.D << endl;
        cout << endl;
    }

    {
        cout << "========================================" << endl;
        cout << "  Test 3: MCMC with fewer iterations" << endl;
        cout << "========================================" << endl;
        HGWR::Options options_mcmc;
        options_mcmc.alpha = 0.1;
        options_mcmc.eps_iter = 1e-6;
        options_mcmc.eps_gradient = 1e-6;
        options_mcmc.max_iters = 1000;
        options_mcmc.max_retries = 10;
        options_mcmc.verbose = 0;
        options_mcmc.ml_type = 2;

        HGWR alg_mcmc(G, X, Z, y, u, group, kernel, bw, options_mcmc);
        alg_mcmc.set_printer(printer);
        auto res_mcmc = alg_mcmc.fit();

        cout << "  MCMC (2000 iters) Results:" << endl;
        cout << "  bw    = " << res_mcmc.bw << endl;
        cout << "  sigma = " << res_mcmc.sigma << endl;
        cout << "  beta  = " << res_mcmc.beta.t() << endl;
        cout << "  D     = " << endl << res_mcmc.D << endl;
    }

    cout << "========================================" << endl;
    cout << "  ALL TESTS PASSED" << endl;
    cout << "========================================" << endl;

    return 0;
}
