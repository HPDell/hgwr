#include "hlmgwr.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <armadillo>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include <boost/program_options.hpp>
#include "helper.h"

using namespace std;
using namespace arma;

const double log2pi = log(2.0 * M_PI);

struct ML_D_Params
{
    const field<mat>* Xf;
    const field<vec>* Yf;
    const field<mat>* Zf;
    const vec* beta;
    uword n;
    uword q;
};


/**
 * @brief Estimate $\gamma$.
 * 
 * @param X Equals to $g$
 * @param y Equals to $\bar{y}$
 * @param S Equals to $s$
 * @param u Used to calculate $W$
 * @param bw Bandwidth 
 * @param wn Equals to $N$
 * @param wD Equals to $D$
 * @return mat 
 */
mat fit_gwr(const mat& G, const field<vec>& Yf, const field<mat>& Zf, const mat& D, const mat& u, double bw)
{
    uword ng = G.n_rows, k = G.n_cols, q = Zf(0).n_cols;
    mat beta(ng, k, arma::fill::zeros), D_inv = D.i();
    mat Vig(ng, k, arma::fill::zeros);
    vec Viy(ng, arma::fill::zeros);
    for (int i = 0; i < ng; i++)
    {
        const mat& Yi = Yf(i);
        const mat& Zi = Zf(i);
        uword ndata = Zi.n_rows;
        mat Vi_inv = eye(ndata, ndata) - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
        Vig.row(i) = accu(Vi_inv) * G.row(i);
        Viy(i) = as_scalar(sum(Vi_inv, 0) * Yi);
    }
    /// Calibrate for each gorup.
    for (int i = 0; i < ng; i++)
    {
        mat d_u = u.each_row() - u.row(i);
        vec d2 = sum(d_u % d_u, 1);
        double b2 = vec(sort(d2))[(int)bw];
        vec wW = exp(- d2 / (2.0 * b2));
        mat GtWVG = G.t() * (Vig.each_col() % wW);
        mat GtWVy = G.t() * (Viy % wW);
        beta.row(i) = solve(GtWVG, GtWVy).t();
    }
    return beta;
}

vec fit_gls(const field<mat>& Xf, const field<vec>& Yf, const field<mat>& Zf, const mat& D)
{
    uword ngroup = Xf.n_rows, p = Xf(0).n_cols;
    mat XtWX(p, p, arma::fill::zeros);
    vec XtWY(p, arma::fill::zeros);
    mat D_inv = D.i();
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf(i);
        const mat& Yi = Yf(i);
        const mat& Zi = Zf(i);
        uword ndata = Zi.n_rows;
        mat Vi_inv = eye(ndata, ndata) - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
        XtWX += Xi.t() * Vi_inv * Xi;
        XtWY += Xi.t() * Vi_inv * Yi;
    }
    return solve(XtWX, XtWY);
}

double loglikelihood_ml(const field<mat>& Xf, const field<vec>& Yf, const field<mat>& Zf, const mat& D, const vec& beta, const uword& ndata)
{
    uword ngroup = Xf.n_rows;
    mat D_inv = D.i();
    double L1 = 0.0, L2 = 0.0;
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf(i);
        const vec& Yi = Yf(i);
        const mat& Zi = Zf(i);
        uword ndata = Zi.n_rows;
        mat Ii = eye<mat>(ndata, ndata);
        mat Vi = ((Zi * D) * Zi.t()) + Ii;
        double detVi, sign_detVi;
        log_det(detVi, sign_detVi, Vi);
        mat Vi_inv = Ii - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
        vec Ri = Yi - Xi * beta;
        L1 += as_scalar(Ri.t() * Vi_inv * Ri);
        L2 += detVi;
    }
    double LL = - (ndata / 2.0) * log(L1) - 0.5 * L2 - 0.5 - 0.5 * log2pi + (ndata / 2.0) * log(ndata);
    return LL;
}

mat loglikelihood_ml_d(const field<mat>& Xf, const field<vec>& Yf, const field<mat>& Zf, const mat& D, const vec& beta, const uword& ndata)
{
    mat J(1, 1, arma::fill::zeros), ZtViZ(arma::size(D), arma::fill::zeros), D_inv = D.i();
    uword ngroup = Xf.n_rows;
    field<mat> Kf(ngroup);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf(i);
        const mat& Yi = Yf(i);
        const mat& Zi = Zf(i);
        uword ndata = Zi.n_rows;
        mat Vi = Zi * D * Zi.t() + eye(ndata, ndata);
        mat Vi_inv = eye(ndata, ndata) - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
        vec Ri = Yi - Xi * beta;
        Kf(i) = Zi.t() * Vi_inv * Ri;
        mat Ji = Ri.t() * Vi_inv * Ri;
        ZtViZ += Zi.t() * Vi_inv * Zi;
    }
    mat KJKt(arma::size(D), arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Ki = Kf(i);
        mat KJKt_i = Ki * J * Ki.t();
        KJKt += KJKt_i;
    }
    mat dL_D = (ndata / 2.0) * KJKt - 0.5 * ZtViZ;
    return dL_D;
}

double loglikelihood_ml_gsl(const gsl_vector* v, void* p)
{
    ML_D_Params* params = (ML_D_Params*)p;
    const field<mat>* Xf = params->Xf;
    const field<vec>* Yf = params->Yf;
    const field<mat>* Zf = params->Zf;
    const vec* beta = params->beta;
    const uword n = params->n;
    const uword q = params->q;
    size_t ntarget = q * (q + 1) / 2;
    vec D_tri(ntarget, arma::fill::zeros);
    for (size_t i = 0; i < ntarget; i++)
    {
        D_tri(i) = gsl_vector_get(v, i);
    }
    mat D(q, q, arma::fill::zeros);
    D(trimatl_ind(size(D))) = D_tri;
    D(trimatu_ind(size(D))) = D_tri;
    return loglikelihood_ml(*Xf, *Yf, *Zf, D, *beta, n);
}

void loglikelihood_ml_d_gsl(const gsl_vector* v, void* p, gsl_vector *df)
{
    ML_D_Params* params = (ML_D_Params*)p;
    const field<mat>* Xf = params->Xf;
    const field<vec>* Yf = params->Yf;
    const field<mat>* Zf = params->Zf;
    const vec* beta = params->beta;
    const uword n = params->n;
    const uword q = params->q;
    size_t ntarget = q * (q + 1) / 2;
    vec D_tri(ntarget, arma::fill::zeros);
    for (size_t i = 0; i < ntarget; i++)
    {
        D_tri(i) = gsl_vector_get(v, i);
    }
    mat D(q, q, arma::fill::zeros);
    D(trimatl_ind(size(D))) = D_tri;
    D(trimatu_ind(size(D))) = D_tri;
    vec dL_D = loglikelihood_ml_d(*Xf, *Yf, *Zf, D, *beta, n)(trimatl_ind(size(D)));
    for (uword i = 0; i < ntarget; i++)
    {
        gsl_vector_set(df, i, dL_D(i));
    }
}

void loglikelihood_ml_fd_gsl(const gsl_vector* v, void* p, double *f, gsl_vector *df)
{
    *f = loglikelihood_ml_gsl(v, p);
    loglikelihood_ml_d_gsl(v, p, df);
}

mat fit_D(const field<mat>& Xf, const field<vec>& Yf, const field<mat>& Zf, const mat& D, const vec& beta, const uword& ndata, const double& alpha, const double& eps, const size_t& max_iters)
{
    uword q = D.n_cols, ntarget = q * (q + 1) / 2;
    ML_D_Params* params = new ML_D_Params();
    params->Xf = &Xf;
    params->Yf = &Yf;
    params->Zf = &Zf;
    params->beta = &beta;
    params->n = ndata;
    params->q = q;
    gsl_multimin_function_fdf minex_fun;
    minex_fun.n = ntarget;
    minex_fun.f = loglikelihood_ml_gsl;
    minex_fun.df = loglikelihood_ml_d_gsl;
    minex_fun.fdf = loglikelihood_ml_fd_gsl;
    minex_fun.params = params;
    gsl_vector *target = gsl_vector_alloc(ntarget), *step_size = gsl_vector_alloc(ntarget);
    uvec D_tril_idx = trimatl_ind(arma::size(D)), D_triu_idx = trimatu_ind(arma::size(D));
    vec D_tril_vec = D(D_tril_idx);
    for (uword i = 0; i < ntarget; i++)
    {
        gsl_vector_set(target, i, D_tril_vec(i));
        gsl_vector_set(step_size, i, alpha);
    }
    gsl_multimin_fdfminimizer *minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, ntarget);
    gsl_multimin_fdfminimizer_set(minimizer, &minex_fun, target, alpha, eps);
    size_t iter = 0;
    int status;
    do
    {
        status = gsl_multimin_fdfminimizer_iterate(minimizer);
        if (status) break;
        status = gsl_multimin_test_gradient(minimizer->gradient, eps);
    } while (status == GSL_CONTINUE && (++iter) < max_iters);
    delete params;
    vec D_tri(arma::size(D_tril_idx));
    for (uword i = 0; i < ntarget; i++)
    {
        D_tri(i) = gsl_vector_get(minimizer->x, i);
    }
    mat D1(arma::size(D));
    D1(D_tril_idx) = D_tri;
    D1(D_triu_idx) = D_tri;
    return D1;
}

mat fit_mu(const field<mat>& Xf, const field<vec>& Yf, const field<mat>& Zf, const vec& beta, const mat& D)
{
    uword ngroup = Xf.n_rows, q = Zf(0).n_cols;
    mat D_inv = D.i(), mu(ngroup, q, arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf(i);
        const mat& Yi = Yf(i);
        const mat& Zi = Zf(i);
        uword ndata = Zi.n_rows;
        mat Vi = Zi * D * Zi.t() + eye(ndata, ndata);
        mat Vi_inv = eye(ndata, ndata) - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
        vec Ri = Yi - Xi * beta;
        mu.row(i) = (D * Zi.t() * Vi_inv * Ri).t();
    }
    return mu;
}

HLMGWRParams backfitting_maximum_likelihood(const HLMGWRArgs& args, double alpha, double eps_iter, double eps_gradient, size_t max_iters, bool verbose) 
{
    int prescition = (int)log10(1 / eps_iter);
    //===============
    // Prepare Matrix
    //===============
    const mat& G = args.G;
    const mat& Z = args.Z;
    const mat& X = args.X;
    const vec& y = args.y;
    const mat& u = args.u;
    const uvec& group = args.group;
    double bw = args.bw;
    uword ngroup = G.n_rows, ndata = X.n_rows;
    uword nvg = G.n_cols, nvx = X.n_cols, nvz = Z.n_cols;
    mat gamma(ngroup, nvg, arma::fill::zeros);
    vec beta(nvx, arma::fill::zeros);
    mat mu(ngroup, nvz, arma::fill::zeros);
    mat D(nvz, nvz, arma::fill::eye);
    mat gmap(ndata, ngroup, arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        gmap.col(i) = conv_to<vec>::from(group == i);
    }
    field<mat> Zf(ngroup), Xf(ngroup);
    field<vec> Yf(ngroup), Ygf(ngroup), Yhf(ngroup);
    for (uword i = 0; i < ngroup; i++)
    {
        uvec ind = find(group == i);
        Yf(i) = y.rows(ind);
        Xf(i) = X.rows(ind);
        Zf(i) = Z.rows(ind);
    }
    //============
    // Backfitting
    //============
    double rss = 0.0, rss0 = 0.0, diff = DBL_MAX;
    for (size_t iter = 0; iter < max_iters && diff > eps_iter; iter++)
    {
        rss0 = rss;
        //--------------------
        // Initial Guess for M
        //--------------------
        for (uword i = 0; i < ngroup; i++)
        {
            Ygf(i) = Yf(i) - Xf(i) * beta;
        }
        gamma = fit_gwr(G, Ygf, Zf, D, u, bw);
        vec hatMg = sum(G % gamma, 1);
        vec hatM = hatMg.rows(group);
        vec yh = y - hatM;
        for (uword i = 0; i < ngroup; i++)
        {
            Yhf(i) = Yf(i) - sum(G.row(i) % gamma.row(i));
        }
        //----------------------------------------------
        // Generalized Least Squared Estimation for beta
        //----------------------------------------------
        beta = fit_gls(Xf, Yf, Zf, D);
        //------------------------------------
        // Maximum Likelihood Estimation for D
        //------------------------------------
        D = fit_D(Xf, Yhf, Zf, D, beta, ndata, alpha, eps_gradient, max_iters);
        mu = fit_mu(Xf, Yhf, Zf, beta, D);
        //------------------------------
        // Calculate Termination Measure
        //------------------------------
        vec yhat = yh - (X * beta) - sum(Z % (mu.rows(group)), 1);
        vec residual = yhat % yhat;
        rss = sum(residual);
        diff = abs(rss - rss0);
        if (verbose)
        {
            std::cout << "RSS: " << fixed << setprecision(prescition) << rss << ", diff: " << diff << endl;
        }
    }
    return { gamma, beta, mu, D };
}

int main(int argc, char *argv[])
{
/// Command Line Options
    boost::program_options::options_description desc("Gradient Descent Solution for HLMGWR");
    desc.add_options()
        ("data-dir,d", boost::program_options::value<string>(), "Data directory")
        ("bandwidth,b", boost::program_options::value<double>(), "Bandwidth")
        ("alpha,a", boost::program_options::value<double>()->default_value(0.01), "Learning speed")
        ("eps-iter,e", boost::program_options::value<double>()->default_value(1e-6, "1e-6"), "Coverage threshold")
        ("eps-gradient,g", boost::program_options::value<double>()->default_value(1e-6, "1e-6"), "Minimize Log-likelihood threshold")
        ("maxiters,m", boost::program_options::value<size_t>()->default_value(1e6, "1e6"), "Maximum iteration")
        ("verbose,v", "Print algorithm details.")
        ("help,h", "Print help.");
    boost::program_options::variables_map var_map;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), var_map);
    boost::program_options::notify(var_map);
    double alpha, eps_iter, eps_gradient;
    size_t max_iters;
    bool verbose;
    if (var_map.count("help") > 0)
    {
        cout << desc << endl;
        return 0;
    }
    if (var_map.count("bandwidth") <= 0)
    {
        cout << "Bandwidth must be specified!" << endl;
        return 1;
    }
    string data_dir;
    if (var_map.count("data-dir") > 0) data_dir = var_map["data-dir"].as<string>();
    else 
    {
        cout << "Argument data-dir must be specified!" << endl;
        return 2;
    }
    if (var_map.count("alpha") > 0) alpha = var_map["alpha"].as<double>();
    if (var_map.count("eps-iter") > 0) eps_iter = var_map["eps-iter"].as<double>();
    if (var_map.count("eps-gradient") > 0) eps_gradient = var_map["eps-gradient"].as<double>();
    if (var_map.count("maxiters") > 0) max_iters = var_map["maxiters"].as<size_t>();
    if (var_map.count("verbose") > 0) verbose = true;
    double bw = var_map["bandwidth"].as<double>();
    /// solve
    mat G,X,Z,u;
    vec y;
    uvec group;
    // Read Data
    X.load(arma::csv_name(string(data_dir) + "/hlmgwr_x.csv"));
    G.load(arma::csv_name(string(data_dir) + "/hlmgwr_g.csv"));
    Z.load(arma::csv_name(string(data_dir) + "/hlmgwr_z.csv"));
    u.load(arma::csv_name(string(data_dir) + "/hlmgwr_u.csv"));
    y.load(arma::csv_name(string(data_dir) + "/hlmgwr_y.csv"));
    group.load(arma::csv_name(string(data_dir) + "/hlmgwr_group.csv"));
    HLMGWRArgs alg_args = { G, X, Z, y, u, group, bw };
    HLMGWRParams alg_params = backfitting_maximum_likelihood(alg_args, alpha, eps_iter, eps_gradient, max_iters, verbose);
    // Diagnostic
    const mat &gamma = alg_params.gamma, &beta = alg_params.beta, &mu = alg_params.mu, &D = alg_params.D;
    uword ngroup = G.n_rows, ndata = y.n_rows;
    vec yhat(ndata, arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        uvec ind = find(group == i);
        yhat(ind) = as_scalar(G.row(i) * gamma.row(i).t()) + X.rows(ind) * beta + Z.rows(ind) * mu.row(i).t();
    }
    vec residual = y - yhat;
    vec deviation = y - mean(y);
    double rss = 1 - sum(residual % residual) / sum(deviation % deviation);
    cout << "Rsquared: " << rss << endl;
    // Save coefficients
    alg_params.gamma.save(arma::csv_name(string(data_dir) + "/hlmgwr_hat_gamma.csv"));
    alg_params.beta.save(arma::csv_name(string(data_dir) + "/hlmgwr_hat_beta.csv"));
    alg_params.mu.save(arma::csv_name(string(data_dir) + "/hlmgwr_hat_mu.csv"));
    alg_params.D.save(arma::csv_name(string(data_dir) + "/hlmgwr_hat_D.csv"));
}