#include "hlmgwr.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>
#include <string>
#include <utility>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>

using namespace std;
using namespace arma;
using namespace hgwr;

const double log2pi = log(2.0 * M_PI);

namespace
{

constexpr double D_JITTER = 1e-10;
constexpr double ML_LINE_TOL = 0.1;
constexpr double ML_BAD_OBJECTIVE = 1e100;

mat theta_to_cholesky(const vec& theta, const uword q)
{
    mat L(q, q, arma::fill::zeros);
    uword k = 0;
    for (uword col = 0; col < q; ++col)
    {
        for (uword row = col; row < q; ++row, ++k)
        {
            if (row == col)
            {
                if (!std::isfinite(theta(k)) || theta(k) > 350.0 || theta(k) < -350.0) return mat();
                L(row, col) = std::exp(theta(k));
            }
            else
            {
                L(row, col) = theta(k);
            }
        }
    }
    return L;
}

mat theta_to_covariance(const vec& theta, const uword q)
{
    mat L = theta_to_cholesky(theta, q);
    if (L.is_empty() || !L.is_finite()) return mat();
    mat D = L * L.t() + D_JITTER * eye<mat>(q, q);
    return D.is_finite() ? D : mat();
}

vec covariance_to_theta(const mat& D)
{
    const uword q = D.n_rows;
    mat L;
    mat adjusted = 0.5 * (D + D.t()) - D_JITTER * eye<mat>(q, q);
    if (!chol(L, adjusted, "lower"))
    {
        adjusted = 0.5 * (D + D.t()) + D_JITTER * eye<mat>(q, q);
        if (!chol(L, adjusted, "lower")) L = eye<mat>(q, q);
    }
    vec theta(q * (q + 1) / 2, arma::fill::zeros);
    uword k = 0;
    for (uword col = 0; col < q; ++col)
    {
        for (uword row = col; row < q; ++row, ++k)
        {
            theta(k) = row == col ? std::log(std::max(L(row, col), std::sqrt(D_JITTER))) : L(row, col);
        }
    }
    return theta;
}

bool stable_covariance_inverse(const mat& D, const mat& Z, mat& V_inv, double* log_det_V = nullptr)
{
    mat L;
    if (!D.is_finite() || !chol(L, 0.5 * (D + D.t()), "lower")) return false;
    mat U = Z * L;
    mat middle = eye<mat>(D.n_rows, D.n_cols) + U.t() * U;
    mat middle_chol;
    if (!middle.is_finite() || !chol(middle_chol, middle, "lower")) return false;
    mat solved;
    if (!solve(solved, middle, U.t(), solve_opts::likely_sympd) || !solved.is_finite()) return false;
    V_inv = eye<mat>(Z.n_rows, Z.n_rows) - U * solved;
    V_inv = 0.5 * (V_inv + V_inv.t());
    if (!V_inv.is_finite()) return false;
    if (log_det_V != nullptr)
    {
        *log_det_V = 2.0 * sum(log(middle_chol.diag()));
        if (!std::isfinite(*log_det_V)) return false;
    }
    return true;
}

vec theta_gradient_from_D(const mat& objective_gradient_D, const vec& theta, const uword q)
{
    mat L = theta_to_cholesky(theta, q);
    if (L.is_empty()) return vec();
    mat gradient_L = (objective_gradient_D + objective_gradient_D.t()) * L;
    vec gradient(theta.n_elem, arma::fill::zeros);
    uword k = 0;
    for (uword col = 0; col < q; ++col)
    {
        for (uword row = col; row < q; ++row, ++k)
        {
            gradient(k) = gradient_L(row, col) * (row == col ? L(row, col) : 1.0);
        }
    }
    return gradient;
}

vec gsl_to_arma(const gsl_vector* value, const uword n)
{
    vec result(n);
    for (uword i = 0; i < n; ++i) result(i) = gsl_vector_get(value, i);
    return result;
}

void arma_to_gsl(const vec& value, gsl_vector* result)
{
    for (uword i = 0; i < value.n_elem; ++i) gsl_vector_set(result, i, value(i));
}

double gradient_norm(const gsl_vector* gradient)
{
    double squared = 0.0;
    for (size_t i = 0; i < gradient->size; ++i)
    {
        double value = gsl_vector_get(gradient, i);
        squared += value * value;
    }
    return std::sqrt(squared);
}

} // namespace

double HGWR::bw_criterion_cv(double bw, void* params)
{
    BwSelectionArgs* args = (BwSelectionArgs*)params;
    const mat& Vig = args->Vig.get();
    const vec& Viy = args->Viy.get();
    const mat& G = args->G.get();
    const mat& u = args->u.get();
    const mat* Ygf = args->Ygf;
    const mat* Zf = args->Zf;
    const mat& mu = args->mu.get();
    const size_t ngroup = Viy.n_rows;
    /// Calibrate for each gorup.
    double cv = 0.0;
    for (size_t i = 0; i < ngroup; i++)
    {
        mat d_u = u.each_row() - u.row(i);
        vec d = sqrt(sum(d_u % d_u, 1));
        double b = actual_bw(d, bw);
        vec wW = (*args->kernel)(d % d, b * b);
        wW(i) = 0;
        mat GtWVG = (G.each_col() % wW).t() * Vig;
        mat GtWVy = (G.each_col() % wW).t() * Viy;
        try
        {
            vec gammai = inv(GtWVG) * GtWVy;
            vec hat_ygi = as_scalar(G.row(i) * gammai) + Zf[i] * mu.row(i).t();
            vec residual = Ygf[i] - hat_ygi;
            cv += sum(residual % residual);
        }
        catch(const std::exception& e)
        {
            (*(args->printer))(string("Error occurred when calculating CV value in bandwidth optimisation: ") + e.what() + "\n");
            return DBL_MAX;
        }
    }
    return cv;
}

double HGWR::bw_criterion_aic(double bw, void* params)
{
    BwSelectionArgs* args = (BwSelectionArgs*)params;
    const mat& Vig = args->Vig.get();
    const vec& Viy = args->Viy.get();
    const mat& G = args->G.get();
    const mat& u = args->u.get();
    const mat* Ygf = args->Ygf;
    const mat* Zf = args->Zf;
    const mat& mu = args->mu.get();
    const mat& rVsigma = args->rVsigma.get();
    const uvec& group = args->group.get();
    const size_t ngroup = Viy.n_rows;
    /// Calibrate for each gorup.
    double rss = 0.0;
    double trS = 0.0;
    for (size_t i = 0; i < ngroup; i++)
    {
        mat d_u = u.each_row() - u.row(i);
        vec d = sqrt(sum(d_u % d_u, 1));
        double b = actual_bw(d, bw);
        vec wW = (*args->kernel)(d % d, b * b);
        mat GtW = trans(G.each_col() % wW);
        mat GtWVG = GtW * Vig;
        mat GtWVy = GtW * Viy;
        try
        {
            mat GtWVG_inv = inv(GtWVG);
            vec gammai = GtWVG_inv * GtWVy;
            uvec igroup = find(group == i);
            // mat GtWe = GtW.cols(group);
            // mat si = G.rows(group.rows(find(group == i))) * GtWVG_inv * (GtWe.each_row() % rVsigma);
            mat si_left = G.rows(group.rows(igroup)) * GtWVG_inv * GtW.col(i);
            mat si = si_left * rVsigma.cols(igroup);
            trS += trace(si);
            vec hat_ygi = as_scalar(G.row(i) * gammai) + Zf[i] * mu.row(i).t();
            vec residual = Ygf[i] - hat_ygi;
            rss += sum(residual % residual);
        }
        catch(const std::exception& e)
        {
            (*(args->printer))(string("Error occurred when calculating AIC value in bandwidth optimisation: ") + e.what());
            return DBL_MAX;
        }
    }
    double n = double(ngroup);
    double aic = n * log(rss / n) + n * log(2 * arma::datum::pi) + n + trS;
    return aic;
}

int HGWR::bw_optimisation(double lower, double upper, const BwSelectionArgs* args)
{
    gsl_set_error_handler([](const char* reason, const char* file, int line, int gsl_errno)
    {
        (void)reason;
        (void)file;
        (void)line;
        (void)gsl_errno;
    });
    gsl_function func;
    func.params = (void*)args;
    func.function = bw_criterion;
    gsl_min_fminimizer* minimizer = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    const double R = (sqrt(5)-1)/2;
    double m = lower + R * (upper - lower);
    int status = gsl_min_fminimizer_set(minimizer, &func, m, lower, upper);
    if (status == GSL_EINVAL)
    {
        if (bw > 0)
        {
            if (verbose > 0) pcout("Bandwidth optimisation failed. Use last value: " + to_string(bw) + "\n");
        }
        else
        {
            bw = gsl_min_fminimizer_x_minimum(minimizer);
            if (verbose > 0) pcout("Bandwidth optimisation failed to initialise. Use default value: " + to_string(bw) + "\n");
        }
        return GSL_EINVAL;
    }
    size_t iter = 0;
    do
    {
        status = gsl_min_fminimizer_iterate(minimizer);
        m = gsl_min_fminimizer_x_minimum(minimizer);
        lower = gsl_min_fminimizer_x_lower(minimizer);
        upper = gsl_min_fminimizer_x_upper(minimizer);
        status = gsl_min_test_interval(lower, upper, 1e-4, 0.0);
        if (verbose > 1)
        {
            double fm = gsl_min_fminimizer_f_minimum(minimizer);
            pcout(string("xL: ") + to_string(lower) + "; xU: " + to_string(upper) + "; x: " + to_string(m) + "; f: " + to_string(fm) + "\r");
        }
    } while (status == GSL_CONTINUE && (++iter) < max_bw_iters);
    if (status == GSL_SUCCESS)
    {
        bw = m;
        double fm = gsl_min_fminimizer_f_minimum(minimizer);
        if (verbose > 1) pcout("\n");
        if (verbose > 0) pcout(string("bw: ") + to_string(bw) + "; f: " + to_string(fm) + "\n");
    }
    else
    {
        if (verbose > 0) pcout("Bandwidth optimisation failed. Use last value: " + to_string(bw) + "\n");
    }
    gsl_min_fminimizer_free(minimizer);
    gsl_set_error_handler_off();
    return status;
}

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
void HGWR::fit_gwr(const bool t_test, const bool f_test)
{
    uword k = G.n_cols;//, q = Zf[0].n_cols;
    gamma.fill(arma::fill::zeros);
    if (t_test) gamma_se.fill(arma::fill::zeros);
    unique_ptr<mat[]> Vf = make_unique<mat[]>(ngroup);
    mat Vig(ngroup, k, arma::fill::zeros);
    vec Viy(ngroup, arma::fill::zeros);
    vec Yg(ngroup, arma::fill::zeros);
    rowvec rVsigma = rowvec(ndata, arma::fill::zeros);
    rowvec Vig_var(ngroup, arma::fill::zeros);
    for (size_t i = 0; i < ngroup; i++)
    {
        const mat& Yi = Ygf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv;
        if (!stable_covariance_inverse(D, Zi, Vi_inv)) throw runtime_error("Unable to factor group covariance in GWR fit.");
        uword nidata = Zi.n_rows;
        if (f_test || t_test) Vf[i] = Zi * D * Zi.t() + eye(Zi.n_rows, Zi.n_rows);
        rowvec Visigma = ones(1, nidata) * Vi_inv;
        Vig.row(i) = Visigma * ones(nidata, 1) * G.row(i);
        Viy(i) = as_scalar(Visigma * Yi);
        rVsigma(find(group == i)) = Visigma;
        if (t_test) Vig_var(i) = as_scalar(Visigma * Vf[i] * Visigma.t());
    }
    /// Check whether need to optimize bw
    if (bw_optim)
    {
        BwSelectionArgs args { Vig, Viy, G, u, Ygf.get(), Zf.get(), mu, rVsigma, group, gwr_kernel, Printer };
        if (verbose > 1) {
            args.printer = pcout;
        }
        uword extra = (kernel == KernelType::BISQUARED) ? 1 : 0;
        double upper = double(ngroup - 1), lower = double(k + 2 + extra);
        bw_optimisation(lower, upper, &args);
    }
    /// Calibrate for each gorup.
    trS = { 0.0, 0.0 };
    trQ = { 0.0, 0.0 };
    unique_ptr<mat[]> Qf = make_unique<mat[]>(ngroup);
    for (uword j = 0; j < ngroup; j++)
    {
        Qf[j].resize(size(Vf[j]));
        Qf[j].fill(0.0);
    }
    for (size_t i = 0; i < ngroup; i++)
    {
        mat d_u = u.each_row() - u.row(i);
        vec d = sqrt(sum(d_u % d_u, 1));
        double b = actual_bw(d, bw);
        vec wW = (*gwr_kernel)(d % d, b * b);
        mat GtW = (G.each_col() % wW).t();
        mat GtWVG = GtW * Vig;
        mat GtWVy = GtW * Viy;
        mat GtWVG_inv = inv(GtWVG);
        vec gammai = GtWVG_inv * GtWVy;
        gamma.row(i) = trans(gammai);
        mat Ci = GtWVG_inv * GtW;
        if (t_test) gamma_se.row(i) = trans(sum((Ci.each_row() % Vig_var) % Ci, 1));
        uvec igroup = find(group == i);
        uword nidata = igroup.n_elem;
        // mat GtWe = GtW.cols(group);
        // mat si = G.rows(group.rows(find(group == i))) * GtWVG_inv * (GtWe.each_row() % rVsigma);
        mat si_left = (repelem(G.row(i), nidata, 1) * Ci).eval().cols(group);
        mat si = si_left.each_row() % rVsigma;
        trS(0) += trace(si.cols(igroup));
        trS(1) += trace(si * si.t());
        if (f_test)
        {
            mat ei(nidata, ndata, arma::fill::zeros);
            ei.cols(igroup) = eye(nidata, nidata);
            mat pi = ei - si;
            for (uword j = 0; j < ngroup; j++)
            {
                mat pij = pi.cols(group_span[j]);
                Qf[j] += pij.t() * pij;
            }
        }
        vec hat_ygi = as_scalar(G.row(i) * gammai) + Zf[i] * mu.row(i).t();
    }
    if (f_test)
    {
        for (uword j = 0; j < ngroup; j++)
        {
            Qf[j] *= Vf[j];
            trQ(0) += trace(Qf[j]);
            trQ(1) += trace(Qf[j] * Qf[j]);
        }
    }
    if (t_test)
    {
        gamma_se = sigma * sqrt(gamma_se);
    }
}

vec HGWR::fit_gls()
{
    uword p = Xf[0].n_cols;
    mat XtWX(p, p, arma::fill::zeros);
    vec XtWY(p, arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yhf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv;
        if (!stable_covariance_inverse(D, Zi, Vi_inv)) throw runtime_error("Unable to factor group covariance in GLS fit.");
        XtWX += Xi.t() * Vi_inv * Xi;
        XtWY += Xi.t() * Vi_inv * Yi;
    }
    return solve(XtWX, XtWY);
}

double loglikelihood(const mat* Xf, const vec* Yf, const mat* Zf, const size_t ngroup, const mat& D, const vec& beta, const uword& ndata)
{
    double L1 = 0.0, L2 = 0.0, n = (double)ndata;
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const vec& Yi = Yf[i];
        const mat& Zi = Zf[i];
        double detVi = 0.0;
        mat Vi_inv;
        if (!stable_covariance_inverse(D, Zi, Vi_inv, &detVi))
        {
            return -std::numeric_limits<double>::infinity();
        }
        vec Ri = Yi - Xi * beta;
        L1 += as_scalar(Ri.t() * Vi_inv * Ri);
        L2 += detVi;
    }
    if (!(L1 > 0.0) || !std::isfinite(L1) || !std::isfinite(L2))
    {
        return -std::numeric_limits<double>::infinity();
    }
    double LL = - (n / 2.0) * log(L1) - 0.5 * L2 - 0.5 - 0.5 * log2pi + (n / 2.0) * log(n);
    return LL;
}

void loglikelihood_d(const mat* Xf, const vec* Yf, const mat* Zf, const size_t ngroup, const mat& D, const vec& beta, const uword& ndata, mat& d_D)
{
    mat ZtViZ(arma::size(D), arma::fill::zeros);
    mat KKt(arma::size(D), arma::fill::zeros);
    double J = 0.0, n = (double)ndata;
    // field<mat> Kf(ngroup);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv;
        if (!stable_covariance_inverse(D, Zi, Vi_inv))
        {
            d_D.fill(std::numeric_limits<double>::quiet_NaN());
            return;
        }
        vec Ri = Yi - Xi * beta;
        mat Ki = Zi.t() * Vi_inv * Ri;
        KKt += Ki * Ki.t();
        ZtViZ += Zi.t() * Vi_inv * Zi;
        J += as_scalar(Ri.t() * Vi_inv * Ri);
    }
    if (!(J > 0.0) || !std::isfinite(J))
    {
        d_D.fill(std::numeric_limits<double>::quiet_NaN());
        return;
    }
    mat KJKt = KKt / J;
    d_D = ((- n / 2.0) * (-KJKt) - 0.5 * ZtViZ);
}

void loglikelihood_d(const mat* Xf, const vec* Yf, const mat* Zf, const size_t ngroup, const mat& D, const vec& beta, const uword& ndata, mat& d_D, mat& d_beta)
{
    mat ZtViZ(arma::size(D), arma::fill::zeros);
    mat KKt(arma::size(D), arma::fill::zeros), G(arma::size(beta), arma::fill::zeros);
    double J = 0.0, n = (double)ndata;
    // field<mat> Kf(ngroup);
    field<mat> Kf(ngroup), Gf(ngroup);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv;
        if (!stable_covariance_inverse(D, Zi, Vi_inv))
        {
            d_D.fill(std::numeric_limits<double>::quiet_NaN());
            d_beta.fill(std::numeric_limits<double>::quiet_NaN());
            return;
        }
        vec Ri = Yi - Xi * beta;
        mat Ki = Zi.t() * Vi_inv * Ri;
        KKt += Ki * Ki.t();
        G += Xi.t() * Vi_inv * Ri;
        ZtViZ += Zi.t() * Vi_inv * Zi;
        J += as_scalar(Ri.t() * Vi_inv * Ri);
    }
    if (!(J > 0.0) || !std::isfinite(J))
    {
        d_D.fill(std::numeric_limits<double>::quiet_NaN());
        d_beta.fill(std::numeric_limits<double>::quiet_NaN());
        return;
    }
    mat KJKt = KKt / J;
    mat GJ = G / J;
    d_D = ((- n / 2.0) * (-KJKt) - 0.5 * ZtViZ);
    d_beta = n * GJ;
}

double ml_gsl_f_D(const gsl_vector* v, void* p)
{
    ML_Params* params = (ML_Params*)p;
    const uword ntarget = params->q * (params->q + 1) / 2;
    vec theta = gsl_to_arma(v, ntarget);
    mat D = theta_to_covariance(theta, params->q);
    if (D.is_empty()) return ML_BAD_OBJECTIVE;
    double logL = loglikelihood(params->Xf, params->Yf, params->Zf,
        params->ngroup, D, *params->beta, params->n);
    return std::isfinite(logL) ? -logL / double(params->n) : ML_BAD_OBJECTIVE;
}

double ml_gsl_f_D_beta(const gsl_vector* v, void* pparams)
{
    ML_Params* params = (ML_Params*)pparams;
    const uword ntheta = params->q * (params->q + 1) / 2;
    vec values = gsl_to_arma(v, params->p + ntheta);
    vec beta = values.head(params->p);
    mat D = theta_to_covariance(values.tail(ntheta), params->q);
    if (D.is_empty()) return ML_BAD_OBJECTIVE;
    double logL = loglikelihood(params->Xf, params->Yf, params->Zf,
        params->ngroup, D, beta, params->n);
    return std::isfinite(logL) ? -logL / double(params->n) : ML_BAD_OBJECTIVE;
}

void ml_gsl_df_D(const gsl_vector* v, void* p, gsl_vector *df)
{
    ML_Params* params = (ML_Params*)p;
    const uword ntarget = params->q * (params->q + 1) / 2;
    vec theta = gsl_to_arma(v, ntarget);
    mat D = theta_to_covariance(theta, params->q);
    mat dL_D;
    vec gradient(ntarget, arma::fill::zeros);
    if (!D.is_empty())
    {
        loglikelihood_d(params->Xf, params->Yf, params->Zf, params->ngroup,
            D, *params->beta, params->n, dL_D);
        if (dL_D.is_finite())
        {
            gradient = theta_gradient_from_D(-dL_D / double(params->n), theta, params->q);
        }
    }
    arma_to_gsl(gradient, df);
}

void ml_gsl_df_D_beta(const gsl_vector* v, void* pparams, gsl_vector *df)
{
    ML_Params* params = (ML_Params*)pparams;
    const uword ntheta = params->q * (params->q + 1) / 2;
    vec values = gsl_to_arma(v, params->p + ntheta);
    vec beta = values.head(params->p), theta = values.tail(ntheta);
    mat D = theta_to_covariance(theta, params->q);
    mat dL_D;
    vec dL_beta;
    vec gradient(params->p + ntheta, arma::fill::zeros);
    if (!D.is_empty())
    {
        loglikelihood_d(params->Xf, params->Yf, params->Zf, params->ngroup,
            D, beta, params->n, dL_D, dL_beta);
        if (dL_D.is_finite() && dL_beta.is_finite())
        {
            gradient.head(params->p) = -dL_beta / double(params->n);
            gradient.tail(ntheta) = theta_gradient_from_D(
                -dL_D / double(params->n), theta, params->q);
        }
    }
    arma_to_gsl(gradient, df);
}

void ml_gsl_fdf_D(const gsl_vector* v, void* p, double *f, gsl_vector *df)
{
    *f = ml_gsl_f_D(v, p);
    ml_gsl_df_D(v, p, df);
}

void ml_gsl_fdf_D_beta(const gsl_vector* v, void* p, double *f, gsl_vector *df)
{
    *f = ml_gsl_f_D_beta(v, p);
    ml_gsl_df_D_beta(v, p, df);
}

double hgwr::ml_gradient_relative_error(
    const mat& X,
    const mat& Z,
    const vec& y,
    const uvec& group,
    const vec& beta,
    const mat& D,
    bool include_beta,
    double step
)
{
    uvec normalized_group = group - group.min();
    const uword ngroup = normalized_group.max() + 1;
    unique_ptr<mat[]> Xf = make_unique<mat[]>(ngroup);
    unique_ptr<mat[]> Zf = make_unique<mat[]>(ngroup);
    unique_ptr<vec[]> Yf = make_unique<vec[]>(ngroup);
    for (uword i = 0; i < ngroup; ++i)
    {
        uvec index = find(normalized_group == i);
        Xf[i] = X.rows(index);
        Zf[i] = Z.rows(index);
        Yf[i] = y.rows(index);
    }
    vec beta_copy = beta;
    ML_Params params {
        Xf.get(), Yf.get(), Zf.get(), &beta_copy, ngroup,
        X.n_rows, X.n_cols, Z.n_cols
    };
    vec theta = covariance_to_theta(D);
    vec target_values = include_beta ? join_cols(beta, theta) : theta;
    gsl_vector* target = gsl_vector_alloc(target_values.n_elem);
    gsl_vector* analytic_gsl = gsl_vector_alloc(target_values.n_elem);
    arma_to_gsl(target_values, target);
    if (include_beta) ml_gsl_df_D_beta(target, &params, analytic_gsl);
    else ml_gsl_df_D(target, &params, analytic_gsl);
    vec analytic = gsl_to_arma(analytic_gsl, target_values.n_elem);
    vec numeric(target_values.n_elem, arma::fill::zeros);
    for (uword i = 0; i < target_values.n_elem; ++i)
    {
        double original = gsl_vector_get(target, i);
        gsl_vector_set(target, i, original + step);
        double upper = include_beta ? ml_gsl_f_D_beta(target, &params) : ml_gsl_f_D(target, &params);
        gsl_vector_set(target, i, original - step);
        double lower = include_beta ? ml_gsl_f_D_beta(target, &params) : ml_gsl_f_D(target, &params);
        gsl_vector_set(target, i, original);
        numeric(i) = (upper - lower) / (2.0 * step);
    }
    gsl_vector_free(analytic_gsl);
    gsl_vector_free(target);
    vec denominator = arma::max(ones<vec>(analytic.n_elem), arma::max(abs(analytic), abs(numeric)));
    return max(abs(analytic - numeric) / denominator);
}

double HGWR::fit_D(ML_Params* params)
{
    const uword ntarget = D.n_cols * (D.n_cols + 1) / 2;
    vec start = covariance_to_theta(D), best = start;
    gsl_vector* target = gsl_vector_alloc(ntarget);
    arma_to_gsl(start, target);
    gsl_multimin_function_fdf function;
    function.f = ml_gsl_f_D;
    function.df = ml_gsl_df_D;
    function.fdf = ml_gsl_fdf_D;
    function.n = ntarget;
    function.params = params;
    gsl_multimin_fdfminimizer* minimizer = gsl_multimin_fdfminimizer_alloc(
        gsl_multimin_fdfminimizer_vector_bfgs2, ntarget);
    gsl_set_error_handler_off();
    int status = gsl_multimin_fdfminimizer_set(minimizer, &function, target, alpha, ML_LINE_TOL);
    double initial = status == GSL_SUCCESS && std::isfinite(minimizer->f)
        ? minimizer->f : ML_BAD_OBJECTIVE;
    double best_objective = initial;
    size_t iter = 0;
    while (status == GSL_SUCCESS && iter < max_iters)
    {
        status = gsl_multimin_fdfminimizer_iterate(minimizer);
        ++iter;
        if (status != GSL_SUCCESS) break;
        if (std::isfinite(minimizer->f) && minimizer->f < best_objective)
        {
            best_objective = minimizer->f;
            best = gsl_to_arma(minimizer->x, ntarget);
        }
        if (verbose > 1)
        {
            pcout("ML iter: " + to_string(iter) + ", objective: " + to_string(minimizer->f)
                + ", gradient: " + to_string(gradient_norm(minimizer->gradient)) + "\n");
        }
        int gradient_status = gsl_multimin_test_gradient(minimizer->gradient, eps_gradient);
        if (gradient_status == GSL_SUCCESS)
        {
            status = GSL_SUCCESS;
            break;
        }
        if (gradient_status != GSL_CONTINUE)
        {
            status = gradient_status;
            break;
        }
    }
    if (iter >= max_iters && status == GSL_SUCCESS) status = GSL_CONTINUE;
    mat candidate = theta_to_covariance(best, D.n_cols);
    if (!candidate.is_empty() && std::isfinite(best_objective) && best_objective <= initial + 1e-12)
    {
        D = candidate;
    }
    ml_iterations += iter;
    ml_status = status;
    ml_converged = status == GSL_SUCCESS;
    if (!ml_converged) ++ml_failures;
    gsl_multimin_fdfminimizer_free(minimizer);
    gsl_vector_free(target);
    return best_objective;
}

double HGWR::fit_D_beta(ML_Params* params)
{
    const uword ntheta = D.n_cols * (D.n_cols + 1) / 2;
    const uword ntarget = beta.n_elem + ntheta;
    vec start = join_cols(beta, covariance_to_theta(D)), best = start;
    gsl_vector* target = gsl_vector_alloc(ntarget);
    arma_to_gsl(start, target);
    gsl_multimin_function_fdf function;
    function.f = ml_gsl_f_D_beta;
    function.df = ml_gsl_df_D_beta;
    function.fdf = ml_gsl_fdf_D_beta;
    function.n = ntarget;
    function.params = params;
    gsl_multimin_fdfminimizer* minimizer = gsl_multimin_fdfminimizer_alloc(
        gsl_multimin_fdfminimizer_vector_bfgs2, ntarget);
    gsl_set_error_handler_off();
    int status = gsl_multimin_fdfminimizer_set(minimizer, &function, target, alpha, ML_LINE_TOL);
    double initial = status == GSL_SUCCESS && std::isfinite(minimizer->f)
        ? minimizer->f : ML_BAD_OBJECTIVE;
    double best_objective = initial;
    size_t iter = 0;
    while (status == GSL_SUCCESS && iter < max_iters)
    {
        status = gsl_multimin_fdfminimizer_iterate(minimizer);
        ++iter;
        if (status != GSL_SUCCESS) break;
        if (std::isfinite(minimizer->f) && minimizer->f < best_objective)
        {
            best_objective = minimizer->f;
            best = gsl_to_arma(minimizer->x, ntarget);
        }
        if (verbose > 1)
        {
            pcout("ML iter: " + to_string(iter) + ", objective: " + to_string(minimizer->f)
                + ", gradient: " + to_string(gradient_norm(minimizer->gradient)) + "\n");
        }
        int gradient_status = gsl_multimin_test_gradient(minimizer->gradient, eps_gradient);
        if (gradient_status == GSL_SUCCESS)
        {
            status = GSL_SUCCESS;
            break;
        }
        if (gradient_status != GSL_CONTINUE)
        {
            status = gradient_status;
            break;
        }
    }
    if (iter >= max_iters && status == GSL_SUCCESS) status = GSL_CONTINUE;
    vec candidate_beta = best.head(beta.n_elem);
    mat candidate_D = theta_to_covariance(best.tail(ntheta), D.n_cols);
    if (!candidate_D.is_empty() && candidate_beta.is_finite() && std::isfinite(best_objective)
        && best_objective <= initial + 1e-12)
    {
        D = candidate_D;
        beta = candidate_beta;
    }
    ml_iterations += iter;
    ml_status = status;
    ml_converged = status == GSL_SUCCESS;
    if (!ml_converged) ++ml_failures;
    gsl_multimin_fdfminimizer_free(minimizer);
    gsl_vector_free(target);
    return best_objective;
}

void HGWR::fit_mu()
{
    mu.fill(arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yhf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv;
        if (!stable_covariance_inverse(D, Zi, Vi_inv)) throw runtime_error("Unable to factor group covariance in random-effect fit.");
        vec Ri = Yi - Xi * beta;
        mu.row(i) = (D * Zi.t() * Vi_inv * Ri).t();
    }
}

double HGWR::fit_sigma()
{
    double sigma2 = 0.0;
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yhf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv;
        if (!stable_covariance_inverse(D, Zi, Vi_inv)) throw runtime_error("Unable to factor group covariance in residual-scale fit.");
        mat Ri = Yi - Xi * beta;
        sigma2 += as_scalar(Ri.t() * Vi_inv * Ri);
    }
    return sqrt(sigma2 / (double)ndata);
}

HGWR::Parameters HGWR::fit(const bool f_test)
{
    //===============
    // Prepare Matrix
    //===============
    int prescition = (int)log10(1 / eps_iter);
    double tss = sum((y - mean(y)) % (y - mean(y)));
    gamma = mat(ngroup, nvg, arma::fill::zeros);
    gamma_se = mat(ngroup, nvg, arma::fill::zeros);
    beta = vec(nvx, arma::fill::zeros);
    mu = mat(ngroup, nvz, arma::fill::zeros);
    D = mat(nvz, nvz, arma::fill::eye);
    ml_converged = false;
    ml_status = GSL_CONTINUE;
    ml_iterations = 0;
    ml_failures = 0;
    Zf = make_unique<arma::mat[]>(ngroup);
    Xf = make_unique<arma::mat[]>(ngroup);
    Yf = make_unique<arma::vec[]>(ngroup);
    Ygf = make_unique<arma::vec[]>(ngroup);
    Yhf = make_unique<arma::vec[]>(ngroup);
    uvec group_size(ngroup);
    group_span.resize(ngroup);
    for (uword i = 0; i < ngroup; i++)
    {
        uvec ind = find(group == i);
        group_size(i) = ind.n_elem;
        Yf[i] = y.rows(ind);
        Xf[i] = X.rows(ind);
        Zf[i] = Z.rows(ind);
        Yhf[i] = y.rows(ind);
    }
    uvec group_to = cumsum(group_size);
    uvec group_from = group_to - group_size;
    group_to = group_to - 1;
    transform(group_from.begin(), group_from.end(), group_to.begin(), group_span.begin(), [](uword from, uword to)
    {
        return span(from, to);
    });
    //----------------------------------------------
    // Generalized Least Squared Estimation for beta
    //----------------------------------------------
    beta = fit_gls();
    fit_mu();
    //============
    // Backfitting
    //============
    size_t retry = 0, iterations = 0;
    double rss = DBL_MAX, rss0 = DBL_MAX, diff = DBL_MAX;
    double mlf = DBL_MAX, mlf0 = DBL_MAX, relative_objective_diff = DBL_MAX;
    for (size_t iter = 0; relative_objective_diff > eps_iter && iter < max_iters && retry < max_retries; iter++)
    {
        iterations = iter + 1;
        rss0 = rss;
        //--------------------
        // Initial Guess for M
        //--------------------
        for (uword i = 0; i < ngroup; i++)
        {
            Ygf[i] = Yf[i] - Xf[i] * beta;
        }
        fit_gwr();
        vec hatMg = sum(G % gamma, 1);
        vec hatM = hatMg.rows(group);
        vec yh = y - hatM;
        for (uword i = 0; i < ngroup; i++)
        {
            Yhf[i] = Yf[i] - sum(G.row(i) % gamma.row(i));
        }
        //------------------------------------
        // Maximum Likelihood Estimation for D
        //------------------------------------
        mlf0 = mlf;
        ML_Params ml_params = { Xf.get(), Yhf.get(), Zf.get(), &beta, ngroup, ndata, nvx, nvz };
        switch (ml_type)
        {
        case 1:
            ml_params.beta = nullptr;
            beta = fit_gls();
            mlf = fit_D_beta(&ml_params);
            break;
        default:
            mlf = fit_D(&ml_params);
            beta = fit_gls();
            break;
        }
        fit_mu();
        //------------------------------
        // Calculate Termination Measure
        //------------------------------
        vec yhat = yh - (X * beta) - sum(Z % (mu.rows(group)), 1);
        vec residual = yhat % yhat;
        rss = sum(residual);
        diff = rss - rss0;
        relative_objective_diff = std::abs(mlf - mlf0) / std::max(std::abs(mlf0), 1.0);
        if (mlf < mlf0)
        {
            if (retry > 0) retry = 0;
        }
        else if (iter > 0 && relative_objective_diff > std::max(std::sqrt(eps_iter), 1e-4)) retry++;
        if (verbose > 0)
        {
            ostringstream sout;
            sout << fixed << setprecision(prescition) << "Iter: " << iter;
            if (bw_optim) sout << ", " << "Bw: " << bw;
            sout << ", " << "RSS: " << rss;
            if (abs(diff) < DBL_MAX) sout << ", " << "dRSS: " << diff;
            sout << ", " << "R2: " << (1 - rss / tss);
            sout << ", " << "-loglik/n: " << mlf;
            if (std::isfinite(mlf0)) sout << ", " << "relative objective change: " << relative_objective_diff;
            if (retry > 0) sout << ", " << "Retry: " << retry;
            sout << endl;
            pcout(sout.str());
        }
        (*(this->pcancel))();
    }
    sigma = fit_sigma();
    if (verbose > 0) pcout("Re-fit GLSW effects for f test\n");
    for (uword i = 0; i < ngroup; i++)
    {
        Ygf[i] = Yf[i] - Xf[i] * beta;
    }
    fit_gwr(true, f_test);
    //============
    // Diagnostic
    //============
    loglik = - mlf * double(ndata);
    calc_var_beta();
    vec D_eigenvalues;
    const bool eig_ok = eig_sym(D_eigenvalues, 0.5 * (D + D.t()));
    const double min_eigen_D = eig_ok ? D_eigenvalues.min() : datum::nan;
    const bool parameters_finite = gamma.is_finite() && beta.is_finite() && mu.is_finite()
        && D.is_finite() && std::isfinite(sigma) && std::isfinite(loglik)
        && std::isfinite(min_eigen_D) && min_eigen_D > 0.0;
    const bool converged = std::isfinite(relative_objective_diff) && relative_objective_diff <= eps_iter
        && ml_converged && parameters_finite;
    return { gamma, beta, mu, D, sigma, bw, iterations, retry, converged,
        ml_converged, ml_status, ml_iterations, ml_failures, min_eigen_D };
}

void HGWR::calc_var_beta()
{
    mat XtViX(X.n_cols, X.n_cols, arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv;
        if (!stable_covariance_inverse(D, Zi, Vi_inv)) throw runtime_error("Unable to factor group covariance in fixed-effect variance fit.");
        XtViX += Xi.t() * Vi_inv * Xi;
    }
    var_beta = diagvec(XtViX.i());
}

std::vector<arma::vec4> HGWR::test_glsw()
{
    if (verbose > 0) pcout("Preparing f test\n");
    uword ng = gamma.n_cols;
    double nd = double(ndata);
    double df2 = trQ(0) * trQ(0) / trQ(1);
    unique_ptr<mat[]> Vf = make_unique<mat[]>(ngroup);
    unique_ptr<mat[]> GVGf = make_unique<mat[]>(ngroup);
    unique_ptr<mat[]> GVf = make_unique<mat[]>(ngroup);
    for (size_t i = 0; i < ngroup; i++)
    {
        uvec ind = find(group == i);
        const mat& Zi = Zf[i];
        Vf[i] = Zi * D * Zi.t() + eye(Zi.n_rows, Zi.n_rows);
        mat Vi_inv;
        if (!stable_covariance_inverse(D, Zi, Vi_inv)) throw runtime_error("Unable to factor group covariance in GLSW test.");
        uword nidata = Zi.n_rows;
        GVf[i] = G.row(i).t() * ones(1, nidata) * Vi_inv;
        GVGf[i] = (GVf[i] * ones(nidata, 1) * G.row(i));
    }
    vec nw(ngroup, arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        nw(i) = double(Zf[i].n_rows);
    }
    vector<vec4> results;
    for (uword k = 0; k < ng; k++)
    {
        if (verbose > 0) pcout("Doing f test for effect " + to_string(k) + "\n");
        double sum_gk = sum(gamma.col(k) % nw);
        double sum_gk2 = sum(gamma.col(k) % gamma.col(k) % nw);
        double vk2 = (sum_gk2 - sum_gk * sum_gk / nd) / nd;
        vec c(ndata, arma::fill::zeros);
        for (uword i = 0; i < ngroup; i++)
        {
            double ni = double(GVf[i].n_cols);
            mat d_u = u.each_row() - u.row(i);
            vec d = sqrt(sum(d_u % d_u, 1));
            double fb = actual_bw(d, bw);
            vec w = (*gwr_kernel)(d % d, fb * fb);
            mat GWVG(ng, ng, arma::fill::zeros), GWV(ng, ndata, arma::fill::zeros);
            for (size_t j = 0; j < ngroup; j++)
            {
                GWVG += (w[j] * GVGf[j]);
                GWV.cols(find(group == j)) = w[j] * GVf[j];
            }
            mat Cit = GWV.t() * GWVG.i().t();
            vec bi = Cit.col(k);
            c += bi * double(ni);
        }
        unique_ptr<mat[]> Bf = make_unique<mat[]>(ngroup);
        for (size_t j = 0; j < ngroup; j++)
        {
            Bf[j].resize(size(Vf[j]));
            Bf[j].fill(0.0);
        }
        for (uword i = 0; i < ngroup; i++)
        {
            double ni = double(GVf[i].n_cols);
            mat d_u = u.each_row() - u.row(i);
            vec d = sqrt(sum(d_u % d_u, 1));
            double fb = actual_bw(d, bw);
            vec w = (*gwr_kernel)(d % d, fb * fb);
            mat GWVG(ng, ng, arma::fill::zeros), GWV(ng, ndata, arma::fill::zeros);
            for (size_t j = 0; j < ngroup; j++)
            {
                GWVG += GVGf[j] * w[j];
                GWV.cols(find(group == j)) = GVf[j] * w[j];
            }
            mat Cit = GWV.t() * GWVG.i().t();
            vec bi = Cit.col(k);
            for (uword j = 0; j < ngroup; j++)
            {
                vec bij = bi.rows(group_span[j]);
                vec cij = c.rows(group_span[j]);
                Bf[j] += bij * bij.t() * ni - cij * bij.t() * ni / nd;
            }
        }
        double trB = 0.0, trB2 = 0.0;
        for (size_t j = 0; j < ngroup; j++)
        {
            Bf[j] *= Vf[j] / nd;
            trB += trace(Bf[j]);
            trB2 += trace(Bf[j] * Bf[j]);
        }
        double fv = vk2 / trB / (sigma * sigma);
        double df1 = trB * trB / trB2;
        double pv = gsl_cdf_fdist_Q(fv, df1, df2);
        vec4 result = { fv, df1, df2, pv };
        results.push_back(result);
    }
    return results;
}
