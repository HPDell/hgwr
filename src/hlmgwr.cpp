#include "hlmgwr.h"
#include <sstream>
#include <iomanip>
#include <string>
#include <utility>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace arma;
using namespace hgwr;

const double log2pi = log(2.0 * M_PI);

/**
 * @brief Leave-one-group-out CV score for adaptive GWR bandwidth selection.
 *
 * Uses cached group-to-group distances from `BwSelectionArgs` and the current
 * transformed response summaries.
 *
 * @param bw Candidate adaptive bandwidth in nearest-neighbour units.
 * @param params Pointer to `HGWR::BwSelectionArgs`.
 * @return Cross-validation residual sum of squares.
 */
double HGWR::bw_criterion_cv(double bw, void* params)
{
    BwSelectionArgs* args = (BwSelectionArgs*)params;
    const mat& Vig = args->Vig.get();
    const vec& Viy = args->Viy.get();
    const mat& G = args->G.get();
    const mat& distance = args->distance.get();
    const mat& distance2 = args->distance2.get();
    const mat* Ygf = args->Ygf;
    const mat* Zf = args->Zf;
    const mat& mu = args->mu.get();
    const size_t ngroup = Viy.n_rows;
    /// Calibrate for each gorup.
    double cv = 0.0;
    for (size_t i = 0; i < ngroup; i++)
    {
        vec d = distance.col(i);
        double b = actual_bw(d, bw);
        vec wW = (*args->kernel)(distance2.col(i), b * b);
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

/**
 * @brief AIC-style score for adaptive GWR bandwidth selection.
 *
 * Uses cached group-to-group distances and an approximate trace of the local
 * smoother matrix.
 *
 * @param bw Candidate adaptive bandwidth in nearest-neighbour units.
 * @param params Pointer to `HGWR::BwSelectionArgs`.
 * @return AIC-style criterion value.
 */
double HGWR::bw_criterion_aic(double bw, void* params)
{
    BwSelectionArgs* args = (BwSelectionArgs*)params;
    const mat& Vig = args->Vig.get();
    const vec& Viy = args->Viy.get();
    const mat& G = args->G.get();
    const mat& distance = args->distance.get();
    const mat& distance2 = args->distance2.get();
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
        vec d = distance.col(i);
        double b = actual_bw(d, bw);
        vec wW = (*args->kernel)(distance2.col(i), b * b);
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

/**
 * @brief Select an adaptive bandwidth by minimizing the configured criterion.
 *
 * The bounds and candidate values are expressed as nearest-neighbour counts.
 * The search fails explicitly when no finite candidate value is available.
 *
 * @param lower Lower nearest-neighbour bound.
 * @param upper Upper nearest-neighbour bound.
 * @param args Reusable summaries and cached distances for the criterion.
 * @return GSL status code.
 */
int HGWR::bw_optimisation(double lower, double upper, const BwSelectionArgs* args)
{
    const int first = static_cast<int>(std::ceil(lower));
    const int last = static_cast<int>(std::floor(upper));
    double best_value = DBL_MAX;
    int best_bandwidth = -1;
    bw_evaluations = 0;

    for (int candidate = first; candidate <= last; ++candidate)
    {
        double value = bw_criterion(double(candidate), const_cast<BwSelectionArgs*>(args));
        ++bw_evaluations;
        if (verbose > 1)
        {
            pcout(
                string("bw: ") + to_string(candidate) +
                "; criterion: " + to_string(value) + "\n"
            );
        }
        if (std::isfinite(value) && value < best_value)
        {
            best_value = value;
            best_bandwidth = candidate;
        }
    }

    if (best_bandwidth < 0)
    {
        bw_objective = arma::datum::nan;
        if (verbose > 0) pcout("Bandwidth grid search failed: no finite criterion value.\n");
        return GSL_EFAILED;
    }

    bw = double(best_bandwidth);
    bw_objective = best_value;
    if (verbose > 0)
    {
        pcout(
            string("bw: ") + to_string(bw) +
            "; criterion: " + to_string(best_value) + "\n"
        );
    }
    return GSL_SUCCESS;
}

/**
 * @brief Precompute pairwise group distance matrices.
 *
 * `distance` stores Euclidean distances and `distance2` stores squared
 * distances between group coordinates. Both are reused by GWR fitting and
 * bandwidth selection. OpenMP parallelism is enabled only when
 * `ENABLE_OPENMP` is defined by the build system.
 */
void HGWR::ensure_distance_cache()
{
    if (distance_cache_ready &&
        distance.n_rows == ngroup && distance.n_cols == ngroup &&
        distance2.n_rows == ngroup && distance2.n_cols == ngroup)
    {
        return;
    }

    distance.set_size(ngroup, ngroup);
    distance2.set_size(ngroup, ngroup);

#ifdef ENABLE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long i = 0; i < static_cast<long long>(ngroup); i++)
    {
        mat d_u = u.each_row() - u.row(static_cast<uword>(i));
        vec d2 = sum(d_u % d_u, 1);
        distance2.col(static_cast<uword>(i)) = d2;
        distance.col(static_cast<uword>(i)) = sqrt(d2);
    }

    distance_cache_ready = true;
}

/**
 * @brief Estimate group-level spatially weighted effects.
 *
 * This step fits one local GLSW coefficient vector per group using the current
 * fixed/random-effect residuals, random-effect covariance, kernel, and
 * adaptive bandwidth. Optional diagnostic paths calculate standard errors and
 * F-test matrices.
 *
 * @param t_test If true, calculate standard errors for GLSW effects.
 * @param f_test If true, calculate matrices needed by GLSW F tests.
 */
void HGWR::fit_gwr(const bool t_test, const bool f_test)
{
    ensure_distance_cache();
    uword k = G.n_cols;//, q = Zf[0].n_cols;
    mat D_inv = D.i();
    gamma.fill(arma::fill::zeros);
    if (t_test) gamma_se.fill(arma::fill::zeros);
    unique_ptr<mat[]> Vf;
    if (f_test || t_test) Vf = make_unique<mat[]>(ngroup);
    mat Vig(ngroup, k, arma::fill::zeros);
    vec Viy(ngroup, arma::fill::zeros);
    rowvec rVsigma = rowvec(ndata, arma::fill::zeros);
    rowvec Vig_var(ngroup, arma::fill::zeros);
#ifdef ENABLE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long i0 = 0; i0 < static_cast<long long>(ngroup); i0++)
    {
        uword i = static_cast<uword>(i0);
        const mat& Yi = Ygf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv = woodbury_eye(D_inv, Zi);
        uword nidata = Zi.n_rows;
        if (f_test || t_test) Vf[i] = Zi * D * Zi.t() + eye(Zi.n_rows, Zi.n_rows);
        rowvec Visigma = ones(1, nidata) * Vi_inv;
        Vig.row(i) = Visigma * ones(nidata, 1) * G.row(i);
        Viy(i) = as_scalar(Visigma * Yi);
        rVsigma(group_span[i]) = Visigma;
        if (t_test) Vig_var(i) = as_scalar(Visigma * Vf[i] * Visigma.t());
    }
    /// Check whether need to optimize bw
    if (bw_optim)
    {
        BwSelectionArgs args { Vig, Viy, G, u, distance, distance2, Ygf.get(), Zf.get(), mu, rVsigma, group, gwr_kernel, Printer };
        if (verbose > 1) {
            args.printer = pcout;
        }
        uword extra = (kernel == KernelType::BISQUARED) ? 1 : 0;
        double upper = double(ngroup - 1), lower = double(k + 2 + extra);
        bw_lower = lower;
        bw_upper = upper;
        bw_optimizer_status = bw_optimisation(lower, upper, &args);
        if (bw_optimizer_status != GSL_SUCCESS) ++bw_optimizer_failures;
    }
    /// Calibrate for each gorup.
    trS = { 0.0, 0.0 };
    trQ = { 0.0, 0.0 };
    unique_ptr<mat[]> Qf;
    if (f_test)
    {
        Qf = make_unique<mat[]>(ngroup);
        for (uword j = 0; j < ngroup; j++)
        {
            Qf[j].resize(size(Vf[j]));
            Qf[j].fill(0.0);
        }
    }
    for (size_t i = 0; i < ngroup; i++)
    {
        vec d = distance.col(i);
        double b = actual_bw(d, bw);
        vec wW = (*gwr_kernel)(distance2.col(i), b * b);
        mat GtW = (G.each_col() % wW).t();
        mat GtWVG = GtW * Vig;
        mat GtWVy = GtW * Viy;
        mat GtWVG_inv = inv(GtWVG);
        vec gammai = GtWVG_inv * GtWVy;
        gamma.row(i) = trans(gammai);
        mat Ci = GtWVG_inv * GtW;
        if (t_test) gamma_se.row(i) = trans(sum((Ci.each_row() % Vig_var) % Ci, 1));
        span igroup = group_span[i];
        uword nidata = Zf[i].n_rows;
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
    mat D_inv = D.i();
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yhf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv = woodbury_eye(D_inv, Zi);
        XtWX += Xi.t() * Vi_inv * Xi;
        XtWY += Xi.t() * Vi_inv * Yi;
    }
    return solve(XtWX, XtWY);
}

double loglikelihood(const mat* Xf, const vec* Yf, const mat* Zf, const size_t ngroup, const mat& D, const vec& beta, const uword& ndata)
{
    mat D_inv = D.i();
    double L1 = 0.0, L2 = 0.0, n = (double)ndata;
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const vec& Yi = Yf[i];
        const mat& Zi = Zf[i];
        mat Vi = ((Zi * D) * Zi.t()) + eye<mat>(Zi.n_rows, Zi.n_rows);
        double detVi, sign_detVi;
        log_det(detVi, sign_detVi, Vi);
        mat Vi_inv = HGWR::woodbury_eye(D_inv, Zi);
        vec Ri = Yi - Xi * beta;
        L1 += as_scalar(Ri.t() * Vi_inv * Ri);
        L2 += detVi;
    }
    double LL = - (n / 2.0) * log(L1) - 0.5 * L2 - 0.5 - 0.5 * log2pi + (n / 2.0) * log(n);
    return LL;
}

void loglikelihood_d(const mat* Xf, const vec* Yf, const mat* Zf, const size_t ngroup, const mat& D, const vec& beta, const uword& ndata, mat& d_D)
{
    mat ZtViZ(arma::size(D), arma::fill::zeros), D_inv = D.i();
    mat KKt(arma::size(D), arma::fill::zeros);
    double J = 0.0, n = (double)ndata;
    // field<mat> Kf(ngroup);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv = HGWR::woodbury_eye(D_inv, Zi);
        vec Ri = Yi - Xi * beta;
        mat Ki = Zi.t() * Vi_inv * Ri;
        KKt += Ki * Ki.t();
        ZtViZ += Zi.t() * Vi_inv * Zi;
        J += as_scalar(Ri.t() * Vi_inv * Ri);
    }
    mat KJKt = KKt / J;
    d_D = ((- n / 2.0) * (-KJKt) - 0.5 * ZtViZ);
}

void loglikelihood_d(const mat* Xf, const vec* Yf, const mat* Zf, const size_t ngroup, const mat& D, const vec& beta, const uword& ndata, mat& d_D, mat& d_beta)
{
    mat ZtViZ(arma::size(D), arma::fill::zeros), D_inv = D.i();
    mat KKt(arma::size(D), arma::fill::zeros), G(arma::size(beta), arma::fill::zeros);
    double J = 0.0, n = (double)ndata;
    // field<mat> Kf(ngroup);
    field<mat> Kf(ngroup), Gf(ngroup);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yf[i];
        const mat& Zi = Zf[i];
        mat Vi_inv = HGWR::woodbury_eye(D_inv, Zi);
        vec Ri = Yi - Xi * beta;
        mat Ki = Zi.t() * Vi_inv * Ri;
        KKt += Ki * Ki.t();
        G += Xi.t() * Vi_inv * Ri;
        ZtViZ += Zi.t() * Vi_inv * Zi;
        J += as_scalar(Ri.t() * Vi_inv * Ri);
    }
    mat KJKt = KKt / J;
    mat GJ = G / J;
    d_D = ((- n / 2.0) * (-KJKt) - 0.5 * ZtViZ);
    d_beta = n * GJ;
}

double ml_gsl_f_D(const gsl_vector* v, void* p)
{
    ML_Params* params = (ML_Params*)p;
    const mat* Xf = params->Xf;
    const vec* Yf = params->Yf;
    const mat* Zf = params->Zf;
    const vec* beta = params->beta;
    const size_t ngroup = params->ngroup;
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
    D = D.t();
    D(trimatl_ind(size(D))) = D_tri;
    double logL = loglikelihood(Xf, Yf, Zf, ngroup, D, *beta, n);
    return -logL / double(n);
}

double ml_gsl_f_D_beta(const gsl_vector* v, void* pparams)
{
    ML_Params* params = (ML_Params*)pparams;
    const mat* Xf = params->Xf;
    const vec* Yf = params->Yf;
    const mat* Zf = params->Zf;
    const size_t ngroup = params->ngroup;
    const uword n = params->n;
    const uword p = params->p;
    const uword q = params->q;
    size_t ntarget = p + q * (q + 1) / 2;
    vec D_tri(q * (q + 1) / 2, arma::fill::zeros), beta(p, arma::fill::zeros);
    for (size_t i = 0; i < p; i++)
    {
        beta(i) = gsl_vector_get(v, i);
    }
    for (size_t i = p; i < ntarget; i++)
    {
        D_tri(i - p) = gsl_vector_get(v, i);
    }
    mat D(q, q, arma::fill::zeros);
    D(trimatl_ind(size(D))) = D_tri;
    D = D.t();
    D(trimatl_ind(size(D))) = D_tri;
    double logL = loglikelihood(Xf, Yf, Zf, ngroup, D, beta, n);
    return -logL / double(n);
}

void ml_gsl_df_D(const gsl_vector* v, void* p, gsl_vector *df)
{
    ML_Params* params = (ML_Params*)p;
    const mat* Xf = params->Xf;
    const vec* Yf = params->Yf;
    const mat* Zf = params->Zf;
    const vec* beta = params->beta;
    const size_t ngroup = params->ngroup;
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
    D = D.t();
    D(trimatl_ind(size(D))) = D_tri;
    mat dL_D;
    loglikelihood_d(Xf, Yf, Zf, ngroup, D, *beta, n, dL_D);
    dL_D = -dL_D / double(n);
    vec dL_D_tri = dL_D(trimatl_ind(size(D)));
    for (uword i = 0; i < ntarget; i++)
    {
        gsl_vector_set(df, i, dL_D_tri(i));
    }
}

void ml_gsl_df_D_beta(const gsl_vector* v, void* pparams, gsl_vector *df)
{
    ML_Params* params = (ML_Params*)pparams;
    const mat* Xf = params->Xf;
    const vec* Yf = params->Yf;
    const mat* Zf = params->Zf;
    const size_t ngroup = params->ngroup;
    const uword n = params->n;
    const uword p = params->p;
    const uword q = params->q;
    size_t ntarget = p + q * (q + 1) / 2;
    vec D_tri(q * (q + 1) / 2, arma::fill::zeros), beta(p, arma::fill::zeros);
    for (size_t i = 0; i < p; i++)
    {
        beta(i) = gsl_vector_get(v, i);
    }
    for (size_t i = p; i < ntarget; i++)
    {
        D_tri(i - p) = gsl_vector_get(v, i);
    }
    mat D(q, q, arma::fill::zeros);
    D(trimatl_ind(size(D))) = D_tri;
    D = D.t();
    D(trimatl_ind(size(D))) = D_tri;
    mat dL_D;
    vec dL_beta;
    loglikelihood_d(Xf, Yf, Zf, ngroup, D, beta, n, dL_D, dL_beta);
    dL_D = -dL_D / double(n);
    dL_beta = -dL_beta / double(n);
    vec dL_D_tri = dL_D(trimatl_ind(size(D)));
    for (size_t i = 0; i < p; i++)
    {
        gsl_vector_set(df, i, dL_beta(i));
    }
    for (uword i = p; i < ntarget; i++)
    {
        gsl_vector_set(df, i, dL_D_tri(i - p));
    }
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

arma::mat covariance_from_cholesky(const gsl_vector* v, arma::uword offset, arma::uword q)
{
    arma::mat L(q, q, arma::fill::zeros);
    arma::uword k = offset;
    for (arma::uword i = 0; i < q; ++i)
    {
        for (arma::uword j = 0; j <= i; ++j, ++k)
        {
            double value = gsl_vector_get(v, k);
            L(i, j) = (i == j) ? std::exp(value) : value;
        }
    }
    return L * L.t();
}

void set_cholesky_start(gsl_vector* target, arma::uword offset, const arma::mat& D)
{
    arma::mat L;
    arma::mat D_sym = 0.5 * (D + D.t());
    if (!arma::chol(L, D_sym, "lower"))
    {
        L = arma::eye(D.n_rows, D.n_cols);
    }
    arma::uword k = offset;
    for (arma::uword i = 0; i < D.n_rows; ++i)
    {
        for (arma::uword j = 0; j <= i; ++j, ++k)
        {
            double value = (i == j) ? std::log(std::max(L(i, i), 1.0e-8)) : L(i, j);
            gsl_vector_set(target, k, value);
        }
    }
}

double ml_gsl_f_cholesky_D(const gsl_vector* v, void* p)
{
    ML_Params* params = static_cast<ML_Params*>(p);
    try
    {
        arma::mat D = covariance_from_cholesky(v, 0, params->q);
        double value = -loglikelihood(
            params->Xf, params->Yf, params->Zf, params->ngroup,
            D, *params->beta, params->n
        ) / double(params->n);
        return std::isfinite(value) ? value : DBL_MAX;
    }
    catch (...)
    {
        return DBL_MAX;
    }
}

double ml_gsl_f_cholesky_D_beta(const gsl_vector* v, void* p)
{
    ML_Params* params = static_cast<ML_Params*>(p);
    try
    {
        arma::vec beta(params->p);
        for (arma::uword i = 0; i < params->p; ++i)
        {
            beta(i) = gsl_vector_get(v, i);
        }
        arma::mat D = covariance_from_cholesky(v, params->p, params->q);
        double value = -loglikelihood(
            params->Xf, params->Yf, params->Zf, params->ngroup,
            D, beta, params->n
        ) / double(params->n);
        return std::isfinite(value) ? value : DBL_MAX;
    }
    catch (...)
    {
        return DBL_MAX;
    }
}

double HGWR::fit_D(ML_Params* params)
{
    const arma::uword q = D.n_cols;
    const arma::uword ntarget = q * (q + 1) / 2;
    gsl_multimin_function objective;
    objective.n = ntarget;
    objective.f = ml_gsl_f_cholesky_D;
    objective.params = params;

    const double start_scales[] = { 1.0, 0.1, 4.0 };
    double best_value = DBL_MAX;
    double best_measure = DBL_MAX;
    arma::mat best_D;
    bool found = false;
    int best_status = GSL_CONTINUE;

    for (size_t start_id = 0; start_id < 3; ++start_id)
    {
        gsl_vector* target = gsl_vector_alloc(ntarget);
        gsl_vector* step = gsl_vector_alloc(ntarget);
        arma::mat start_D = (start_id == 0)
            ? D
            : arma::eye(q, q) * start_scales[start_id];
        set_cholesky_start(target, 0, start_D);
        gsl_vector_set_all(step, std::max(alpha, DBL_EPSILON));

        gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(
            gsl_multimin_fminimizer_nmsimplex2, ntarget
        );
        int status = gsl_multimin_fminimizer_set(minimizer, &objective, target, step);
        size_t iter = 0;
        double measure = DBL_MAX;
        if (status == GSL_SUCCESS)
        {
            do
            {
                status = gsl_multimin_fminimizer_iterate(minimizer);
                if (status != GSL_SUCCESS) break;
                measure = gsl_multimin_fminimizer_size(minimizer);
                status = gsl_multimin_test_size(measure, eps_gradient);
                ++iter;
            }
            while (status == GSL_CONTINUE && iter < std::max<size_t>(1000, max_iters));
        }

        optimizer_iterations += iter;
        if (start_id > 0) ++optimizer_restarts;
        double value = minimizer->fval;
        if (status == GSL_SUCCESS && std::isfinite(value) && value < best_value)
        {
            best_value = value;
            best_measure = measure;
            best_D = covariance_from_cholesky(minimizer->x, 0, q);
            best_status = status;
            found = true;
        }

        gsl_multimin_fminimizer_free(minimizer);
        gsl_vector_free(step);
        gsl_vector_free(target);
    }

    optimizer_status = best_status;
    optimizer_measure = best_measure;
    if (!found)
    {
        ++optimizer_failures;
        return arma::datum::nan;
    }
    D = best_D;
    return best_value;
}

double HGWR::fit_D_beta(ML_Params* params)
{
    const arma::uword p = beta.n_rows;
    const arma::uword q = D.n_cols;
    const arma::uword ntarget = p + q * (q + 1) / 2;
    gsl_multimin_function objective;
    objective.n = ntarget;
    objective.f = ml_gsl_f_cholesky_D_beta;
    objective.params = params;

    const double start_scales[] = { 1.0, 0.1, 4.0 };
    double best_value = DBL_MAX;
    double best_measure = DBL_MAX;
    arma::mat best_D;
    arma::vec best_beta;
    bool found = false;
    int best_status = GSL_CONTINUE;

    for (size_t start_id = 0; start_id < 3; ++start_id)
    {
        gsl_vector* target = gsl_vector_alloc(ntarget);
        gsl_vector* step = gsl_vector_alloc(ntarget);
        for (arma::uword i = 0; i < p; ++i)
        {
            gsl_vector_set(target, i, beta(i));
        }
        arma::mat start_D = (start_id == 0)
            ? D
            : arma::eye(q, q) * start_scales[start_id];
        set_cholesky_start(target, p, start_D);
        gsl_vector_set_all(step, std::max(alpha, DBL_EPSILON));

        gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(
            gsl_multimin_fminimizer_nmsimplex2, ntarget
        );
        int status = gsl_multimin_fminimizer_set(minimizer, &objective, target, step);
        size_t iter = 0;
        double measure = DBL_MAX;
        if (status == GSL_SUCCESS)
        {
            do
            {
                status = gsl_multimin_fminimizer_iterate(minimizer);
                if (status != GSL_SUCCESS) break;
                measure = gsl_multimin_fminimizer_size(minimizer);
                status = gsl_multimin_test_size(measure, eps_gradient);
                ++iter;
            }
            while (status == GSL_CONTINUE && iter < std::max<size_t>(1000, max_iters));
        }

        optimizer_iterations += iter;
        if (start_id > 0) ++optimizer_restarts;
        double value = minimizer->fval;
        if (status == GSL_SUCCESS && std::isfinite(value) && value < best_value)
        {
            best_value = value;
            best_measure = measure;
            best_D = covariance_from_cholesky(minimizer->x, p, q);
            best_beta.set_size(p);
            for (arma::uword i = 0; i < p; ++i)
            {
                best_beta(i) = gsl_vector_get(minimizer->x, i);
            }
            best_status = status;
            found = true;
        }

        gsl_multimin_fminimizer_free(minimizer);
        gsl_vector_free(step);
        gsl_vector_free(target);
    }

    optimizer_status = best_status;
    optimizer_measure = best_measure;
    if (!found)
    {
        ++optimizer_failures;
        return arma::datum::nan;
    }
    D = best_D;
    beta = best_beta;
    return best_value;
}

void HGWR::fit_mu()
{
    mat D_inv = D.i();
    mu.fill(arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yhf[i];
        const mat& Zi = Zf[i];
        uword ndata = Zi.n_rows;
        mat Vi = Zi * D * Zi.t() + eye(ndata, ndata);
        mat Vi_inv = woodbury_eye(D_inv, Zi);
        vec Ri = Yi - Xi * beta;
        mu.row(i) = (D * Zi.t() * Vi_inv * Ri).t();
    }
}

double HGWR::fit_sigma()
{
    mat D_inv = D.i();
    double sigma2 = 0.0;
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yhf[i];
        const mat& Zi = Zf[i];
        uword ndata = Zi.n_rows;
        mat Vi = Zi * D * Zi.t() + eye(ndata, ndata);
        mat Vi_inv = woodbury_eye(D_inv, Zi);
        mat Ri = Yi - Xi * beta;
        sigma2 += as_scalar(Ri.t() * Vi_inv * Ri);
    }
    return sqrt(sigma2 / (double)ndata);
}

arma::mat make_spd(
    const arma::mat& A,
    double jitter = 1e-8,
    size_t* correction_count = nullptr
)
{
    arma::mat B = 0.5 * (A + A.t());

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, B);

    double min_eig = eigval.min();
    if (min_eig < jitter)
    {
        B += arma::eye(B.n_rows, B.n_cols) * (jitter - min_eig);
        if (correction_count != nullptr) ++(*correction_count);
    }

    return 0.5 * (B + B.t());
}

arma::mat safe_inv_sympd(
    const arma::mat& A,
    double jitter = 1e-8,
    size_t* correction_count = nullptr,
    size_t* jitter_count = nullptr,
    size_t* pseudoinverse_count = nullptr
)
{
    arma::mat B;
    arma::mat A_spd = make_spd(A, jitter, correction_count);

    bool ok = arma::inv_sympd(B, A_spd);
    if (ok) return B;

    for (int k = 0; k < 8; ++k)
    {
        double j = jitter * std::pow(10.0, k + 1);
        ok = arma::inv_sympd(B, A_spd + arma::eye(A.n_rows, A.n_cols) * j);
        if (ok)
        {
            if (jitter_count != nullptr) ++(*jitter_count);
            return B;
        }
    }

    if (pseudoinverse_count != nullptr) ++(*pseudoinverse_count);
    return arma::pinv(A_spd);
}

arma::mat safe_chol_lower(
    const arma::mat& A,
    double jitter = 1e-8,
    size_t* correction_count = nullptr,
    size_t* jitter_count = nullptr
)
{
    arma::mat L;
    arma::mat A_spd = make_spd(A, jitter, correction_count);

    bool ok = arma::chol(L, A_spd, "lower");
    if (ok) return L;

    for (int k = 0; k < 8; ++k)
    {
        double j = jitter * std::pow(10.0, k + 1);
        ok = arma::chol(L, A_spd + arma::eye(A.n_rows, A.n_cols) * j, "lower");
        if (ok)
        {
            if (jitter_count != nullptr) ++(*jitter_count);
            return L;
        }
    }

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A_spd);

    eigval.transform([jitter](double x)
    {
        return std::sqrt(std::max(x, jitter));
    });

    return eigvec * arma::diagmat(eigval);
}

int unit_intercept_column(const arma::mat& design, double tolerance = 1.0e-10)
{
    for (arma::uword column = 0; column < design.n_cols; ++column)
    {
        if (arma::max(arma::abs(design.col(column) - 1.0)) <= tolerance)
        {
            return static_cast<int>(column);
        }
    }
    return -1;
}

double relative_matrix_change(const arma::mat& current, const arma::mat& previous)
{
    return arma::norm(current - previous, "fro") /
        (1.0 + arma::norm(previous, "fro"));
}

double relative_scalar_change(double current, double previous)
{
    return std::abs(current - previous) / (1.0 + std::abs(previous));
}

void center_local_intercept(
    arma::mat& gamma,
    arma::vec& beta,
    const arma::vec& weights,
    int local_intercept,
    int fixed_intercept
)
{
    if (local_intercept < 0 || fixed_intercept < 0) return;
    double shift = arma::dot(weights, gamma.col(local_intercept)) / arma::sum(weights);
    gamma.col(local_intercept) -= shift;
    beta(fixed_intercept) += shift;
}

HGWR::Parameters HGWR::fit(const bool f_test)
{
    int precision = (int)log10(1 / eps_iter);
    double tss = sum((y - mean(y)) % (y - mean(y)));
    converged = false;
    stop_reason = "running";
    outer_iterations = 0;
    total_inner_updates = 0;
    final_inner_iterations = 0;
    inner_converged = false;
    all_inner_converged = false;
    inner_failures = 0;
    max_inner_iterations_used = 0;
    max_parameter_change = arma::datum::inf;
    optimizer_status = GSL_SUCCESS;
    optimizer_iterations = 0;
    optimizer_measure = arma::datum::nan;
    optimizer_failures = 0;
    optimizer_restarts = 0;
    bw_optimizer_status = GSL_SUCCESS;
    bw_optimizer_failures = 0;
    bw_objective = arma::datum::nan;
    bw_evaluations = 0;
    spd_corrections = 0;
    inverse_jitter_uses = 0;
    pseudoinverse_uses = 0;
    cholesky_jitter_uses = 0;
    gamma = mat(ngroup, nvg, arma::fill::zeros);
    gamma_se = mat(ngroup, nvg, arma::fill::zeros);
    beta = vec(nvx, arma::fill::zeros);
    mu = mat(ngroup, nvz, arma::fill::zeros);
    D = mat(nvz, nvz, arma::fill::eye);
    sigma = 1.0;
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
    group_weights = arma::conv_to<arma::vec>::from(group_size);
    uvec group_to = cumsum(group_size);
    uvec group_from = group_to - group_size;
    group_to = group_to - 1;
    transform(group_from.begin(), group_from.end(), group_to.begin(), group_span.begin(), [](uword from, uword to)
    {
        return span(from, to);
    });
    const int local_intercept = unit_intercept_column(G);
    const int fixed_intercept = unit_intercept_column(X);

    beta = fit_gls();
    fit_mu();

    size_t retry = 0;
    double rss = DBL_MAX;
    double rss_previous = DBL_MAX;
    double mlf = arma::datum::nan;
    for (size_t iter = 0; iter < max_iters; ++iter)
    {
        arma::mat gamma_previous = gamma;
        arma::vec beta_previous = beta;
        arma::mat mu_previous = mu;
        arma::mat D_previous = D;
        double sigma_previous = sigma;
        rss_previous = rss;

        for (uword i = 0; i < ngroup; i++)
        {
            Ygf[i] = Yf[i] - Xf[i] * beta;
        }
        fit_gwr();
        if (bw_optim && bw_optimizer_status != GSL_SUCCESS)
        {
            stop_reason = "bandwidth_optimizer_failed";
            outer_iterations = iter + 1;
            break;
        }
        center_local_intercept(
            gamma, beta, group_weights, local_intercept, fixed_intercept
        );

        vec hatMg = sum(G % gamma, 1);
        vec hatM = hatMg.rows(group);
        vec yh = y - hatM;
        for (uword i = 0; i < ngroup; i++)
        {
            Yhf[i] = Yf[i] - sum(G.row(i) % gamma.row(i));
        }
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
        if (!std::isfinite(mlf))
        {
            stop_reason = "ml_optimizer_failed";
            outer_iterations = iter + 1;
            break;
        }
        fit_mu();
        sigma = fit_sigma();

        vec yhat = yh - (X * beta) - sum(Z % (mu.rows(group)), 1);
        vec residual = yhat % yhat;
        rss = sum(residual);
        double rss_change = std::isfinite(rss_previous)
            ? relative_scalar_change(rss, rss_previous)
            : arma::datum::inf;
        max_parameter_change = std::max({
            relative_matrix_change(gamma, gamma_previous),
            relative_matrix_change(beta, beta_previous),
            relative_matrix_change(mu, mu_previous),
            relative_matrix_change(D, D_previous),
            relative_scalar_change(sigma, sigma_previous)
        });
        outer_iterations = iter + 1;

        if (rss < rss_previous)
        {
            if (retry > 0) retry = 0;
        }
        else if (iter > 0)
        {
            ++retry;
        }
        if (verbose > 0)
        {
            ostringstream sout;
            sout << fixed << setprecision(precision) << "Iter: " << iter;
            if (bw_optim) sout << ", " << "Bw: " << bw;
            sout << ", " << "RSS: " << rss;
            sout << ", " << "max change: " << max_parameter_change;
            sout << ", " << "R2: " << (1 - rss / tss);
            sout << ", " << "-loglik/n: " << mlf;
            if (retry > 0) sout << ", " << "Retry: " << retry;
            sout << endl;
            pcout(sout.str());
        }
        (*(this->pcancel))();

        const double parameter_tolerance = std::sqrt(eps_iter);
        if (rss_change <= eps_iter && max_parameter_change <= parameter_tolerance)
        {
            converged = true;
            stop_reason = "converged";
            break;
        }
        if (retry >= max_retries)
        {
            stop_reason = "retry_limit";
            break;
        }
    }
    if (!converged && stop_reason == "running")
    {
        stop_reason = "max_iters";
    }

    if (verbose > 0) pcout("Calculate GLSW diagnostics at the returned state\n");
    for (uword i = 0; i < ngroup; i++)
    {
        Ygf[i] = Yf[i] - Xf[i] * beta;
    }
    arma::mat gamma_returned = gamma;
    bool bw_optim_returned = bw_optim;
    bw_optim = false;
    fit_gwr(true, f_test);
    gamma = gamma_returned;
    bw_optim = bw_optim_returned;
    loglik = - mlf * double(ndata);
    calc_var_beta();
    return { gamma, beta, mu, D, sigma, bw };
}

/**
 * @brief Fit HGWR with back-fitting and an inner conditional MCMC estimator.
 *
 * Each outer iteration first updates GLSW effects, then runs a Gibbs-style
 * MCMC block for beta, group random effects, the relative random-effect
 * covariance D, and residual variance sigma^2. Conditional retained means
 * become the next back-fitting state; the outer loop is not a joint sampler.
 *
 * @param f_test If true, calculate GLSW F-test diagnostics after fitting.
 * @param niters Total MCMC iterations in each back-fitting iteration.
 * @param nburnin Number of initial MCMC iterations discarded as burn-in.
 * @return Estimated HGWR parameters.
 */
HGWR::Parameters HGWR::fit_mcmc_backfitting(const bool f_test)
{
    MonteCarloOptions options;
    return fit_mcmc_backfitting(f_test, options);
}

HGWR::Parameters HGWR::fit_mcmc_backfitting(
    const bool f_test,
    const MonteCarloOptions& mc_options
)
{
    const std::size_t niters = mc_options.iters;
    const std::size_t nburnin = mc_options.conditional_mode ? 0 : mc_options.burnin;
    if (niters == 0 || (!mc_options.conditional_mode && niters <= nburnin))
    {
        throw std::runtime_error("MCMC error: niters must be larger than nburnin.");
    }
    if (mc_options.beta_prior_variance <= 0.0 ||
        mc_options.d_prior_df_offset <= 1.0 ||
        mc_options.d_prior_mean <= 0.0 ||
        mc_options.sigma_prior_shape <= 0.0 ||
        mc_options.sigma_prior_scale <= 0.0 ||
        mc_options.inner_tolerance <= 0.0)
    {
        throw std::runtime_error("Monte Carlo prior parameters must be positive, with d_prior_df_offset > 1.");
    }

    int precision = int(std::log10(1.0 / eps_iter));
    double tss = arma::sum((y - arma::mean(y)) % (y - arma::mean(y)));
    converged = false;
    stop_reason = "running";
    outer_iterations = 0;
    total_inner_updates = 0;
    final_inner_iterations = 0;
    inner_converged = false;
    all_inner_converged = mc_options.conditional_mode;
    inner_failures = 0;
    max_inner_iterations_used = 0;
    max_parameter_change = arma::datum::inf;
    optimizer_status = GSL_SUCCESS;
    optimizer_iterations = 0;
    optimizer_measure = arma::datum::nan;
    optimizer_failures = 0;
    optimizer_restarts = 0;
    bw_optimizer_status = GSL_SUCCESS;
    bw_optimizer_failures = 0;
    bw_objective = arma::datum::nan;
    bw_evaluations = 0;
    spd_corrections = 0;
    inverse_jitter_uses = 0;
    pseudoinverse_uses = 0;
    cholesky_jitter_uses = 0;
    last_beta_draws.reset();
    last_D_draws.reset();
    last_sigma2_draws.reset();

    // -------------------------
    // 1. Initialise parameters
    // -------------------------
    gamma = arma::mat(ngroup, nvg, arma::fill::zeros);
    gamma_se = arma::mat(ngroup, nvg, arma::fill::zeros);

    beta = arma::vec(nvx, arma::fill::zeros);
    mu = arma::mat(ngroup, nvz, arma::fill::zeros);
    D = arma::mat(nvz, nvz, arma::fill::eye);
    sigma = 1.0;

    // -------------------------
    // 2. Prepare grouped data
    // -------------------------
    Zf = std::make_unique<arma::mat[]>(ngroup);
    Xf = std::make_unique<arma::mat[]>(ngroup);
    Yf = std::make_unique<arma::vec[]>(ngroup);
    Ygf = std::make_unique<arma::vec[]>(ngroup);
    Yhf = std::make_unique<arma::vec[]>(ngroup);

    arma::uvec group_size(ngroup);
    group_span.resize(ngroup);

    for (arma::uword j = 0; j < ngroup; ++j)
    {
        arma::uvec ind = arma::find(group == j);

        group_size(j) = ind.n_elem;

        Yf[j] = y.rows(ind);
        Xf[j] = X.rows(ind);
        Zf[j] = Z.rows(ind);

        Ygf[j] = y.rows(ind);
        Yhf[j] = y.rows(ind);
    }
    group_weights = arma::conv_to<arma::vec>::from(group_size);

    arma::uvec group_to = arma::cumsum(group_size);
    arma::uvec group_from = group_to - group_size;
    group_to = group_to - 1;

    std::transform(
        group_from.begin(),
        group_from.end(),
        group_to.begin(),
        group_span.begin(),
        [](arma::uword from, arma::uword to)
        {
            return arma::span(from, to);
        }
    );

    // -------------------------------------------------
    // 3. Initial beta using simple GLS under D = I
    // -------------------------------------------------
    beta = fit_gls();

    // Initial mu under the same parameterisation.
    fit_mu();

    const int local_intercept = unit_intercept_column(G);
    const int fixed_intercept = unit_intercept_column(X);

    // Initial sigma2 from residual.
    double sigma2_cur = 1.0;
    {
        double rss_init = 0.0;

        for (arma::uword j = 0; j < ngroup; ++j)
        {
            arma::vec rj = Yf[j] - Xf[j] * beta - Zf[j] * mu.row(j).t();
            rss_init += arma::dot(rj, rj);
        }

        if (std::isfinite(rss_init) && rss_init > 0.0)
        {
            sigma2_cur = rss_init / double(ndata);
            sigma = std::sqrt(sigma2_cur);
        }
    }

    // -------------------------
    // 4. MCMC hyperparameters
    // -------------------------

    // beta | sigma2 ~ N(beta0, sigma2 B0)
    const double tau_beta = mc_options.beta_prior_variance;
    arma::vec beta0(nvx, arma::fill::zeros);
    arma::mat B0_inv = arma::eye(nvx, nvx) / tau_beta;

    // D ~ Inv-Wishart(nu0, S0)
    // Use a weak but proper prior.
    const double nu0 = double(nvz) + mc_options.d_prior_df_offset;
    arma::mat D_prior_mean = arma::eye(nvz, nvz) * mc_options.d_prior_mean;
    arma::mat S0 = (nu0 - double(nvz) - 1.0) * D_prior_mean;
    S0 = make_spd(S0);

    // sigma2 ~ Inv-Gamma(a0, b0)
    const double a0 = mc_options.sigma_prior_shape;
    const double b0 = mc_options.sigma_prior_scale;

    // -------------------------
    // 5. Backfitting loop
    // -------------------------
    double rss = DBL_MAX;
    double rss_prev = DBL_MAX;
    double rel_diff = DBL_MAX;
    double mlf = 0.0;

    for (size_t bf_iter = 0; bf_iter < max_iters; bf_iter++)
    {
        arma::mat gamma_previous = gamma;
        arma::vec beta_previous = beta;
        arma::mat mu_previous = mu;
        arma::mat D_previous = D;
        double sigma_previous = sigma;
        rss_prev = rss;

        // =====================================================
        // Step A. Estimate gamma by original GWR-like estimator
        // =====================================================
        //
        // Yg = y - X beta
        //
        for (arma::uword j = 0; j < ngroup; ++j)
        {
            Ygf[j] = Yf[j] - Xf[j] * beta;
        }

        fit_gwr(false, false);
        if (bw_optim && bw_optimizer_status != GSL_SUCCESS)
        {
            stop_reason = "bandwidth_optimizer_failed";
            outer_iterations = bf_iter + 1;
            break;
        }
        center_local_intercept(
            gamma, beta, group_weights, local_intercept, fixed_intercept
        );

        // =====================================================
        // Step B. Construct Yh = y - G gamma
        // =====================================================
        //
        // For each group j:
        //   Yh_j = Y_j - G_j gamma_j
        //
        for (arma::uword j = 0; j < ngroup; ++j)
        {
            double gj = arma::as_scalar(G.row(j) * gamma.row(j).t());
            Yhf[j] = Yf[j] - gj;
        }

        // =====================================================
        // Step C. MCMC block for beta, D, sigma2, mu
        // =====================================================
        //
        // Conditional model:
        //   Yh_j = X_j beta + Z_j mu_j + e_j
        //   mu_j ~ N(0, sigma2 D)
        //   e_j  ~ N(0, sigma2 I)
        //
        arma::vec beta_cur = beta;
        arma::mat D_cur = make_spd(D, 1.0e-8, &spd_corrections);
        arma::mat mu_cur = mu;

        // Use previous sigma as initial value.
        sigma2_cur = sigma * sigma;
        if (!std::isfinite(sigma2_cur) || sigma2_cur <= 0.0)
        {
            sigma2_cur = 1.0;
        }

        size_t nkeep = 0;

        arma::vec beta_sum(nvx, arma::fill::zeros);
        arma::mat D_sum(nvz, nvz, arma::fill::zeros);
        arma::mat mu_sum(ngroup, nvz, arma::fill::zeros);
        double sigma2_sum = 0.0;
        // Increase the retained Monte Carlo budget during early backfitting
        // iterations, but cap it to keep a non-stopping fit computationally
        // bounded. Reaching the outer tolerance remains an operational stop,
        // not a claim of stochastic fixed-point convergence.
        const size_t budget_increments = std::min<std::size_t>(bf_iter, 4);
        const size_t outer_niters = mc_options.conditional_mode
            ? niters
            : niters + budget_increments * (niters - nburnin);
        arma::mat beta_draws;
        arma::mat D_draws;
        arma::vec sigma2_draws;
        if (mc_options.save_draws && !mc_options.conditional_mode)
        {
            const size_t expected = outer_niters - nburnin;
            beta_draws.set_size(expected, nvx);
            D_draws.set_size(expected, nvz * nvz);
            sigma2_draws.set_size(expected);
        }

        bool current_inner_converged = false;
        size_t current_inner_iterations = 0;

        for (size_t mcmc_iter = 0; mcmc_iter < outer_niters; ++mcmc_iter)
        {
            arma::vec beta_inner_previous = beta_cur;
            arma::mat mu_inner_previous = mu_cur;
            arma::mat D_inner_previous = D_cur;
            double sigma2_inner_previous = sigma2_cur;
            // -------------------------------------------------
            // C1. Sample mu_j | beta, D, sigma2, Yh
            // -------------------------------------------------
            arma::mat D_inv = safe_inv_sympd(
                D_cur, 1.0e-8, &spd_corrections, &inverse_jitter_uses,
                &pseudoinverse_uses
            );
            arma::mat mu_new(ngroup, nvz, arma::fill::zeros);

            for (arma::uword j = 0; j < ngroup; ++j)
            {
                const arma::mat& Xj = Xf[j];
                const arma::mat& Zj = Zf[j];
                const arma::vec& Yhj = Yhf[j];

                arma::mat C_mu = safe_inv_sympd(
                    Zj.t() * Zj + D_inv, 1.0e-8, &spd_corrections,
                    &inverse_jitter_uses, &pseudoinverse_uses
                );
                arma::vec m_mu = C_mu * Zj.t() * (Yhj - Xj * beta_cur);

                arma::vec draw_mu = m_mu;
                if (!mc_options.conditional_mode)
                {
                    arma::mat L_mu = safe_chol_lower(
                        sigma2_cur * C_mu, 1.0e-8, &spd_corrections,
                        &cholesky_jitter_uses
                    );
                    draw_mu += L_mu * arma::randn(nvz);
                }

                mu_new.row(j) = draw_mu.t();
            }

            mu_cur = mu_new;

            // -------------------------------------------------
            // C2. Sample beta | mu, sigma2, Yh
            // -------------------------------------------------
            arma::mat XtX(nvx, nvx, arma::fill::zeros);
            arma::vec XtY(nvx, arma::fill::zeros);

            for (arma::uword j = 0; j < ngroup; ++j)
            {
                const arma::mat& Xj = Xf[j];
                const arma::mat& Zj = Zf[j];
                const arma::vec& Yhj = Yhf[j];

                arma::vec muj = mu_cur.row(j).t();
                arma::vec y_tilde = Yhj - Zj * muj;

                XtX += Xj.t() * Xj;
                XtY += Xj.t() * y_tilde;
            }

            arma::mat C_beta = safe_inv_sympd(
                XtX + B0_inv, 1.0e-8, &spd_corrections,
                &inverse_jitter_uses, &pseudoinverse_uses
            );
            arma::vec m_beta = C_beta * (XtY + B0_inv * beta0);

            beta_cur = m_beta;
            if (!mc_options.conditional_mode)
            {
                arma::mat L_beta = safe_chol_lower(
                    sigma2_cur * C_beta, 1.0e-8, &spd_corrections,
                    &cholesky_jitter_uses
                );
                beta_cur += L_beta * arma::randn(nvx);
            }

            // -------------------------------------------------
            // C3. Sample D | mu, sigma2
            // -------------------------------------------------
            //
            // Since mu_j ~ N(0, sigma2 D),
            // the sufficient statistic is:
            //   sum(mu_j mu_j') / sigma2
            //
            arma::mat SS(nvz, nvz, arma::fill::zeros);

            for (arma::uword j = 0; j < ngroup; ++j)
            {
                arma::vec muj = mu_cur.row(j).t();
                SS += muj * muj.t();
            }

            arma::mat S_post = make_spd(
                S0 + SS / sigma2_cur, 1.0e-8, &spd_corrections
            );
            double nu_post = nu0 + double(ngroup);

            if (mc_options.conditional_mode)
            {
                D_cur = S_post / (nu_post + double(nvz) + 1.0);
            }
            else
            {
                D_cur = rinvwishart(nu_post, S_post);
            }
            D_cur = make_spd(D_cur, 1.0e-8, &spd_corrections);

            // -------------------------------------------------
            // C4. Sample sigma2 | beta, mu, D, Yh
            // -------------------------------------------------
            D_inv = safe_inv_sympd(
                D_cur, 1.0e-8, &spd_corrections, &inverse_jitter_uses,
                &pseudoinverse_uses
            );

            double rss_mcmc = 0.0;
            double re_quad = 0.0;

            for (arma::uword j = 0; j < ngroup; ++j)
            {
                const arma::mat& Xj = Xf[j];
                const arma::mat& Zj = Zf[j];
                const arma::vec& Yhj = Yhf[j];

                arma::vec muj = mu_cur.row(j).t();
                arma::vec resid = Yhj - Xj * beta_cur - Zj * muj;

                rss_mcmc += arma::dot(resid, resid);
                re_quad += arma::as_scalar(muj.t() * D_inv * muj);
            }

            arma::vec beta_diff = beta_cur - beta0;
            double beta_quad = arma::as_scalar(beta_diff.t() * B0_inv * beta_diff);

            double a_post = a0 + 0.5 * double(ndata + ngroup * nvz + nvx);
            double b_post = b0 + 0.5 * (rss_mcmc + re_quad + beta_quad);

            sigma2_cur = mc_options.conditional_mode
                ? b_post / (a_post + 1.0)
                : rinvgamma(a_post, b_post);

            if (!std::isfinite(sigma2_cur) || sigma2_cur <= 0.0)
            {
                sigma2_cur = std::max(1.0e-8, rss_mcmc / double(ndata));
            }

            // -------------------------------------------------
            // C5. Store draws from the current conditional hierarchical chain
            // -------------------------------------------------
            ++current_inner_iterations;
            ++total_inner_updates;

            double inner_change = std::max({
                relative_matrix_change(beta_cur, beta_inner_previous),
                relative_matrix_change(mu_cur, mu_inner_previous),
                relative_matrix_change(D_cur, D_inner_previous),
                relative_scalar_change(sigma2_cur, sigma2_inner_previous)
            });
            if (mc_options.conditional_mode && mcmc_iter >= 4 &&
                inner_change <= mc_options.inner_tolerance)
            {
                current_inner_converged = true;
            }

            if ((!mc_options.conditional_mode && mcmc_iter >= nburnin) ||
                (mc_options.conditional_mode &&
                 (current_inner_converged || mcmc_iter + 1 == outer_niters)))
            {
                beta_sum += beta_cur;
                D_sum += D_cur;
                mu_sum += mu_cur;
                sigma2_sum += sigma2_cur;
                if (mc_options.save_draws && !mc_options.conditional_mode)
                {
                    beta_draws.row(nkeep) = beta_cur.t();
                    D_draws.row(nkeep) = arma::vectorise(D_cur).t();
                    sigma2_draws(nkeep) = sigma2_cur;
                }
                ++nkeep;
            }

            if (current_inner_converged) break;

            if (verbose > 1)
            {
                std::ostringstream sout;
                sout << "BF iter=" << bf_iter
                     << ", MCMC iter=" << mcmc_iter
                     << ", rss_mcmc=" << rss_mcmc
                     << ", sigma2=" << sigma2_cur
                     << "\n";
                pcout(sout.str());
            }

            (*(pcancel))();
        }

        if (nkeep == 0)
        {
            throw std::runtime_error("MCMC error: no conditional samples retained.");
        }

        // =====================================================
        // Step D. Conditional retained means become new Backfitting values
        // =====================================================
        beta = beta_sum / double(nkeep);
        D = make_spd(
            D_sum / double(nkeep), 1.0e-8, &spd_corrections
        );
        mu = mu_sum / double(nkeep);
        sigma = std::sqrt(sigma2_sum / double(nkeep));
        final_inner_iterations = current_inner_iterations;
        inner_converged = current_inner_converged;
        max_inner_iterations_used = std::max(
            max_inner_iterations_used, current_inner_iterations
        );
        if (mc_options.conditional_mode && !current_inner_converged)
        {
            all_inner_converged = false;
            ++inner_failures;
        }
        if (mc_options.save_draws && !mc_options.conditional_mode)
        {
            last_beta_draws = beta_draws;
            last_D_draws = D_draws;
            last_sigma2_draws = sigma2_draws;
        }

        // =====================================================
        // Step E. Calculate full HGWR residual
        // =====================================================
        arma::vec fitted_glsw_group = arma::sum(G % gamma, 1);
        arma::vec fitted_glsw_sample = fitted_glsw_group.rows(group);

        arma::vec fitted_fixed = X * beta;
        arma::vec fitted_random = arma::sum(Z % mu.rows(group), 1);

        arma::vec resid_full = y - fitted_glsw_sample - fitted_fixed - fitted_random;

        rss = arma::dot(resid_full, resid_full);

        if (bf_iter == 0 || !std::isfinite(rss_prev))
        {
            rel_diff = DBL_MAX;
        }
        else
        {
            rel_diff = std::abs(rss - rss_prev) / (rss_prev + 1.0e-12);
        }
        max_parameter_change = std::max({
            relative_matrix_change(gamma, gamma_previous),
            relative_matrix_change(beta, beta_previous),
            relative_matrix_change(mu, mu_previous),
            relative_matrix_change(D, D_previous),
            relative_scalar_change(sigma, sigma_previous)
        });
        outer_iterations = bf_iter + 1;

        // Existing marginal/profile likelihood diagnostic.
        mlf = -loglikelihood(Xf.get(), Yhf.get(), Zf.get(), ngroup, D, beta, ndata) / double(ndata);

        if (verbose > 0)
        {
            std::ostringstream sout;
            sout << std::fixed << std::setprecision(precision)
                 << "BF-MCMC Iter: " << bf_iter;

            if (bw_optim)
            {
                sout << ", Bw: " << bw;
            }

            sout << ", RSS: " << rss
                 << ", rel_dRSS: " << rel_diff
                 << ", max change: " << max_parameter_change
                 << ", R2: " << (1.0 - rss / tss)
                 << ", sigma: " << sigma
                 << ", -loglik/n: " << mlf
                 << ", beta: " << beta.t()
                 << "\n";

            pcout(sout.str());
        }

        (*(this->pcancel))();

        const double parameter_tolerance = std::sqrt(eps_iter);
        if (rel_diff <= eps_iter && max_parameter_change <= parameter_tolerance)
        {
            converged = true;
            stop_reason = "tolerance_reached";
            break;
        }
    }

    if (!converged && stop_reason == "running")
    {
        stop_reason = "max_iters";
    }

    // =========================================================
    // 6. Final re-fit gamma for t-test / f-test diagnostics
    // =========================================================
    if (verbose > 0) pcout("Calculate GLSW diagnostics at the returned state\n");

    for (arma::uword j = 0; j < ngroup; ++j)
    {
        Ygf[j] = Yf[j] - Xf[j] * beta;
    }

    arma::mat gamma_returned = gamma;
    bool bw_optim_returned = bw_optim;
    bw_optim = false;
    fit_gwr(true, f_test);
    gamma = gamma_returned;
    bw_optim = bw_optim_returned;

    // =========================================================
    // 7. Diagnostics
    // =========================================================
    loglik = -mlf * double(ndata);

    // calc_var_beta() is still based on the marginal GLS formula.
    // It is an approximate point-estimate diagnostic. The saved conditional
    // draws do not provide joint HGWR uncertainty for beta or gamma.
    calc_var_beta();

    return { gamma, beta, mu, D, sigma, bw };
}

void HGWR::calc_var_beta()
{
    mat D_inv = D.i(), XtViX(X.n_cols, X.n_cols, arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Zi = Zf[i];
        uword ndata = Zi.n_rows;
        mat Vi = Zi * D * Zi.t() + eye(ndata, ndata);
        mat Vi_inv = woodbury_eye(D_inv, Zi);
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
    mat D_inv = D.i();
    unique_ptr<mat[]> Vf = make_unique<mat[]>(ngroup);
    unique_ptr<mat[]> GVGf = make_unique<mat[]>(ngroup);
    unique_ptr<mat[]> GVf = make_unique<mat[]>(ngroup);
    for (size_t i = 0; i < ngroup; i++)
    {
        uvec ind = find(group == i);
        const mat& Zi = Zf[i];
        Vf[i] = Zi * D * Zi.t() + eye(Zi.n_rows, Zi.n_rows);
        mat Vi_inv = woodbury_eye(D_inv, Zi);
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
        double pv = arma::datum::nan;
        if (std::isfinite(fv) && fv >= 0.0 &&
            std::isfinite(df1) && df1 > 0.0 &&
            std::isfinite(df2) && df2 > 0.0)
        {
            gsl_error_handler_t* previous_handler = gsl_set_error_handler_off();
            pv = gsl_cdf_fdist_Q(fv, df1, df2);
            gsl_set_error_handler(previous_handler);
        }
        vec4 result = { fv, df1, df2, pv };
        results.push_back(result);
    }
    return results;
}
