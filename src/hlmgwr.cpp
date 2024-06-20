#ifdef HGWRR_RCPP
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#include "hlmgwr.h"
#include <sstream>
#include <iomanip>
#include <string>
#include <utility>
#include <armadillo>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace arma;

const double log2pi = log(2.0 * M_PI);

struct ML_Params
{
    const mat* Xf;
    const vec* Yf;
    const mat* Zf;
    const vec* beta;
    size_t ngroup;
    uword n;
    uword p;
    uword q;
};

inline vec gwr_kernel_gaussian2(vec dist2, double bw2)
{
    return exp(- dist2 / (2.0 * bw2));
}

inline vec gwr_kernel_bisquare2(vec dist2, double bw2)
{
    return ((1 - dist2 / bw2) % (1 - dist2 / bw2)) % (dist2 < bw2);
}

const GWRKernelFunctionSquared gwr_kernel_functions[] = {
    gwr_kernel_gaussian2,
    gwr_kernel_bisquare2
};

double criterion_bw(const double bw, const HLMGWRBWArgs* args, const bool verbose, const PrintFunction pcout)
{
    const mat G = args->G;
    const mat Vig = args->Vig;
    const vec Viy = args->Viy;
    const mat u = args->u;
    const size_t ngroup = Viy.n_rows;
    /// Calibrate for each gorup.
    GWRKernelFunctionSquared gwr_kernel = args->gwr_kernel;
    double cv = 0.0;
    for (size_t i = 0; i < ngroup; i++)
    {
        mat d_u = u.each_row() - u.row(i);
        vec d2 = sum(d_u % d_u, 1);
        double b2 = vec(sort(d2))[(int)bw - 1];
        vec wW = (*gwr_kernel)(d2, b2);
        wW(i) = 0;
        mat GtWVG = (G.each_col() % wW).t() * Vig;
        mat GtWVy = (G.each_col() % wW).t() * Viy;
        try
        {
            vec bi = solve(GtWVG, GtWVy);
            double yhi = as_scalar(Vig.row(i) * bi);
            double residual = Viy(i) - yhi;
            cv += residual * residual;
        }
        catch(const std::exception& e)
        {
            return DBL_MAX;
        }
    }
    if (verbose)
    {
        ostringstream sout;
        sout << "bw: " << bw << "; " << "cv: " << cv << "\r";
        pcout(sout.str());
    }
    return cv;
}

double golden_selection(const double lower, const double upper, const bool adaptive, const HLMGWRBWArgs* args, const bool verbose, const PrintFunction pcout)
{
    double xU = upper, xL = lower;
    bool adaptBw = adaptive;
    const double eps = 1e-4;
    const double R = (sqrt(5)-1)/2;
    int iter = 0;
    double d = R * (xU - xL);
    double x1 = adaptBw ? floor(xL + d) : (xL + d);
    double x2 = adaptBw ? round(xU - d) : (xU - d);
    double f1 = criterion_bw(x1, args, verbose, pcout);
    double f2 = criterion_bw(x2, args, verbose, pcout);
    double d1 = f2 - f1;
    double xopt = f1 < f2 ? x1 : x2;
    double ea = 100;
    while ((fabs(d) > eps) && (fabs(d1) > eps) && iter < ea)
    {
        d = R * d;
        if (f1 < f2)
        {
            xL = x2;
            x2 = x1;
            x1 = adaptBw ? round(xL + d) : (xL + d);
            f2 = f1;
            f1 = criterion_bw(x1, args, verbose, pcout);
        }
        else
        {
            xU = x1;
            x1 = x2;
            x2 = adaptBw ? floor(xU - d) : (xU - d);
            f1 = f2;
            f2 = criterion_bw(x2, args, verbose, pcout);
        }
        iter = iter + 1;
        xopt = (f1 < f2) ? x1 : x2;
        d1 = f2 - f1;
    }
    if (verbose)
    {
        pcout("\n");
    }
    return xopt;
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
tuple<mat, double> fit_gwr(
    const mat& G, const vec* Yf, const mat* Zf, const size_t ngroup,
    const mat& D, const mat& u, const double bw, const GWRKernelType kernel,
    const size_t verbose, const PrintFunction pcout
)
{
    uword k = G.n_cols;//, q = Zf[0].n_cols;
    mat beta(ngroup, k, arma::fill::zeros), D_inv = D.i();
    mat Vig(ngroup, k, arma::fill::zeros);
    vec Viy(ngroup, arma::fill::zeros);
    for (size_t i = 0; i < ngroup; i++)
    {
        const mat& Yi = Yf[i];
        const mat& Zi = Zf[i];
        uword ndata = Zi.n_rows;
        mat Vi_inv = eye(ndata, ndata) - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
        mat Visigma = ones(1, ndata) * Vi_inv;
        Vig.row(i) = Visigma * ones(ndata, 1) * G.row(i);
        Viy(i) = as_scalar(Visigma * Yi);
    }
    GWRKernelFunctionSquared gwr_kernel = gwr_kernel_functions[kernel];
    /// Check whether need to optimize bw
    double b = bw;
    if (bw == 0.0)
    {
        HLMGWRBWArgs bw_args { G, Vig, Viy, u, gwr_kernel };
        double upper = ngroup, lower = k + 1;
        b = golden_selection(lower, upper, true, &bw_args, verbose > 1, pcout);
    }
    /// Calibrate for each gorup.
    for (size_t i = 0; i < ngroup; i++)
    {
        mat d_u = u.each_row() - u.row(i);
        vec d2 = sum(d_u % d_u, 1);
        double b2 = vec(sort(d2))[(int)b - 1];
        vec wW = (*gwr_kernel)(d2, b2);
        mat GtWVG = (G.each_col() % wW).t() * Vig;
        mat GtWVy = (G.each_col() % wW).t() * Viy;
        beta.row(i) = solve(GtWVG, GtWVy).t();
    }
    return make_tuple(beta, b);
}

vec fit_gls(const mat* Xf, const vec* Yf, const mat* Zf, const size_t ngroup, const mat& D)
{
    uword p = Xf[0].n_cols;
    mat XtWX(p, p, arma::fill::zeros);
    vec XtWY(p, arma::fill::zeros);
    mat D_inv = D.i();
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yf[i];
        const mat& Zi = Zf[i];
        uword ndata = Zi.n_rows;
        mat Vi_inv = eye(ndata, ndata) - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
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
        uword nidata = Zi.n_rows;
        mat Ii = eye<mat>(nidata, nidata);
        mat Vi = ((Zi * D) * Zi.t()) + Ii;
        double detVi, sign_detVi;
        log_det(detVi, sign_detVi, Vi);
        mat Vi_inv = Ii - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
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
        uword nidata = Zi.n_rows;
        mat Vi_inv = eye(nidata, nidata) - Zi * inv(D_inv + Zi.t() * Zi) * Zi.t();
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
        uword nidata = Zi.n_rows;
        mat Vi_inv = eye(nidata, nidata) - Zi * inv(D_inv + Zi.t() * Zi) * Zi.t();
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
        gsl_vector_set(df, i, dL_D(i));
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
    *f = ml_gsl_f_D(v, p);
    ml_gsl_df_D_beta(v, p, df);
}

double fit_D(mat& D, const ML_Params* params, const double alpha, const double eps, const size_t max_iters, const bool verbose, const PrintFunction pcout)
{
    int precision = int(log10(1.0 / eps));
    uword q = D.n_cols, ntarget = q * (q + 1) / 2;
    gsl_multimin_function_fdf minex_fun;
    minex_fun.n = ntarget;
    minex_fun.f = ml_gsl_f_D;
    minex_fun.df = ml_gsl_df_D;
    minex_fun.fdf = ml_gsl_fdf_D;
    minex_fun.params = (void*)params;
    gsl_vector *target = gsl_vector_alloc(ntarget), *step_size = gsl_vector_alloc(ntarget);
    uvec D_tril_idx = trimatl_ind(arma::size(D)), D_triu_idx = trimatu_ind(arma::size(D));
    vec D_tril_vec = D(D_tril_idx);
    for (uword i = 0; i < ntarget; i++)
    {
        gsl_vector_set(target, i, D_tril_vec(i));
        gsl_vector_set(step_size, i, alpha);
    }
    gsl_vector *x0 = gsl_vector_alloc(ntarget);
    gsl_vector_memcpy(x0, target);
    gsl_multimin_fdfminimizer *minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, ntarget);
    gsl_multimin_fdfminimizer_set(minimizer, &minex_fun, target, alpha, eps);
    if (verbose)
    {
        ostringstream sout;
        sout << setprecision(precision) << fixed << minimizer->x->data[0];
        for (size_t i = 1; i < D_tril_vec.n_elem; i++)
        {
            sout << "," << minimizer->x->data[i];
        }
        sout << ";";
        sout << setprecision(precision) << fixed << minimizer->gradient->data[0] << ",";
        for (size_t i = 1; i < D_tril_vec.n_elem; i++)
        {
            sout << "," << minimizer->gradient->data[i];
        }
        sout << ";";
        sout << minimizer->f << '\r';
        pcout(sout.str());
    }
    size_t iter = 0;
    int status;
    do
    {
        gsl_vector_memcpy(x0, minimizer->x);
        status = gsl_multimin_fdfminimizer_iterate(minimizer);
        if (verbose)
        {
            ostringstream sout;
            sout << setprecision(precision) << fixed << minimizer->x->data[0];
            for (size_t i = 1; i < D_tril_vec.n_elem; i++)
            {
                sout << "," << minimizer->x->data[i];
            }
            sout << ";";
            sout << setprecision(precision) << fixed << minimizer->gradient->data[0];
            for (size_t i = 1; i < D_tril_vec.n_elem; i++)
            {
                sout << "," << minimizer->gradient->data[i];
            }
            sout << ";";
            sout << minimizer->f << '\r';
            pcout(sout.str());
        }
        if (status || gsl_isnan(minimizer->f)) break;
        status = gsl_multimin_test_gradient(minimizer->gradient, eps);
    } while (status == GSL_CONTINUE && (++iter) < max_iters);
    if (verbose) 
    {
        pcout("\n");
    }
    if (!gsl_isnan(minimizer->f))
    {
        mat D1 = mat(arma::size(D), arma::fill::eye);
        vec D_tri(arma::size(D_tril_idx));
        for (uword i = 0; i < ntarget; i++)
        {
            D_tri(i) = gsl_vector_get(minimizer->x, i);
        }
        D1(D_tril_idx) = D_tri;
        D1 = D1.t();
        D1(D_tril_idx) = D_tri;
        D = D1;
    }
    return minimizer->f;
}

void fit_D_beta(mat& D, vec& beta, const ML_Params* params, const double alpha, const double eps, const size_t max_iters, const bool verbose, const PrintFunction pcout)
{
    int precision = int(log10(1.0 / eps));
    uword p = beta.n_rows, q = D.n_cols, ntarget = p + q * (q + 1) / 2;
    gsl_multimin_function_fdf minex_fun;
    minex_fun.n = ntarget;
    minex_fun.f = ml_gsl_f_D_beta;
    minex_fun.df = ml_gsl_df_D_beta;
    minex_fun.fdf = ml_gsl_fdf_D_beta;
    minex_fun.params = (void*)params;
    gsl_vector *target = gsl_vector_alloc(ntarget);//, *step_size = gsl_vector_alloc(ntarget);
    for (uword i = 0; i < p; i++)
    {
        gsl_vector_set(target, i, beta(i));
    }
    uvec D_tril_idx = trimatl_ind(arma::size(D)), D_triu_idx = trimatu_ind(arma::size(D));
    vec D_tril_vec = D(D_tril_idx);
    for (uword i = p; i < ntarget; i++)
    {
        uword e = i - p;
        gsl_vector_set(target, i, D_tril_vec(e));
    }
    gsl_vector *x0 = gsl_vector_alloc(ntarget);
    gsl_vector_memcpy(x0, target);
    gsl_multimin_fdfminimizer *minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, ntarget);
    gsl_multimin_fdfminimizer_set(minimizer, &minex_fun, target, alpha, eps);
    if (verbose)
    {
        ostringstream sout;
        sout << setprecision(precision) << fixed;
        for (size_t i = 0; i < p; i++)
        {
            sout << minimizer->x->data[i] << ",";
        }
        sout << minimizer->x->data[p];
        for (size_t i = 1; i < D_tril_vec.n_elem; i++)
        {
            sout << "," << minimizer->x->data[p + i];
        }
        sout << ";";
        sout << setprecision(precision) << fixed;
        for (size_t i = 0; i < p; i++)
        {
            sout << minimizer->gradient->data[i] << ",";
        }
        sout << minimizer->x->data[p];
        for (size_t i = 1; i < D_tril_vec.n_elem; i++)
        {
            sout << "," << minimizer->gradient->data[p + i];
        }
        sout << ";"; 
        sout << minimizer->f << '\r';
        pcout(sout.str());
    }
    size_t iter = 0;
    int status;
    do
    {
        gsl_vector_memcpy(x0, minimizer->x);
        status = gsl_multimin_fdfminimizer_iterate(minimizer);
        if (verbose)
        {
            ostringstream sout;
            sout << setprecision(precision) << fixed;
            for (size_t i = 0; i < p; i++)
            {
                sout << minimizer->x->data[i] << ",";
            }
            sout << minimizer->x->data[p];
            for (size_t i = 1; i < D_tril_vec.n_elem; i++)
            {
                sout << "," << minimizer->x->data[p + i];
            }
            sout << ";";
            sout << setprecision(precision) << fixed;
            for (size_t i = 0; i < p; i++)
            {
                sout << minimizer->gradient->data[i] << ",";
            }
            sout << minimizer->x->data[p];
            for (size_t i = 1; i < D_tril_vec.n_elem; i++)
            {
                sout << "," << minimizer->gradient->data[p + i];
            }
            sout << ";"; 
            sout << minimizer->f << '\r';
            pcout(sout.str());
        }
        if (status || gsl_isnan(minimizer->f)) break;
        status = gsl_multimin_test_gradient(minimizer->gradient, eps);
    } while (status == GSL_CONTINUE && (++iter) < max_iters);
    if (verbose) 
    {
        pcout("\n");
    }
    mat D1(arma::size(D), arma::fill::eye);
    vec beta1(arma::size(beta), arma::fill::ones);
    if (!gsl_isnan(minimizer->f))
    {
        vec D_tri(arma::size(D_tril_idx));
        for (uword i = p; i < ntarget; i++)
        {
            D_tri(i - p) = gsl_vector_get(minimizer->x, i);
        }
        D1(D_tril_idx) = D_tri;
        D1 = D1.t();
        D1(D_tril_idx) = D_tri;
        for (uword i = 0; i < p; i++)
        {
            beta1(i) = gsl_vector_get(minimizer->x, i);
        }
    }
    D = D1;
    beta = beta1;
}

mat fit_mu(const mat* Xf, const vec* Yf, const mat* Zf, const size_t ngroup, const vec& beta, const mat& D)
{
    uword q = Zf[0].n_cols;
    mat D_inv = D.i(), mu(ngroup, q, arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yf[i];
        const mat& Zi = Zf[i];
        uword ndata = Zi.n_rows;
        mat Vi = Zi * D * Zi.t() + eye(ndata, ndata);
        mat Vi_inv = eye(ndata, ndata) - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
        vec Ri = Yi - Xi * beta;
        mu.row(i) = (D * Zi.t() * Vi_inv * Ri).t();
    }
    return mu;
}

double fit_sigma(const mat* Xf, const vec* Yf, const mat* Zf, const size_t ngroup, const size_t ndata, const vec& beta, const mat& D)
{
    mat D_inv = D.i();
    double sigma2 = 0.0;
    for (uword i = 0; i < ngroup; i++)
    {
        const mat& Xi = Xf[i];
        const mat& Yi = Yf[i];
        const mat& Zi = Zf[i];
        uword ndata = Zi.n_rows;
        mat Vi = Zi * D * Zi.t() + eye(ndata, ndata);
        mat Vi_inv = eye(ndata, ndata) - Zi * (D_inv + Zi.t() * Zi).i() * Zi.t();
        mat Ri = Yi - Xi * beta;
        sigma2 += as_scalar(Ri.t() * Vi_inv * Ri);
    }
    return sqrt(sigma2 / (double)ndata);
}

HLMGWRParams backfitting_maximum_likelihood(const HLMGWRArgs& args, const HLMGWROptions& options, const PrintFunction pcout)
{
    //================
    // Extract Options
    //================
    double alpha = options.alpha;
    double eps_iter = options.eps_iter;
    double eps_gradient = options.eps_gradient;
    size_t max_iters = options.max_iters;
    size_t max_retries = options.max_retries;
    size_t verbose = options.verbose;
    size_t ml_type = options.ml_type;
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
    double bw = args.bw, bw_optim = bw;
    GWRKernelType kernel = args.kernel;
    uword ngroup = G.n_rows, ndata = X.n_rows;
    uword nvg = G.n_cols, nvx = X.n_cols, nvz = Z.n_cols;
    double tss = sum((y - mean(y)) % (y - mean(y)));
    mat gamma(ngroup, nvg, arma::fill::zeros);
    vec beta(nvx, arma::fill::zeros);
    mat mu(ngroup, nvz, arma::fill::zeros);
    mat D(nvz, nvz, arma::fill::eye);
    mat gmap(ndata, ngroup, arma::fill::zeros);
    for (uword i = 0; i < ngroup; i++)
    {
        gmap.col(i) = conv_to<vec>::from(group == i);
    }
    mat* Zf = new mat[ngroup];
    mat* Xf = new mat[ngroup];
    vec* Yf = new vec[ngroup];
    vec* Ygf = new vec[ngroup];
    vec* Yhf = new vec[ngroup];
    for (uword i = 0; i < ngroup; i++)
    {
        uvec ind = find(group == i);
        Yf[i] = y.rows(ind);
        Xf[i] = X.rows(ind);
        Zf[i] = Z.rows(ind);
    }
    //----------------------------------------------
    // Generalized Least Squared Estimation for beta
    //----------------------------------------------
    beta = fit_gls(Xf, Yf, Zf, ngroup, D);
    //============
    // Backfitting
    //============
    size_t retry = 0;
    double rss = DBL_MAX, rss0 = DBL_MAX, diff = DBL_MAX, mlf = 0.0;
    for (size_t iter = 0; (abs(diff) > eps_iter) && iter < max_iters && retry < max_retries; iter++)
    {
        rss0 = rss;
        //--------------------
        // Initial Guess for M
        //--------------------
        for (uword i = 0; i < ngroup; i++)
        {
            Ygf[i] = Yf[i] - Xf[i] * beta;
        }
        auto gwr_result = fit_gwr(G, Ygf, Zf, ngroup, D, u, bw, kernel, verbose, pcout);
        gamma = get<0>(gwr_result);
        bw_optim = get<1>(gwr_result);
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
        ML_Params ml_params = { Xf, Yhf, Zf, &beta, ngroup, ndata, nvx, nvz };
        switch (ml_type)
        {
        case 0:
            mlf = fit_D(D, &ml_params, alpha, eps_gradient, max_iters, verbose > 1, pcout);
            beta = fit_gls(Xf, Yhf, Zf, ngroup, D);
            break;
        case 1:
            ml_params.beta = nullptr;
            beta = fit_gls(Xf, Yhf, Zf, ngroup, D);
            fit_D_beta(D, beta, &ml_params, alpha, eps_gradient, max_iters, verbose > 1, pcout);
            break;
        default:
            fit_D(D, &ml_params, alpha, eps_gradient, max_iters, verbose > 1, pcout);
            beta = fit_gls(Xf, Yhf, Zf, ngroup, D);
            break;
        }
        mu = fit_mu(Xf, Yhf, Zf, ngroup, beta, D);
        //------------------------------
        // Calculate Termination Measure
        //------------------------------
        vec yhat = yh - (X * beta) - sum(Z % (mu.rows(group)), 1);
        vec residual = yhat % yhat;
        rss = sum(residual);
        diff = rss - rss0;
        if (rss < rss0) 
        {
            if (retry > 0) retry = 0;
        }
        else if (iter > 0) retry++;
        if (verbose > 0)
        {
            ostringstream sout;
            sout << fixed << setprecision(prescition) << "Iter: " << iter;
            if (bw == 0.0) sout << ", " << "Bw: " << bw_optim;
            sout << ", " << "RSS: " << rss;
            if (abs(diff) < DBL_MAX) sout << ", " << "dRSS: " << diff;
            sout << ", " << "R2: " << (1 - rss / tss);
            sout << ", " << "-loglik/n: " << mlf;
            if (retry > 0) sout << ", " << "Retry: " << retry;
            sout << endl;
            pcout(sout.str());
        }
    }
    double sigma = fit_sigma(Xf, Yhf, Zf, ngroup, ndata, beta, D);
    delete[] Zf;
    delete[] Xf;
    delete[] Yf;
    delete[] Ygf;
    delete[] Yhf ;
    return { gamma, beta, mu, D, sigma, bw_optim };
}
