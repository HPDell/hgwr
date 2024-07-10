#ifndef HLMGWR_H
#define HLMGWR_H

#include <string>
#include <memory>
#include <functional>
#include "armadillo_config.h"

namespace hgwr
{

struct ML_Params
{
    arma::mat* Xf;
    arma::vec* Yf;
    arma::mat* Zf;
    arma::vec* beta;
    arma::uword ngroup;
    arma::uword n;
    arma::uword p;
    arma::uword q;
};

class HGWR
{
public:  // Type defs
    enum class KernelType
    {
        GAUSSIAN,
        BISQUARED
    };

    typedef arma::vec (*GWRKernelFunctionSquared)(arma::vec, double);

    static double actual_bw(arma::vec d, double bw)
    {
        arma::vec ds = sort(d);
        double b0 = floor(bw), bx = bw - b0;
        double d0 = ds(int(b0) - 1), d1 = ds(int(b0));
        return d0 + (d1 - d0) * bx;
    }
    
    static arma::vec gwr_kernel_gaussian2(arma::vec dist2, double bw2)
    {
        return exp(- dist2 / (2.0 * bw2));
    }

    
    static arma::vec gwr_kernel_bisquare2(arma::vec dist2, double bw2)
    {
        return ((1 - dist2 / bw2) % (1 - dist2 / bw2)) % (dist2 < bw2);
    }


    /**
     * @brief Calculate $(I+LML^T)^{-1}$
     * 
     * @param M_inv The inverse matrix of $M$
     * @param L The matrix $L$
     * @return arma::mat The inverse $I - L(M^{-1}+L^TL)^-1 L^T$
     */
    static arma::mat woodbury_eye(const arma::mat& M_inv, const arma::mat& L)
    {
        return arma::eye<arma::mat>(L.n_rows, L.n_rows) - L * (M_inv + L.t() * L).i() * L.t();
    }

    typedef void (*PrintFunction)(const std::string&);
    
    static void Printer(const std::string& msg)
    {
        (void)msg;
    }

    struct Options
    {
        double alpha = 0.01;
        double eps_iter = 1e-6;
        double eps_gradient = 1e-6;
        size_t max_iters = (size_t)1e6;
        size_t max_retries = (size_t)10;
        size_t verbose = (size_t)0;
        size_t ml_type = (size_t)0;
    };

    struct Parameters
    {
        arma::mat gamma;
        arma::mat beta;
        arma::mat mu;
        arma::mat D;
        double sigma;
        double bw;
    };

    // using BwSelectionArgs = std::pair<std::reference_wrapper<arma::mat>, std::reference_wrapper<arma::vec>>;
    struct BwSelectionArgs
    {
        std::reference_wrapper<arma::mat> Vig;
        std::reference_wrapper<arma::vec> Viy;
        std::reference_wrapper<arma::mat> G;
        std::reference_wrapper<arma::mat> u;
        std::reference_wrapper<arma::rowvec> rVisigma;
        GWRKernelFunctionSquared kernel;
    };

    static double bw_criterion_cv(double bw, void* params);

    static double bw_criterion_aic(double bw, void* params);

public:
    explicit HGWR(const arma::mat& G, const arma::mat& X, const arma::mat& Z, const arma::vec& y, const arma::mat& u, const arma::uvec& group)
    {
        this->G = G;
        this->X = X;
        this->Z = Z;
        this->y = y;
        this->u = u;
        this->group = group;
        ngroup = G.n_rows;
        ndata = X.n_rows;
        nvg = G.n_cols;
        nvx = X.n_cols;
        nvz = Z.n_cols;
        bw_optim = true;
    }

    explicit HGWR(const arma::mat& G, const arma::mat& X, const arma::mat& Z, const arma::vec& y, const arma::mat& u, const arma::uvec& group, const Options& options)
        : HGWR(G, X, Z, y, u, group)
    {
        this->alpha = options.alpha;
        this->eps_iter = options.eps_iter;
        this->eps_gradient = options.eps_gradient;
        this->max_iters = options.max_iters;
        this->max_retries = options.max_retries;
        this->verbose = options.verbose;
        this->ml_type = options.ml_type;
        bw_optim = true;
    }

    explicit HGWR(const arma::mat& G, const arma::mat& X, const arma::mat& Z, const arma::vec& y, const arma::mat& u, const arma::uvec& group, KernelType kernel)
        : HGWR(G, X, Z, y, u, group)
    {
        this->kernel = kernel;
        bw_optim = true;
    }

    explicit HGWR(const arma::mat& G, const arma::mat& X, const arma::mat& Z, const arma::vec& y, const arma::mat& u, const arma::uvec& group, KernelType kernel, const Options& options)
        : HGWR(G, X, Z, y, u, group, options)
    {
        this->kernel = kernel;
        bw_optim = true;
    }

    explicit HGWR(const arma::mat& G, const arma::mat& X, const arma::mat& Z, const arma::vec& y, const arma::mat& u, const arma::uvec& group, KernelType kernel, double bw)
        : HGWR(G, X, Z, y, u, group, kernel)
    {
        this->bw = bw;
        bw_optim = false;
    }

    explicit HGWR(const arma::mat& G, const arma::mat& X, const arma::mat& Z, const arma::vec& y, const arma::mat& u, const arma::uvec& group, KernelType kernel, double bw, const Options& options)
        : HGWR(G, X, Z, y, u, group, options)
    {
        this->kernel = kernel;
        this->bw = bw;
        bw_optim = false;
    }

    explicit HGWR(const arma::mat& G, const arma::mat& X, const arma::mat& Z, const arma::vec& y, const arma::mat& u, const arma::uvec& group, KernelType kernel, double bw, const Options& options, const PrintFunction printer)
        : HGWR(G, X, Z, y, u, group, kernel, bw, options)
    {
        this->pcout = printer;
    }

public:
    
    const arma::mat& get_G() { return G; }
    void set_G(const arma::mat& value) { G = value; }
    
    const arma::mat& get_X() { return X; }
    void set_X(const arma::mat& value) { X = value; }
    
    const arma::mat& get_Z() { return Z; }
    void set_Z(const arma::mat& value) { Z = value; }
    
    const arma::vec& get_y() { return y; }
    void set_y(const arma::vec& value) { y = value; }
    
    const arma::mat& get_u() { return u; }
    void set_u(const arma::mat& value) { u = value; }
    
    const arma::uvec& get_group() { return group; }
    void set_group(const arma::uvec& value) { group = value; }
    
    double get_bw() { return bw; }
    void set_bw(double value)
    {
        bw = value;
        bw_optim = false;
    }

    bool get_bw_optim() { return bw_optim; }
    void set_bw_optim(bool value) { bw_optim = value; }
    
    KernelType get_kernel() { return kernel; }
    void set_kernel(KernelType value)
    {
        kernel = value;
        switch (kernel)
        {
        case KernelType::BISQUARED:
            gwr_kernel = &gwr_kernel_bisquare2;
            break;
        default:
            gwr_kernel = &gwr_kernel_gaussian2;
            break;
        }
    }

    double get_alpha() { return alpha; }
    void set_alpha(double value) { alpha = value; }
    
    double get_eps_iter() { return eps_iter; }
    void set_eps_iter(double value) { eps_iter = value; }
    
    double get_eps_gradient() { return eps_gradient; }
    void set_eps_gradient(double value) { eps_gradient = value; }
    
    size_t get_max_iters() { return max_iters; }
    void set_max_iters(size_t value) { max_iters = value; }
    
    size_t get_max_retries() { return max_retries; }
    void set_max_retries(size_t value) { max_retries = value; }
    
    size_t get_verbose() { return verbose; }
    void set_verbose(size_t value) { verbose = value; }
    
    size_t get_ml_type() { return ml_type; }
    void set_ml_type(size_t value) { ml_type = value; }
    
    arma::mat get_gamma() { return gamma; }
    
    arma::mat get_beta() { return beta; }
    
    arma::mat get_mu() { return mu; }
    
    arma::mat get_D() { return D; }
    
    double get_sigma() { return sigma; }

    double get_loglik() { return loglik; }

    arma::vec get_trS() { return trS; }

    arma::vec get_var_beta() { return var_beta; }

    void set_printer(PrintFunction printer) { pcout = printer; }

public:
    double criterion_bw(double bw, const BwSelectionArgs& args) { return 0.0; }
    int bw_optimisation(double lower, double upper, const BwSelectionArgs* args);
    double golden_selection(const double lower, const double upper, const bool adaptive, const BwSelectionArgs& args);
    void fit_gwr();
    arma::vec fit_gls();
    double fit_D(ML_Params* params);
    double fit_D_beta(ML_Params* params);
    void fit_mu();
    double fit_sigma();
    Parameters fit();
    void calc_var_beta();

private:
    /* data */
    arma::mat G;
    arma::mat X;
    arma::mat Z;
    arma::vec y;
    arma::mat u;
    arma::uvec group;
    double bw = 0.0;
    bool bw_optim = false;
    KernelType kernel = KernelType::GAUSSIAN;
    GWRKernelFunctionSquared gwr_kernel = &gwr_kernel_gaussian2;

    /* model parameters */
    arma::mat gamma;
    arma::vec beta;
    arma::mat mu;
    arma::mat D;
    double sigma;

    /* options */
    double alpha = 0.01;
    double eps_iter = 1e-6;
    double eps_gradient = 1e-6;
    size_t max_bw_iters = (size_t)100;
    size_t max_iters = (size_t)1e6;
    size_t max_retries = (size_t)10;
    size_t verbose = (size_t)0;
    size_t ml_type = (size_t)0;
    PrintFunction pcout = &Printer;


    /* others */
    std::unique_ptr<arma::mat[]> Zf;
    std::unique_ptr<arma::mat[]> Xf;
    std::unique_ptr<arma::vec[]> Yf;
    std::unique_ptr<arma::vec[]> Ygf;
    std::unique_ptr<arma::vec[]> Yhf;
    arma::uword ngroup;
    arma::uword ndata;
    arma::uword nvg;
    arma::uword nvx;
    arma::uword nvz;

    /* diagnostic information */
    double loglik = 0;
    arma::vec trS;
    arma::vec var_beta;
    double enp = 0;
    double edf = 0;
};
    
}

#endif  // HLMGWR_H
