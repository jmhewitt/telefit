/*
	basic tools to implement models and methods from Rue and Held
 */

// disable assertions
#define EIGEN_NO_DEBUG

#include <RcppEigen.h>


namespace mcstat2 {
namespace glm {

	using namespace Rcpp;

  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::Map;

  // supported glm families
  enum glmfamily { poisson };

  // evaluate glm log-likelihood
  //
  // Parameters:
  //  y - observations
  //  eta - GLM means
  //  n - number of observations
  //  family - specification of likelihood family
  double ll(const double* y, const double* eta, const int n,
    const glmfamily family);

  /*
   Follow Rue and Held (2005) Section 4.4.1 to construct a GMRF approximation
   for the posterior distribution of a latent GMRF with respect to a
   conditionally independent GLM response.  This function only focuses on the
   quadratic expansion of the GLM likelihood function so that the approximation
   can be paired with any number of latent structures.

   Parameters:
     b - (pre-initialized output) array in which to store b vector in eqn. 4.27
     c - (pre-initialized output) array in which to store c vector in eqn. 4.27
     x0 - values around which the Taylor approximation should be centered
     y - observations
     n - GMRF dimension
  */
  void gmrf_approx(double* b, double* c, const double* x0, const double* y,
    int n, const glmfamily family);

  /*
   Implement extension of Rue and Held (2005) Section 4.4.1.  Compute the
   refactored version of the Taylor expansion of the log-likelihood for
   a GLM likelihood
      f(Y | beta, eta_0) =
        \prod_{i=1}^n \prod_{i=1}^t f(y_{ij} | g^{-1}(eta_{ij}))
   where eta_{ij} = x_{ij}^T beta + eta_{0ij}.

   This function returns b and C such that
    log f(Y | beta, eta_0) \approx a_0 + b^T beta - .6 beta^T C beta,
   i.e., the refactored Taylor expansion of the log-likelihood around beta_0.

   Parameters:
    b - (pre-initialized output) array in which to store b vector (p x 1)
    c - (pre-initialized output) array in which to store C matrix in
     column-major format (p x p)
    beta0 - values around which Taylor approximation should be centered (nt x 1)
    eta0 - values of eta_{0ij} (nt x 1)
    y - observations (nt x 1)
    x - covariate matrix (nt x p)
    n, t - number of locations and timepoints
    p - the dimension of p (i.e., number of coefficients in beta)
  */
  void glm_taylor_beta(double* b, double* c, const double* beta0,
    const double* eta0, const double* x, const double*y, int n, int t, int p,
    const glmfamily family);


  /*
   Implement extension of Rue and Held (2005) Section 4.4.1.  Compute the
   refactored version of the Taylor expansion of the log-likelihood for
   a GLM likelihood
      f(Y | beta, eta_0) =
        \prod_{i=1}^n \prod_{i=1}^t f(y_{ij} | g^{-1}(eta_{ij}))
   where eta_{ij} = x_{ij}^T beta + eta_{0ij}.

   This function returns b and C such that
    log f(Y | beta, eta_0) \approx a_0 + b^T eta_0 - .6 beta^T C eta_0,
   i.e., the refactored Taylor expansion of the log-likelihood around eta_0.

   Note: The Hessian for this expansion is a diagonal matrix.

   Parameters:
    b - (pre-initialized output) array in which to store b vector (nt x 1)
    c - (pre-initialized output) array in which to store diagonal hessian (nt)
    beta - values of beta (p x 1)
    eta0 - values of eta_{0ij} around which Taylor approx. is centered (nt x 1)
    y - observations (nt x 1)
    x - covariate matrix (nt x p)
    n, t - number of locations and timepoints
    p - the dimension of p (i.e., number of coefficients in beta)
  */
  void glm_taylor_eta0(double* b, double* c, const double* beta,
    const double* eta0, const double* x, const double*y, int n, int t, int p,
    const glmfamily family);

}}
