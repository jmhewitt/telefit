/*
	basic tools to implement models and methods from Rue and Held
 */

// disable assertions
#define EIGEN_NO_DEBUG

#include <RcppEigen.h>

// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<double> EigenSpMat;


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
		evaluate glm log-likelihood given covariate and random effects information

		Parameters:
		 y - observations
		 eta0 - GLM random effects
		 beta - values of beta (p x 1)
		 x - covariate matrix (nt x p)
		 n, t - number of locations and timepoints
		 p - the dimension of p (i.e., number of coefficients in beta)
		 family - specification of likelihood family
	*/
	double ll(const double* y, const double* eta0, const double* beta,
		const double* x, int n, int t, int p, const glmfamily family);

    /*
     Implement extension of Rue and Held (2005) Section 4.4.1.  Compute a
     Gaussian approximation to the posterior
      f(x|w) \propto \exp(-.5 (x-mu)^T Q (x-mu) + l(x;w))
     where mu=0, Q is sparse, and l(x;w) has a l(x;w) has a quadratic Taylor
     expansion that includes a diagonal Hessian matrix, C.  A sequence of Taylor
     expansions are used to try to derive a Gaussian approximation to f(x|w)
     that is centered near the mode of f(x|w).

     Specifically, this function is tailored to log-likelihoods l(x;w) with the
     structure used in the RESP GLM model.

     Parameters:

       (to build approximation)

       eta0 - initial value of eta0 around which to base expansion (ntx1)
       Q - sparse prior precision matrix (nt x nt)
       it - number of Taylor expansions to use in developing approximation
       mu - (output) mean of Gaussian approximation
       prec_chol - (output) cholesky for precision matrix of approximation

       (to Taylor-expand the likelihood)

       beta - values of beta (p x 1)
       y - observations (nt x 1)
       x - covariate matrix (nt x p)
       n, t - number of locations and timepoints
			 p - the dimension of p (i.e., number of coefficients in beta)
    */
    void gaussian_approx_eta0(const double* eta0, EigenSpMat& Q, int it,
			double* mu, Eigen::SimplicialLLT<EigenSpMat>& prec_chol,
			const double* beta, const double* y, const double* x, int n, int t, int p,
			const glmfamily family);

  /*
   Implement extension of Rue and Held (2005) Section 4.4.1.  Compute a
   Gaussian approximation to the posterior
    f(x|w) \propto \exp(-.5 (x-mu)^T Q (x-mu) + l(x;w))
   where mu=0, Q is diagonal, and l(x;w) has a l(x;w) has a quadratic Taylor
   expansion that includes a dense Hessian matrix, C.  A sequence of Taylor
   expansions are used to try to derive a Gaussian approximation to f(x|w) that
   is centered near the mode of f(x|w).

   Specifically, this function is tailored to log-likelihoods l(x;w) with the
   structure used in the RESP GLM model.

   Parameters:

     (to build approximation)

     beta0 - initial value of beta around which to base expansion (px1)
     Q - diagonal elements of prior covariance matrix for beta (px1)
     p - the dimension of beta (i.e., number of coefficients in beta)
     it - number of Taylor expansions to use in developing approximation
     mu - (output) mean of Gaussian approximation
     prec_chol - (output) cholesky for precision matrix of approximation

     (to Taylor-expand the likelihood)

     eta0 - values of eta_{0ij} (nt x 1)
     y - observations (nt x 1)
     x - covariate matrix (nt x p)
     n, t - number of locations and timepoints
  */
  void gaussian_approx_beta(const double* beta0, const double* Q, int p, int it,
    double* mu, Eigen::LLT<MatrixXd>& prec_chol, const double* eta0,
    const double* y, const double* x, int n, int t, const glmfamily family);

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
    x - covariate matrix (nt x p)
    y - observations (nt x 1)
    n, t - number of locations and timepoints
    p - the dimension of beta (i.e., number of coefficients in beta)
  */
  void glm_taylor_beta(double* b, double* c, const double* beta0,
    const double* eta0, const double* x, const double* y, int n, int t, int p,
    const glmfamily family);


  /*
   Implement extension of Rue and Held (2005) Section 4.4.1.  Compute the
   refactored version of the Taylor expansion of the log-likelihood for
   a GLM likelihood
      f(Y | beta, eta_0) =
        \prod_{i=1}^n \prod_{i=1}^t f(y_{ij} | g^{-1}(eta_{ij}))
   where eta_{ij} = x_{ij}^T beta + eta_{0ij}.

   This function returns b and C such that
    log f(Y | beta, eta_0) \approx a_0 + b^T eta_0 - .6 eta_0^T C eta_0,
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
