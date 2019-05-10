#include "glm_gmrf.h"

using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Dynamic;
using Eigen::DiagonalMatrix;
using Eigen::LLT;
using Eigen::SimplicialLLT;


double mcstat2::glm::ll(const double* y, const double* eta, const int n,
    const mcstat2::glm::glmfamily family) {

    double res = 0;

    switch(family) {
        case poisson:  // exponential link function
            for(int i=0; i<n; i++)
                res += y[i] * eta[i] - std::exp(eta[i]) - std::lgamma(y[i] + 1);
            break;
    }

    return res;
}

void mcstat2::glm::glm_taylor_beta(double* b, double* c, const double* beta0,
  const double* eta0, const double* x, const double* y, int n, int t, int p,
  const glmfamily family) {

  // map inputs/outputs
  int nt = n*t;
  Map<const MatrixXd> xmat(x,nt,p);
  Map<const MatrixXd> eta0mat(eta0,n,t);
  Map<const VectorXd> bmat(beta0,p);
  Map<const VectorXd> ymat(y,nt);
  Map<MatrixXd> H(c,p,p);
  Map<VectorXd> g(b,p);

  // initialize derivatives; they will be built sequentially
  g.setZero();
  H.setZero();

  switch (family) {
    case poisson:

      // pre-compute exponential term
      VectorXd m = xmat * bmat;
      Map<MatrixXd> eta(m.data(),n,t);
      for(int i=0; i<n; i++) {
        for(int j=0; j<t; j++) {
          eta(i,j) = std::exp(eta(i,j) + eta0mat(i,j));
        }
      }

      // loop over timepoints
      for(int j=0; j<t; j++) {

        // extract covariate matrix and data at time j
        MatrixXd Xj = xmat.block(j*n, 0, n, p);
        VectorXd Yj = ymat.segment(j*n, n);

        // loop over locations
        for(int i=0; i<n; i++) {

          // primary loop over regression coefficients
          for(int l=0; l<p; l++) {

            // build gradient
            g(l) += Xj(i,l) * (Yj(i) - eta(i,j));

            // secondary loop over regression coefficients
            for(int m=0; m<p; m++) {

              // build negative hessian
              H(l,m) += Xj(i,l) * Xj(i,m) * eta(i,j);

            }
          }
        }
      }

      // finish objects
      g += H*bmat;

      break;
  }

}

void mcstat2::glm::glm_taylor_eta0(double* b, double* c, const double* beta,
  const double* eta0, const double* x, const double* y, int n, int t, int p,
  const glmfamily family) {

  // map inputs/outputs
  int nt = n*t;
  Map<const MatrixXd> xmat(x,nt,p);
  Map<const MatrixXd> eta0mat(eta0,n,t);
  Map<const VectorXd> bmat(beta,p);
  Map<const VectorXd> ymat(y,nt);
  Map<MatrixXd> H(c,n,t);
  Map<MatrixXd> g(b,n,t);


  switch (family) {
    case poisson:

      // pre-compute exponential term
      VectorXd m = xmat * bmat;
      Map<MatrixXd> eta(m.data(),n,t);
      for(int i=0; i<n; i++) {
        for(int j=0; j<t; j++) {
          eta(i,j) = std::exp(eta(i,j) + eta0mat(i,j));
        }
      }

      // loop over timepoints
      for(int j=0; j<t; j++) {

        // extract data at time j
        VectorXd Yj = ymat.segment(j*n, n);

        // loop over locations
        for(int i=0; i<n; i++) {

            // build gradient
            g(i,j) = Yj(i) - eta(i,j);

            // build negative of hessian
            H(i,j) = eta(i,j);
        }
      }

      // finish objects
      Map<VectorXd> Hm(c,nt);
      Map<VectorXd> gm(b,nt);
      Map<const VectorXd> eta0vec(eta0,nt);
      gm += Hm.asDiagonal() * eta0vec;

      break;
  }

}

void mcstat2::glm::gaussian_approx_beta(const double* beta0, const double* Q,
  int p, int it, double* mu, Eigen::LLT<MatrixXd>& prec_chol,
  const double* eta0, const double* y, const double* x, int n, int t,
  const glmfamily family) {

  // map inputs/outputs
  Map<const VectorXd> Qvec(Q, p);
  Map<VectorXd> mvec(mu, p);

  // temporary objects
  memcpy(mu, beta0, sizeof(double)*p);
  VectorXd b(p);
  MatrixXd C(p,p);

  // iterate to find mean and precision close to mode of target density
  for(int i=0; i<it; i++) {
    // Taylor-expand likelihood
    glm_taylor_beta(b.data(), C.data(), mu, eta0, x, y, n, t, p, family);
    // factor approximation's precision matrix
    C.diagonal() += Qvec;
    prec_chol.compute(C);
    // compute approximation's mean
    mvec = prec_chol.solve(b);
  }
}

void mcstat2::glm::gaussian_approx_eta0(const double* eta0,
  EigenSpMat& Q, int it, double* mu,
  Eigen::SimplicialLLT<EigenSpMat>& prec_chol, const double* beta,
  const double* y, const double* x, int n, int t, int p,
  const glmfamily family) {

  // map inputs/outputs
  int nt = n*t;
  Map<VectorXd> mvec(mu, nt);

  // temporary objects
  memcpy(mu, eta0, sizeof(double)*nt);
  VectorXd b(nt);
  VectorXd C(nt);
  EigenSpMat Qm(Q);

/*
  // TODO: these steps should be common across all calls to the function, so
  // as an optimization routine, we could consider modifying the function so
  // that these values are passed in

  std::vector<double*> diag_Qm(Qm.rows(), NULL);
  std::vector<double*> diag_Q(Qm.rows(), NULL);
  for(int i=0; i<Qm.rows(); i++) {
    diag_Qm[i] = &Qm.coeffRef(i,i);
    diag_Q[i] = &Q.coeffRef(i,i);
  }

  prec_chol.analyzePattern(Qm);

  // update the diagonal entry list
  // check to see if entries are really moved instead?
  // avoiding this check may have diminishing returns
  for(int i=0; i<Qm.rows(); i++) {
    diag_Qm[i] = &Qm.coeffRef(i,i);
  }

  // we spend a lot of time in this step.  coeffref and
  // analyze take up half of the remaining time.  additionally, it looks like
  // the call to factorize may not be recognizing the decomposition

  // END TODO
*/

  // iterate to find mean and precision close to mode of target density
  for(int i=0; i<it; i++) {
    // Taylor-expand likelihood
    glm_taylor_eta0(b.data(), C.data(), beta, mu, x, y, n, t, p, family);
    // factor approximation's precision matrix
    for(int j=0; j<Qm.rows(); j++)
      Qm.coeffRef(j,j) = Q.coeff(j,j) + C[j];
      //*diag_Qm[j] = *diag_Q[j] + C[j];
    prec_chol.compute(Qm);
    // compute approximation's mean
    mvec = prec_chol.solve(b);
  }
}

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
double mcstat2::glm::ll(const double* y, const double* eta0, const double* beta,
  const double* x, int n, int t, int p, const glmfamily family) {

  // map inputs/outputs
  int nt = n*t;
  Map<const MatrixXd> xmat(x,nt,p);
  Map<const VectorXd> eta0vec(eta0,nt);
  Map<const VectorXd> bmat(beta,p);

  // build linear predictor
  VectorXd eta = xmat * bmat + eta0vec;

  // evaluate log-likelihood
  return mcstat2::glm::ll(y, eta.data(), nt, family);
}


//
// Rcpp exports
//

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
List test_gaussian_approx_eta0(const Eigen::Map<Eigen::VectorXd> eta0,
  Eigen::SparseMatrix<double> Q, int it,
  const Eigen::Map<Eigen::VectorXd> beta, const Eigen::Map<Eigen::VectorXd> y,
  const Eigen::Map<Eigen::MatrixXd> x, int n, int t, int p) {

  // initialize output
  VectorXd mu(n*t);
  SimplicialLLT<EigenSpMat> prec_chol;

  mcstat2::glm::gaussian_approx_eta0(eta0.data(), Q, it, mu.data(),
    prec_chol, beta.data(), y.data(), x.data(), n, t, p,
    mcstat2::glm::glmfamily::poisson);

  EigenSpMat L = prec_chol.matrixL();
  return List::create(Named("mu")=mu, Named("L")=L);
}

// [[Rcpp::export]]
List test_gaussian_approx_beta(const Eigen::Map<Eigen::VectorXd> beta0,
  const Eigen::Map<Eigen::VectorXd> Q, int p, int it,
  const Eigen::Map<Eigen::VectorXd> eta0, const Eigen::Map<Eigen::VectorXd> y,
  const Eigen::Map<Eigen::MatrixXd> x, int n, int t) {

  // initialize output
  VectorXd mu(p);
  LLT<MatrixXd> prec_chol;

  mcstat2::glm::gaussian_approx_beta(beta0.data(), Q.data(), p, it, mu.data(),
    prec_chol, eta0.data(), y.data(), x.data(), n, t,
    mcstat2::glm::glmfamily::poisson);

  MatrixXd L = prec_chol.matrixL();
  return List::create(Named("mu")=mu, Named("L")=L);
}

// [[Rcpp::export]]
List test_taylor_beta(Eigen::Map<Eigen::MatrixXd> beta0,
                      Eigen::Map<Eigen::MatrixXd> eta0,
                      Eigen::Map<Eigen::MatrixXd> y,
                      Eigen::Map<Eigen::MatrixXd> x,
                      int n, int t, int p) {

  VectorXd b(p);
  MatrixXd c(p,p);

  mcstat2::glm::glm_taylor_beta(b.data(), c.data(), beta0.data(), eta0.data(),
    x.data(), y.data(), n, t, p, mcstat2::glm::glmfamily::poisson);

  return List::create(Named("b") = b, Named("C") = c);
}

// [[Rcpp::export]]
List test_taylor_eta0(Eigen::Map<Eigen::MatrixXd> beta,
                      Eigen::Map<Eigen::MatrixXd> eta0,
                      Eigen::Map<Eigen::MatrixXd> y,
                      Eigen::Map<Eigen::MatrixXd> x,
                      int n, int t, int p) {

  int nt = n*t;
  VectorXd b(nt);
  VectorXd c(nt);

  mcstat2::glm::glm_taylor_eta0(b.data(), c.data(), beta.data(), eta0.data(),
    x.data(), y.data(), n, t, p, mcstat2::glm::glmfamily::poisson);

  return List::create(Named("b") = b, Named("C") = c);
}


// [[Rcpp::export]]
NumericVector test_ll(NumericVector y, NumericVector lambda) {

	//
	// extract data
	//

	int n = y.size();

	std::vector<double> y_v = Rcpp::as<std::vector<double> >(y);
	std::vector<double> lambda_v = Rcpp::as<std::vector<double> >(lambda);

	double* y_data = &y_v[0];
	double* lambda_data = &lambda_v[0];

	// evaluate likelihood
	return wrap(mcstat2::glm::ll(y_data, lambda_data, n,
    mcstat2::glm::glmfamily::poisson));
}

// [[Rcpp::export]]
double test_ll_alt(Eigen::VectorXd y, Eigen::VectorXd eta0,
  Eigen::VectorXd beta, Eigen::MatrixXd x, int n, int t) {
  return mcstat2::glm::ll(y.data(), eta0.data(), beta.data(),
    x.data(), n, t, x.cols(), mcstat2::glm::glmfamily::poisson);
}
