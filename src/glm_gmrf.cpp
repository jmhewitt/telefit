#include "glm_gmrf.h"

using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;


void mcstat2::glm::gmrf_approx(double* b, double* c, const double* x0,
  const double* y, int n, const mcstat2::glm::glmfamily family) {

  switch(family) {
    case poisson: // exponential link function
      for(int i=0; i<n; i++) {
        double d2 = - std::exp(x0[i]);
        double d1 = y[i] + d2;
        b[i] = d1 - d2 * x0[i];
        c[i] = d2>0 ? 0 : - d2;
      }
      break;
  }
}


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
  Map<const MatrixXd> cmat(c,p,p);
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
      H *= -1;
      g -= H*bmat;

      break;
  }

}


//
// Rcpp exports
//

// [[Rcpp::depends(RcppEigen)]]

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
List test_gmrf_approx(NumericVector y, NumericVector x0) {

  int n = y.size();

  // initialize results
  NumericVector b(n);
  NumericVector c(n);

  //
  // get access to raw data
  //

  std::vector<double> y_v = Rcpp::as<std::vector<double> >(y);
  std::vector<double> x_v = Rcpp::as<std::vector<double> >(x0);
  std::vector<double> b_v = Rcpp::as<std::vector<double> >(b);
  std::vector<double> c_v = Rcpp::as<std::vector<double> >(c);
  double* b_data = &b_v[0];
  double* c_data = &c_v[0];
  double* y_data = &y_v[0];
  double* x_data = &x_v[0];

  // build approximation elements
  mcstat2::glm::gmrf_approx(b_data, c_data, x_data, y_data, n,
    mcstat2::glm::glmfamily::poisson);

  // extract approximation
  for(int i=0; i<n; i++) {
    b[i] = b_data[i];
    c[i] = c_data[i];
  }

  return List::create(Named("b") = b, Named("c") = c);
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
