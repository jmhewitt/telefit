#ifndef _DISTRIBUTIONS_H
#define _DISTRIBUTIONS_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>

namespace mcstat2 {
	
	using namespace Rcpp;
	using namespace arma;
	
	using Eigen::PermutationMatrix;
	using Eigen::SparseMatrix;
	using Eigen::Dynamic;
	using Eigen::MatrixXd;
	
	
	//
	// densities
	//
	
	// log inverse gamma density without normalization constant
	double logdinvgamma_unscaled(double x, double a, double b);
	
	// log beta density without normalization constant
	double logdbeta_unscaled(double x, double a, double b);
	
	// log log-normal distribution without normalization constant
	double loglognormal_unscaled(double x, double mu, double sigma);
	
	
	//
	// samplers
	//
	
	// sample a mean 0 multivariate normal r.v.
	vec mvrnorm(const mat & sigma);
	
	// sample a mean 0 multivariate normal r.v. when given chol(sigma)
	vec mvrnorm_chol(const mat & L);
	
	// sample a mean 0 multivariate normal r.v. with kronecker covariance AxB
	// when given upper cholesky decompositions for A and B
	vec mvrnorm_cholkron(const mat & Ra, const mat & Rb);
	
	/* sample a mean 0 multivariate normal r.v. when given the cholesky
	 decomposition for a sparse precision matrix
	 */
	vec mvrnorm_spchol(const sp_mat & Linv);
	
	// sample a mean 0 multivariate normal r.v. when given solver for sparse precision
	vec mvrnorm_spchol(const SparseMatrix<double> &QL,
					   const PermutationMatrix<Dynamic,Dynamic> &QPinv,
					   int n);
	
	// sample N(0, V \otimes U) when given solvers for sparse precision V^-1 = Q
	// and U
	vec mvrnorm_spcholkron(const SparseMatrix<double> &QL,
						   const PermutationMatrix<Dynamic,Dynamic> &QPinv,
						   int Qn,
						   const MatrixXd &UA,
						   int Un);
	
	// sample a multivariate normal r.v.
	vec mvrnorm(const vec & mu, const mat & sigma);
	
	// sample a multivariate normal r.v., optionally when given chol(sigma)
	vec mvrnorm_chol(const vec & mu, const mat & L);
	
	// sample a multivariate normal r.v. using a precision or covariance matrix
	vec mvrnorm(const vec &mu, const mat &sigma, bool precision);
	
	// sample a wishart matrix using bartlett's decompostion
	mat rwishart(const mat &V, double n);
	
	// sample an inverse wishart matrix
	mat rinvwishart(const mat &V, double n);
	
	// sample a mean-zero multivariate normal r.v. x ~ N(0, Sigma).
	// Set precision=true if Q = inv(Sigma) is provided instead.
	mat mvrnorm(mat & Sigma, int nSamples, bool precision);
	
	// sample a multivariate normal r.v. with a bayesian posterior-type form
	// x ~ N( Sigma * y, Sigma ).  Set precision=true if Q = inv(Sigma) is
	// provided instead.
	mat mvrnorm_post(vec & y, mat & Sigma, int nSamples, bool precision);
	
	// sample a multivariate normal r.v. with both a kronecker and bayesian
	// posterior-type form x ~ N( Sigma * y, Sigma ) with Sigma = A x B.
	// Set precision=true if Qa = inv(A), Qb = inv(B) is provided instead.
	mat mvrnorm_postKron(vec & y, mat & A, mat & B, int nSamples,
						 bool precision);
	
}



#endif
