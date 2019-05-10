#ifndef _DISTRIBUTIONS_H
#define _DISTRIBUTIONS_H

// disable assertions
#define EIGEN_NO_DEBUG

#include <RcppArmadillo.h>
#include <RcppEigen.h>

namespace mcstat2 {

	using namespace Rcpp;
	using namespace arma;

	using Eigen::PermutationMatrix;
	using Eigen::SparseMatrix;
	using Eigen::Dynamic;
	using Eigen::MatrixXd;
	using Eigen::VectorXd;
	using Eigen::LLT;
	using Eigen::Lower;
	using Eigen::SimplicialLLT;

	// declares a column-major sparse matrix type of double
	typedef Eigen::SparseMatrix<double> EigenSpMat;

	//
	// densities
	//

	// log inverse gamma density without normalization constant
	double logdinvgamma_unscaled(double x, double a, double b);

	// log beta density without normalization constant
	double logdbeta_unscaled(double x, double a, double b);

	// log log-normal distribution without normalization constant
	double loglognormal_unscaled(double x, double mu, double sigma);

  // log inverse gamma density with normalization constant
  double ldinvgamma(double x, double a, double b);

  // log uniform density with normalization constant
  double ldunif(double x, double a, double b);

  /*
   log density for intrinsic GMRF with mean 0 and precision given by
	 kron(Q1, Q2).  Q1 is a sparse n x n matrix with rank n-k.  Q1 takes the
	 form q * R, where q is a scaling constant and R is the sparse structural
	 form for Q1.  Q2 is an m x m precision matrix for a dense covariance
	 matrix Sigma.

	 Reference: Basic IGMRF density is eqn. 3.18 in Rue and Held (2005)

	 Parameters:
	   x - vector at which to evaluate density
		 m - dimension of Q2
		 n - dimension of Q1
		 k - rank deficiency of Q1
		 R - structural form of Q1
		 q - scale for Q1
		 ldetR - log of the generalized determinant of R
		 L - cholesky decomposition of Sigma
   */
  double ldigmrfSpD(double* x, int m, int n, int k,
		const SparseMatrix<double>& R, double q, double ldetR,
		const LLT<MatrixXd, Lower>& L);

	/*
	 log density for intrinsic GMRF with mean 0 and precision given by Q.  Q takes
	 the form z * R, where q is a scaling constant and R is the sparse structural
	 form for Q.

	 Reference: Basic IGMRF density is eqn. 3.18 in Rue and Held (2005

	 Parameters:
	  x - vector at which to evaluate density
		n - dimension of Q
		k - rank deficiency of Q
		R - structural form of Q
		q - scale for Q
		ldetR - log of the generalized determinant of R
	*/
	double ldigmrfSp(double* x, int n, int k, const SparseMatrix<double>& R ,
		double q, double ldetR);


	/*
		log density for intrinsic GMRF with mean mu and kronecker-structured
		precision given by Qa x Qb.  Qa is the precision for a dense covariance
		matrix, and Qb = k * R where k is a constant and R is the structural form
		for an intrinsic GMRF dependence structure.

		Parameters:
			x - vector at which to evaluate density
			mu - density location parameter
			lltSigmaA - LLT decomposition for SigmaA = inv(Qa)
			k - scale for Qb
			df - rank deficiency of Qb
			Qab - prior precision matrix Qa x Qb
	*/
	double ldigmrfKron(const VectorXd& x, const VectorXd& mu,
		const LLT<MatrixXd>& lltSigmaA, double k, int df,
		const SparseMatrix<double> Qab);


	/*
	 	log density for a multivariate normal r.v. x ~ N(mu, Sigma) where
	  1) Sigma has a sparse precision matrix Q
	  2) An Eigen::SimplicialLLT decomposition of Q is provided

		Parameters:
		  x - value at which to evaluate log density
			mu - mean
			prec - Eigen::SimplicialLLT decomposition of Q
	*/
	double ldmvrnorm_spchol(const VectorXd& x, const VectorXd& mu,
		const SimplicialLLT<EigenSpMat>& prec);

	/*
	  log density for a multivariate normal r.v. x ~ N(mu, Sigma) where
		an Eigen::LLT decomposition of inv(Sigma) is provided
	*/
	double ldmvrnorm_chol(const VectorXd& x, const VectorXd& mu,
		const LLT<MatrixXd>& prec);


	//
	// probabilities
	//

	vec qintnorm(const vec & breaks, double mu, double sigma);


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

	// sample a mean zero multivariate normal r.v. x ~ N(0, Sigma) where
	//  1) Sigma has a sparse precision matrix Q
	//  2) An Eigen::SimplicialLLT decomposition of Q is provided
	VectorXd mvrnorm_spchol(const SimplicialLLT<EigenSpMat>& prec);

	// sample a mean zero multivariate normal r.v. x ~ N(0, Sigma) where
	//  an Eigen::LLT decomposition of inv(Sigma) is provided
	VectorXd mvrnorm_precchol(const LLT<MatrixXd>& prec_chol);

}



#endif
