#ifndef _telefit_MCSTAT_H
#define _telefit_MCSTAT_H


#include <RcppArmadillo.h>


namespace mcstat {
	
	using namespace arma;
	
	
	//
	// efficient linear algebra routines
	//
	
	// evaluate kron(A,B) * C without storing kron(A,B)
	mat dgemkmm(mat A, mat B, mat C);
	mat dsemkmm(mat A, mat B, SpMat<double> C);
	
	
	//
	// transformations
	//
	
	inline double logit(double x) { return log( x / (1.0 - x) ); }
	inline double invlogit(double x) {
		double expX = exp(x);
		return isinf(expX)!=0 ? 1 : expX / (1.0 + expX);
	}
	

	//
	// densities
	//
	
	// log inverse gamma density, without normalization constant
	inline double logdinvgamma_unscaled(double x, double a, double b) {
		return - (a + 1.0) * log(x) - b/x;
	}
	
	
	//
	// proposal functions
	//
	
	inline double logitProposal(double x, double min_x, double max_x, double sd) {
		double w = max_x - min_x;
		return invlogit( logit((x - min_x)/w) + R::rnorm(0, sd) ) * w + min_x;
	}
	
	inline double logProposal(double x, double sd) {
		return exp( log(x) + R::rnorm(0, sd) );
	}
	
	
	//
	// proposal jacobians
	//
	
	inline double loglogJacobian(double x) { return -log( std::abs(x) ); }
	inline double loglogitJacobian(double x) { return -log( std::abs(x*(1.0-x)) ); }
	
	
	//
	// samplers
	//
	
	// sample a mean 0 multivariate normal r.v.
	inline vec mvrnorm(const mat & sigma) {
		return chol(sigma, "lower") * randn<vec>(sigma.n_rows, 1);
	}
	
	// sample a multivariate normal r.v.
	inline vec mvrnorm(const vec & mu, const mat & sigma) {
		return mu + chol(sigma, "lower") * randn<vec>(sigma.n_rows, 1);
	}
	
	// sample a wishart matrix using bartlett's decompostion
	inline mat rwishart(mat V, int n) {
		// Params:
		//  n - degrees of freedom
		//  V - (symmetric) scale matrix
		
		int p = V.n_rows;
		mat A = mat(p, p, fill::zeros);
		
		// fill diagonal
		for(int i=0; i<p; i++)
			A(i,i) = sqrt(R::rchisq(n-i));
		
		// fill lower-triangular portion of matrix
		for(int i=1; i<p; i++)
			for(int j=0; j<i; j++)
				A(i,j) = R::rnorm(0,1);
		
		mat C = chol(V, "lower") * trimatl(A);
		return C * C.t();
	}
	
	// sample an inverse wishart matrix
	inline mat rinvwishart(mat V, int n) {
		return inv_sympd( rwishart(inv_sympd(V), n) );
	}
	
	
	//
	// utility functions
	//

	
	/// call 'message' from C++
	// inline Rcpp::Function msg(\"message\");
	// usage: msg(std::string(\"Hello, txt: \") + txt);
}


#endif