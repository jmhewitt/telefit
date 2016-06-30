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
	
	
	//
	// utility functions
	//

	
	/// call 'message' from C++
	// inline Rcpp::Function msg(\"message\");
	// usage: msg(std::string(\"Hello, txt: \") + txt);
}


#endif