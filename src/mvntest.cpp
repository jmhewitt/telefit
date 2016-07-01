#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;


RcppExport SEXP _mvnShort (SEXP L_, SEXP R_) {

	
	// draw from a matrix normal
	
	
	mat L = as<mat>(L_);
	mat R = as<mat>(R_);
	
	int r = R.n_rows;
	int ns = L.n_rows;
	
	
	return(wrap( vectorise(chol(R, "lower") * randn<mat>(r,ns) * chol(L, "upper")) ));
}


RcppExport SEXP _mvnLong (SEXP L_, SEXP R_) {
	
	// draw from a kronecker

	
	mat L = as<mat>(L_);
	mat R = as<mat>(R_);
	
	int r = R.n_rows;
	int ns = L.n_rows;
	
	
	return(wrap( kron(chol(L, "lower"), chol(R, "lower")) * randn<vec>(r*ns) ));
}