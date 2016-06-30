#include "covs.h"

using namespace arma;


// build a matern covariance matrix in place
//
//  WARNING: this function assumes cov is a distance matrix---symmetric, and
//			 with diag(cov)=0.  no checks are made.
//
void maternCov(mat & cov, const mat & d, double scale, double range,
			   double smoothness, double nugget ) {
	
	double cst = pow(2.0, 1.0 - smoothness) / R::gammafn(smoothness);
	double cstInv = 1.0 / cst;
	
	// compute elementwise correlations
	int n = cov.n_rows;
	for(int i=0; i<n; i++) {
		
		// diagonal
		cov.at(i,i) = cstInv;
		
		// off-diagonal
		for(int j=0; j<i; j++) {
			double v = d.at(i,j) / range;
			cov.at(i,j) = pow(v, smoothness) * R::bessel_k(v, smoothness, 1.0);
			cov.at(j,i) = cov.at(i,j);
		}
	}
	
	// scale to covariances
	cov = scale * cst * cov;
	
	// add nugget
	if(nugget!=0)
		cov.diag() += nugget;
}



//
// R exports
//


RcppExport SEXP _maternCov(SEXP d, SEXP scale, SEXP range, SEXP smoothness,
						 SEXP nugget) {
	
	using namespace Rcpp;
	
	mat dist = as<mat>(d);
	mat res = mat(dist.n_rows, dist.n_rows, fill::zeros);
	
	maternCov( res, dist, as<double>(scale), as<double>(range),
			   as<double>(smoothness), as<double>(nugget) );
	
	return wrap(res);
}