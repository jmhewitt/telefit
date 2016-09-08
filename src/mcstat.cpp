#include "mcstat.h"

using namespace arma;


// evaluate kron(A,B) * C without storing kron(A,B)
mat mcstat::dgemkmm(mat A, mat B, mat C) {
	int m = A.n_rows;
	int n = A.n_cols;
	int p = B.n_rows;
	int q = B.n_cols;
	int r = C.n_cols;
	
	mat res;
	res = mat(m*p, r, fill::zeros);
	
	mat resBlock;
	resBlock.set_size(q, r);
	for(int i=0; i<m; i++) {
		resBlock.zeros();
		
		for(int j=0; j<n; j++)
			resBlock += A.at(i,j) * C.rows( j*q, (j+1)*q-1 );
		
		res.rows( i*p, (i+1)*p-1 ) += B * resBlock;
	}
	
	return res;
}


// evaluate kron(A,B) * C without storing kron(A,B)
mat mcstat::dsemkmm(mat A, mat B, SpMat<double> C) {
	int m = A.n_rows;
	int n = A.n_cols;
	int p = B.n_rows;
	int q = B.n_cols;
	int r = C.n_cols;
	
	mat res;
	res = mat(m*p, r, fill::zeros);
	
	mat resBlock;
	resBlock.set_size(q, r);
	for(int i=0; i<m; i++) {
		resBlock.zeros();
		
		for(int j=0; j<n; j++)
			resBlock += A.at(i,j) * C.rows( j*q, (j+1)*q-1 );
		
		res.rows( i*p, (i+1)*p-1 ) += B * resBlock;
	}
	
	return res;
}


//
// R exports
//


RcppExport SEXP _dgemkmm(SEXP A, SEXP B, SEXP C) {
	
	using namespace Rcpp;
	
	return wrap( mcstat::dgemkmm(as<mat>(A), as<mat>(B), as<mat>(C)) );
}


RcppExport SEXP _rwishart(SEXP V, SEXP n) {
	
	using namespace Rcpp;
	
	return wrap( mcstat::rwishart(as<mat>(V), as<int>(n)) );
}


RcppExport SEXP _rinvwishart(SEXP V, SEXP n) {
	
	using namespace Rcpp;
	
	return wrap( mcstat::rinvwishart(as<mat>(V), as<int>(n)) );
}