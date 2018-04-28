#include "numAlg.h"

using namespace arma;


mat mcstat2::dgemkmm(mat A, mat B, mat C) {
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

mat mcstat2::dgeikmm(int N, mat A, mat B) {
	int m = A.n_rows;
	int n = A.n_cols;
	int p = B.n_cols;
	
	mat res;
	res = mat(N*m, p, fill::zeros);
	
	for(int i=0; i<N; i++)
		res.rows(i*m, (i+1)*m-1) = A * B.rows(i*n, (i+1)*n-1);
	
	return res;
}


//
// Rcpp Exports
//


RcppExport SEXP _dgeikmm(SEXP _N, SEXP _A, SEXP _B) {
	
	using namespace Rcpp;
	
	int N = as<int>(_N);
	mat A = as<mat>(_A);
	mat B = as<mat>(_B);
	
	return wrap(mcstat2::dgeikmm(N, A, B));
}
