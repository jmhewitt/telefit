#include "distributions.h"

using namespace arma;

using Eigen::PermutationMatrix;
using Eigen::Dynamic;
using Eigen::SparseMatrix;
using Eigen::MatrixXd;

double mcstat2::logdinvgamma_unscaled(double x, double a, double b) {
	return - (a + 1.0) * log(x) - b/x;
}

double mcstat2::logdbeta_unscaled(double x, double a, double b) {
	return (a-1) * log(x) + (b-1) * log(1-x);
}


double mcstat2::loglognormal_unscaled(double x, double mu, double sigma) {
	double logx = log(x);
	return - logx - log(sigma) - .5 * pow(logx - mu, 2.0) / pow(sigma, 2.0);
}

vec mcstat2::qintnorm(const vec & breaks, double mu, double sigma) {
	// compute normal probabilities for intervals specified by breaks
	
	int n = breaks.n_elem + 1;
	
	vec res = vec(n, fill::zeros);
	
	double F0 = 0;
	for(int i=0; i<(n-1); i++) {
		double F = Rf_pnorm5(breaks.at(i), mu, sigma, 1, 0);
		res.at(i) = F - F0;
		F0 = F;
	}
	res.at(n-1) = 1 - F0;
	
	return res;
}

vec mcstat2::mvrnorm(const mat & sigma) {
	return chol(sigma, "lower") * randn<vec>(sigma.n_rows, 1);
}

vec mcstat2::mvrnorm_chol(const mat & R) {
	// Note: This assumes that the cholesky is stored in uppertri
	// format for more efficient vector operations. it also uses in-place evaluation
	
	vec z = randn<vec>(R.n_rows, 1);
	double* z_out = z.memptr() + R.n_rows;
	double* z_mem = z.memptr();
	
	for(int i=R.n_rows-1; i>=0; i--){
		const double* R_mem = R.colptr(i);
		*(--z_out) *= *(R_mem + i);
		for(int j=i-1; j>=0; j--)
			*z_out += *(z_mem+j) * *(R_mem+j);
	}
	
	return z;
}

vec mcstat2::mvrnorm_cholkron(const mat & Ra, const mat & Rb) {
	//TODO: take advantage of Ra and Rb's triangular forms
	return vectorise( Rb.t() * randn<mat>(Rb.n_rows, Ra.n_rows) * Ra);
}

vec mcstat2::mvrnorm_spchol(const sp_mat & L) {
	// redo this by using a sparse backsolve after determining ordering
	// will likely need to use lapack routines for further efficiency
	return L * randn<vec>(L.n_rows, 1);
}

vec mcstat2::mvrnorm_spchol(const SparseMatrix<double> &QL,
						   const PermutationMatrix<Dynamic,Dynamic> &QPinv,
						   int n) {
	// generate independent normal variates
	vec _z = randn(n);
	
	// map variates into eigen
	using Eigen::Map;
	using Eigen::VectorXd;
	Map<VectorXd> z(_z.memptr(), n);
	
	// solve using sparse cholesky to generate sample
	using Eigen::Lower;
	VectorXd _x = QPinv * QL.triangularView<Lower>().transpose().solve(z);
	
	// convert sample to arma vec
	vec x = vec(_x.data(), n);
	
	return x;
}

vec mcstat2::mvrnorm_spcholkron(const SparseMatrix<double> &QL,
							   const PermutationMatrix<Dynamic,Dynamic> &QPinv,
							   int Qn,
							   const MatrixXd &UA,
							   int Un) {
	// generate independent normal variates
	mat _z = randn<mat>(Un, Qn);
	
	// map variates into eigen
	using Eigen::Map;
	using Eigen::MatrixXd;
	Map<MatrixXd> z(_z.memptr(), Un, Qn);
	
	// solve using sparse cholesky to generate sample
	using Eigen::Lower;
	MatrixXd _x = (QPinv * QL.triangularView<Lower>().transpose().solve(
								(UA.triangularView<Lower>() * z).transpose()
							   )).transpose();
	
	// convert sample to arma vec
	vec x = vec(_x.data(), Qn * Un);
	
	return x;
}

vec mcstat2::mvrnorm(const vec & mu, const mat & sigma) {
	return mu + chol(sigma, "lower") * randn<vec>(sigma.n_rows, 1);
}

vec mcstat2::mvrnorm_chol(const vec & mu, const mat & L) {
	return mu + L * randn<vec>(L.n_rows, 1);
}

vec mcstat2::mvrnorm(const vec &mu, const mat &sigma, bool precision) {
	if(precision) {
		return solve(chol(sigma, "upper"), randn<vec>(sigma.n_rows, 1)) + mu;
	} else {
		return mu + chol(sigma, "lower") * randn<vec>(sigma.n_rows, 1);
	}
}

mat mcstat2::rwishart(const mat &V, double n) {
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

mat mcstat2::rinvwishart(const mat &V, double n) {
	return inv_sympd( mcstat2::rwishart(inv_sympd(V), n) );
}


mat mcstat2::mvrnorm(mat & Sigma, int nSamples, bool precision) {
	
	using Eigen::MatrixXd;
	using Eigen::Map;
	using Eigen::LLT;
	using Eigen::Upper;
	using Eigen::Lower;
	
	int n = Sigma.n_rows;
	
	// generate independent normal variates
	GetRNGstate();
	mat _z = randn(n, nSamples);
	PutRNGstate();
	
	// map variates into eigen
	Map<MatrixXd> z(_z.memptr(), n, nSamples);
	
	// internal storage for r.v.s
	MatrixXd _x;
	
	if(precision) {
		
		// read precision matrix into eigen
		Map<MatrixXd> Q(Sigma.memptr(), n, n);
		
		// factor precision matrix; store upper cholesky
		LLT<MatrixXd, Upper> llt(Q);
		
		// compute r.v. with precision matrix Q
		_x = llt.matrixU().solve(z);
		
	} else {
		
		// read covariance matrix into eigen
		Map<MatrixXd> S(Sigma.memptr(), n, n);
		
		// factor covariance matrix; store lower cholesky
		LLT<MatrixXd, Lower> llt(S);
		
		// compute r.v. with covariance matrix S
		_x = llt.matrixL() * z;
	}
	
	// convert samples to arma vec
	mat x = mat(_x.data(), Sigma.n_rows, nSamples);
	return x;
}


mat mcstat2::mvrnorm_post(vec & y, mat & Sigma, int nSamples, bool precision) {
	
	using Eigen::MatrixXd;
	using Eigen::VectorXd;
	using Eigen::Map;
	using Eigen::LLT;
	using Eigen::Upper;
	
	if(precision) {
		
		// read precision matrix into eigen
		int n = Sigma.n_rows;
		Map<MatrixXd> Q(Sigma.memptr(), n, n);
		
		// read mean component into eigen
		Map<MatrixXd> _y(y.memptr(), n, 1);
		
		// factor precision matrix and compute r.v. mean; store upper cholesky
		LLT<MatrixXd, Upper> llt(Q);
		VectorXd mu = llt.solve(_y);
		
		// generate independent normal variates
		GetRNGstate();
		mat _z = randn(n, nSamples);
		PutRNGstate();
		
		// map variates into eigen
		Map<MatrixXd> z(_z.memptr(), n, nSamples);
		
		// compute r.v. with precision matrix Q and mean mu
		MatrixXd _x = llt.matrixU().solve(z);
		_x.colwise() += mu;
		
		// convert to arma vec
		mat x = mat(_x.data(), Sigma.n_rows, nSamples);

		return x;
	} else {
		return zeros(1);
	}
}


mat mcstat2::mvrnorm_postKron(vec & _y, mat & _A, mat & _B, int nSamples,
							  bool precision) {
	using Eigen::MatrixXd;
	using Eigen::VectorXd;
	using Eigen::Map;
	using Eigen::LLT;
	using Eigen::Upper;
	
	if(precision) {
		
		// read precision matrices into eigen
		int nA = _A.n_rows;
		int nB = _B.n_rows;
		Map<MatrixXd> Qa(_A.memptr(), nA, nA);
		Map<MatrixXd> Qb(_B.memptr(), nB, nB);
		
		// read mean component into eigen as matrix
		Map<MatrixXd> Y(_y.memptr(), nB, nA);
		
		// factor precision matrix and compute r.v. mean; store upper choleskys
		LLT<MatrixXd, Upper> lltA(Qa);
		LLT<MatrixXd, Upper> lltB(Qb);
		MatrixXd b = lltB.solve(Y);
		MatrixXd muMat = lltA.solve(b.transpose()).transpose();
		VectorXd mu(Map<VectorXd>(muMat.data(), nA * nB));
		
		// generate independent normal variates
		GetRNGstate();
		mat _z = randn(nB, nA * nSamples);
		PutRNGstate();
		
		// map variates into eigen
		Map<MatrixXd> z(_z.memptr(), nB, nA * nSamples);
		
		// compute r.v. with precision matrix A x B (but in matrix form)
		MatrixXd _x(nB, nA * nSamples);
		for(int i=0; i<nSamples; i++) {
			MatrixXd b = lltB.matrixU().solve(z.block(0, i*nA, nB, nA));
			_x.block(0, i*nA, nB, nA) =
				lltA.matrixU().solve(b.transpose()).transpose();
		}
		
		// reshape to column vectors, add mean mu
		Map<MatrixXd> _xV(_x.data(), nA * nB, nSamples);
		_xV.colwise() += mu;
		
		// convert to arma vec
		mat x = mat(_xV.data(), nA * nB, nSamples);
		
		return x;
	} else {
		return zeros(1);
	}
}




//
// Rcpp exports
//


// [[Rcpp::export]]
arma::mat r_mc2_rinvwishart(arma::mat V, double n) {
	return mcstat2::rinvwishart(V, n);
}


// [[Rcpp::export]]
arma::mat r_mvrnorm_postKron(arma::vec y, arma::mat A, arma::mat B,
							int nSamples, bool precision) {
	return mcstat2::mvrnorm_postKron(y, A, B, nSamples, precision);
}


// [[Rcpp::export]]
arma::mat r_mvrnorm_post(arma::vec y, arma::mat Sigma, int nSamples,
						 bool precision) {
	return mcstat2::mvrnorm_post(y, Sigma, nSamples, precision);
}


// [[Rcpp::export]]
arma::vec r_qintnorm(arma::vec breaks, double mu, double sigma) {
	return mcstat2::qintnorm(breaks, mu, sigma);
}
