/*
	gibbs sampler to forecast using a fitted spatially varying coefficient model.
 */

#include <RcppArmadillo.h>
#include "GibbsSampler.h"
#include "numAlg.h"
#include "covs.h"
#include "distributions.h"

using namespace Rcpp;
using namespace arma;

struct Data {
	mat *X, // each row has local covariates for 1 location + timepoint (Nntxp)
	    *Z, // each column has remote covariates for one timepoint (kxnt)
	    *d; // matrix containing interpoint distances (NxN)
	
	mat *T,	// posterior samples
		*beta,
		*theta;
	vec *sigmasq,
		*sigmasqeps,
		*rho;
};

struct Params {
	vec beta;		   // fixed coefficients
	vec theta;		   // spatially varying coefficients
	//mat T; 		       // local covariance for spatially varying coefficients
	double sigmasq,	   // scale for spatial covariance function
		   sigmasqeps, // scale for spatial nugget function (via param. ex.)
		   rho,		   // range for spatial covariance function
		   nu;		   // spatial covariance smoothness
	
	void set(const Data &samples, int i, int Tdim) {
		beta = samples.beta->row(i).t();
		theta = samples.theta->row(i).t();
		//T = mat(samples.T.row(i));
		//T.reshape(Tdim, Tdim);
		sigmasq = samples.sigmasq->at(i);
		sigmasqeps = samples.sigmasqeps->at(i);
		rho = samples.rho->at(i);
	}
};

struct Consts {
	int nSamples, // number of posterior draws available
		N,		  // number of spatial locations
		nt,	      // number of prediction timepoints
		k;		  // dimension of T
};


class YSampler : public mcstat2::Sampler {

private:
	
	Data *data;
	Params *params;
	Consts *consts;
	int it;
	mat H;
	vec zTheta;
	
public:
	
	YSampler(Data &_data, Params &_params, Consts &_consts) {
		name = "y"; type = VECTOR;
		data = &_data;
		params = &_params;
		consts = &_consts;
	 it = 0;
	 H = mat(consts->N, consts->N, fill::zeros);
	 zTheta = vec(consts->nt * consts->N, fill::zeros);
	}
	
	int getSize() { return consts->N * consts->nt; }
	
	vec sample() {
		
		// extract posterior parameter sample
		params->set(*data, it++, consts->k);
		
		maternCov( H, *(data->d), params->sigmasq, params->rho,
				   params->nu,
				   params->sigmasq * params->sigmasqeps );
		
		// sample spatial noise
		mat w = mcstat2::mvrnorm(H, consts->nt, false);
		
		// compute posterior mean component
		
		for(int i=0; i<consts->nt; i++)
			zTheta.rows(i * consts->N, (i+1) * consts->N - 1) =
				mcstat2::dgeikmm(consts->N, data->Z->col(i).t(),
								 params->theta );
		
		// save MCMC output
		return *(data->X) * params->beta + zTheta + vectorise(w);
	}
};


RcppExport SEXP _svcpredict (SEXP _T, SEXP _beta, SEXP _theta, SEXP _sigmasq,
							 SEXP _sigmasqeps, SEXP _rho, SEXP _Xn, SEXP _Zn,
							 SEXP _d, SEXP nu) {
	
	using namespace Rcpp;
	
	// extract model configuration
	
	Data data = Data();
	Params params = Params();
	Consts consts = Consts();
	
	mat X = as<mat>(_Xn);
	mat Z = as<mat>(_Zn);
	mat d = as<mat>(_d);
	data.X = &X;
	data.Z = &Z;
	data.d = &d;
	
	mat T = as<mat>(_T);
	mat beta = as<mat>(_beta);
	mat theta = as<mat>(_theta);
	vec sigmasq = as<vec>(_sigmasq);
	vec sigmasqeps = as<vec>(_sigmasqeps);
	vec rho = as<vec>(_rho);
	
	data.T = &T;
	data.beta = &beta;
	data.theta = &theta;
	data.sigmasq = &sigmasq;
	data.sigmasqeps = &sigmasqeps;
	data.rho = &rho;
	
	consts.nSamples = data.rho->size();
	consts.N = data.d->n_rows;
	consts.nt = data.Z->n_cols;
	consts.k = sqrt(data.T->n_cols);
	
	params.nu = as<double>(nu);
	
	// instantiate and run samplers
	
	YSampler ys = YSampler(data, params, consts);
	
	mcstat2::GibbsSampler sampler = mcstat2::GibbsSampler();
	sampler.addSampler(ys);
	sampler.run(consts.nSamples);
	
	// return samples
	return sampler.getSamples();
}
