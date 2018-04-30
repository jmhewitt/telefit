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
	mat X, // each row has local covariates for one location + timepoint (Nntxp)
	    Z, // each column has remote covariates for one timepoint (kxnt)
	    d; // matrix containing interpoint distances (NxN)
	
	mat T,	// posterior samples
		beta,
		theta;
	vec sigmasq,
		sigmasqeps,
		rho;
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
		beta = samples.beta.row(i).t();
		theta = samples.theta.row(i).t();
		//T = mat(samples.T.row(i));
		//T.reshape(Tdim, Tdim);
		sigmasq = samples.sigmasq.at(i);
		sigmasqeps = samples.sigmasqeps.at(i);
		rho = samples.rho.at(i);
	}
};

struct Consts {
	int nSamples, // number of posterior draws available
		N,		  // number of spatial locations
		nt,	      // number of prediction timepoints
		k;		  // dimension of T
};

struct Config {
	Data data; Params params; Consts consts;
};

class YSampler : public mcstat2::Sampler {

private:
	
	Config *cfg;
	int it;
	mat H;
	vec zTheta;
	
public:
	
	YSampler(Config &_cfg) { name = "y"; type = VECTOR; cfg = &_cfg;
	 it = 0;
	 H = mat(cfg->consts.N, cfg->consts.N, fill::zeros);
	 zTheta = vec(cfg->consts.nt * cfg->consts.N, fill::zeros);
	}
	
	int getSize() { return cfg->consts.N * cfg->consts.nt; }
	
	vec sample() {
		
		// extract posterior parameter sample
		
		cfg->params.set(cfg->data, it++, cfg->consts.k);
		
		maternCov( H, cfg->data.d, cfg->params.sigmasq, cfg->params.rho,
				   cfg->params.nu,
				   cfg->params.sigmasq * cfg->params.sigmasqeps );
		
		// sample spatial noise
		mat w = mcstat2::mvrnorm(H, cfg->consts.nt, false);
		
		// compute posterior mean component
		
		for(int i=0; i<cfg->consts.nt; i++)
			zTheta.rows(i * cfg->consts.N, (i+1) * cfg->consts.N - 1) =
				mcstat2::dgeikmm(cfg->consts.N, cfg->data.Z.col(i).t(),
								 cfg->params.theta);
		
		// save MCMC output
	    return cfg->data.X * cfg->params.beta + vectorise(w);
	}
};


RcppExport SEXP _svcpredict (SEXP _samples, SEXP Xn, SEXP Zn, SEXP d, SEXP nu) {
	
	using namespace Rcpp;
	
	// extract model configuration
	
	List samples = as<List>(_samples);
	
	Config cfg = Config();
	
	cfg.data.X = as<mat>(Xn);
	cfg.data.Z = as<mat>(Zn);
	cfg.data.d = as<mat>(d);
	
	cfg.data.T = as<mat>(samples["T"]);
	cfg.data.beta = as<mat>(samples["beta"]);
	cfg.data.theta = as<mat>(samples["theta"]);
	cfg.data.sigmasq = as<vec>(samples["sigmasq"]);
	cfg.data.sigmasqeps = as<vec>(samples["sigmasqeps"]);
	cfg.data.rho = as<vec>(samples["rho"]);
	
	cfg.consts.nSamples = cfg.data.rho.size();
	cfg.consts.N = cfg.data.d.n_rows;
	cfg.consts.nt = cfg.data.Z.n_cols;
	cfg.consts.k = sqrt(cfg.data.T.n_cols);
	
	cfg.params.nu = as<double>(nu);
	
	// instantiate and run samplers
	
	YSampler ys = YSampler(cfg);
	
	mcstat2::GibbsSampler sampler = mcstat2::GibbsSampler();
	sampler.addSampler(ys);
	sampler.run(cfg.consts.nSamples);
	
	// return samples
	sampler.getSamples();
	return sampler.getSamples();
}
