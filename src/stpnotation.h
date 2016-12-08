#include <RcppArmadillo.h>
#include "mcstat.h"

#ifndef _telefit_stpnotation_h
#define _telefit_stpnotation_h

using namespace Rcpp;
using namespace arma;

struct Priors {
	mcstat::mvnorm beta;
	mcstat::invgamma sigmasq_y, sigmasq_r, sigmasq_eps, sigmasq_r_eps;
	mcstat::uniform rho_y, rho_r;
	
	Priors() { }
	
	Priors(const vec &beta_mu, const mat &beta_Sigma, double sigmasq_y_shape,
		   double sigmasq_y_rate, double sigmasq_r_shape, double sigmasq_r_rate,
		   double sigmasq_eps_shape, double sigmasq_eps_rate, double rho_y_a,
		   double rho_y_b, double rho_r_a, double rho_r_b,
		   double sigmasq_r_eps_shape, double sigmasq_r_eps_rate) {
		
		beta = mcstat::mvnorm(beta_mu, beta_Sigma);
		sigmasq_y = mcstat::invgamma(sigmasq_y_shape, sigmasq_y_rate);
		sigmasq_r = mcstat::invgamma(sigmasq_r_shape, sigmasq_r_rate);
		sigmasq_r_eps = mcstat::invgamma(sigmasq_r_eps_shape, sigmasq_r_eps_rate);
		sigmasq_eps = mcstat::invgamma(sigmasq_eps_shape, sigmasq_eps_rate);
		rho_y = mcstat::uniform(rho_y_a, rho_y_b);
		rho_r = mcstat::uniform(rho_r_a, rho_r_b);
	}
};

struct Data {
	mat X, Z;
	vec Y;
	
	Data() { }
	
	// intended to be used for forecasting
	Data(const mat &_X, const mat &_Z) {
		X = _X;
		Z = _Z;
	}
	
	// intended to be used for fitting
	Data(const mat &_X, const mat &_Z, const vec &_Y) {
		X = _X;
		Z = _Z;
		Y = _Y;
	}
};

struct Constants {
	mat Dy, Dz_knots, Dz_to_knots;
	int p, ns, nr, nr_knots, nt;
	double smoothness_y, smoothness_r;
	bool localOnly;
	
	Constants() { }
	
	Constants(const mat &_Dy, const mat &_Dz_knots, const mat &_Dz_to_knots,
			  int _p, int _ns, int _nr, int _nr_knots, int _nt,
			  double _smoothness_y, double _smoothness_r, bool _localOnly) {
		Dy = _Dy;
		Dz_knots = _Dz_knots;
		Dz_to_knots = _Dz_to_knots;
		p = _p;
		ns = _ns;
		nr = _nr;
		nr_knots = _nr_knots;
		nt = _nt;
		smoothness_y = _smoothness_y;
		smoothness_r = _smoothness_r;
		localOnly = _localOnly;
	}
};



struct CompositionSamples {
	
	bool return_full_alpha, return_forecast, localOnly;
	
	mat alpha_knots, alpha;
	cube forecast, local, remote;
	
	CompositionSamples(int nSamples, const Constants &consts,
					   bool _return_full_alpha, int nt0=-1) {
		
		localOnly = consts.localOnly;
		
		return_forecast = nt0 > 0 ? true : false;
		return_full_alpha = _return_full_alpha;
		
		if(!localOnly) {
			alpha_knots = mat(nSamples, consts.ns * consts.nr_knots, fill::zeros);
			
			if(return_full_alpha)
				alpha = mat(nSamples, consts.ns * consts.nr, fill::zeros);
		}
		
		if(return_forecast) {
			forecast = cube(consts.ns, nt0, nSamples, fill::zeros);
			if(!localOnly) {
				local = cube(consts.ns, nt0, nSamples, fill::zeros);
				remote = cube(consts.ns, nt0, nSamples, fill::zeros);
			}
		}
	}
	
	List toSummarizedList();
	
};


struct Samples {
	
	mat beta;
	vec sigmasq_y, sigmasq_r, sigmasq_eps, rho_y, rho_r, ll, sigmasq_r_eps;
	
	Samples(Constants &consts, int nSamples) {
		beta = mat(nSamples, consts.p, fill::zeros);
		sigmasq_y = vec(nSamples, fill::zeros);
		sigmasq_eps = vec(nSamples, fill::zeros);
		rho_y = vec(nSamples, fill::zeros);
		ll = vec(nSamples, fill::zeros);
		if(!consts.localOnly) {
			sigmasq_r = vec(nSamples, fill::zeros);
			rho_r = vec(nSamples, fill::zeros);
			sigmasq_r_eps = vec(nSamples, fill::zeros);
		}
	}
	
	Samples(const mat &_beta, const vec &_sigmasq_y, const vec &_sigmasq_r,
			const vec &_sigmasq_eps, const vec &_rho_y, const vec &_rho_r,
			const vec &_ll, const vec &_sigmasq_r_eps) {
		beta = _beta;
		sigmasq_y = _sigmasq_y;
		sigmasq_r = _sigmasq_r;
		sigmasq_eps = _sigmasq_eps;
		rho_y = _rho_y;
		rho_r = _rho_r;
		ll = _ll;
		sigmasq_r_eps = _sigmasq_r_eps;
	}
	
	List toList() {
		return List::create(
			_["beta"] = beta,
			_["sigmasq_y"] = sigmasq_y,
			_["sigmasq_r"] = sigmasq_r,
			_["sigmasq_r_eps"] = sigmasq_r_eps,
			_["sigmasq_eps"] = sigmasq_eps,
			_["rho_y"] = rho_y,
			_["rho_r"] = rho_r,
			_["ll"] = ll
		);
	}
	
};


#endif
