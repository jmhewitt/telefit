#ifndef _RWSAMPLER_H
#define _RWSAMPLER_H

#include <RcppArmadillo.h>
#include "GibbsSampler.h"

namespace mcstat2 {

	// univariate adaptive random walk sampler
	class RWSampler : public Sampler {

	public:
		enum ProposalType { NORMAL, LOG, LOGIT };
		// implement Sampler class' sample method
		vec sample();
		// sample a new parameter value given the current value x0
		double sample(double x0);
		// auto-tune proposal s.d.
		void adapt(double adaptScale, double targetRate);
		// return list with acceptance rate and proposal s.d.
		Rcpp::List toList();

		RWSampler(double t_sd, double t_current, double t_C, double t_alpha,
							std::string nom, ProposalType t_type, double t_L=0,
							double t_U=1) : Sampler(REAL, nom) {
			nSamples = 0;
			accept = 0;
			sd = t_sd;
			current = t_current;
			C = t_C;
			alpha = t_alpha;
			L = t_L;
			U = t_U;
			propType = t_type;
		}

		double getAcceptanceRate();
		double getSd();

private:

	int nSamples;
	double accept, sd, L, U, C, alpha, current;
	ProposalType propType;

protected:

	// return likelihood and prior, no jacobian; gets a copy of
	// proposed value x and current value x0
	virtual double logR_posterior(double x, double x0) = 0;

	// run this code if the proposal is accepted
	virtual void update() {};

};

}

#endif
