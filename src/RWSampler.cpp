#include "RWSampler.h"
#include "transformations.h"



double mcstat2::RWSampler::getAcceptanceRate() { return accept; }

double mcstat2::RWSampler::getSd() { return sd; }

void mcstat2::RWSampler::adapt(double adaptScale, double targetRate) {
	sd *= exp( adaptScale * (accept - targetRate) );
}

Rcpp::List mcstat2::RWSampler::toList() {
	return Rcpp::List::create(
						Rcpp::_["accept"] = accept,
						Rcpp::_["sd"] = sd
						);
}

double mcstat2::RWSampler::sample(double x0) {

	double logR = 0;
	double x = x0;

	switch(propType) {
		case NORMAL:
			x += R::rnorm(0, sd);
			logR = logR_posterior(x, x0);
			break;

		case LOG:
			x = mcstat2::logProposal(x0, sd);
			logR = logR_posterior(x, x0) +
				mcstat2::loglogJacobian(x0) - mcstat2::loglogJacobian(x);
			break;

		case LOGIT:
			x = mcstat2::logitProposal(x0, L, U, sd);
			logR = logR_posterior(x, x0) +
				mcstat2::loglogitJacobian(x0) - mcstat2::loglogitJacobian(x);
			break;
	}

	bool accepted = log(R::runif(0,1)) <= std::min(logR, 0.0);

	if(accepted) {
		update();
		current = x;
	} else {
		x = x0;
	}

	accept += ((accepted ? 1.0 : 0.0) - accept) / (double) (++nSamples);

	adapt( C / sqrt( (double) (nSamples)), alpha);

	return x;
}

arma::vec mcstat2::RWSampler::sample() {
	double x = mcstat2::RWSampler::sample(current);
	arma::vec s = {x};
	return s;
}


//
// Rcpp exports
//

// [[Rcpp::depends(RcppArmadillo)]]

class GammaRWSampler : public mcstat2::RWSampler {

	// Construct a RW Sampler for approximating integrals of a Gamma distribution

	private:

		// gamma distribution parameters and integration constant
		double alpha, beta, logc, am1;

		double ll(double x) {
			return logc + am1 * std::log(x) - beta * x;
		}

		double logR_posterior(double x, double x0) {
			return ll(x) - ll(x0);
		}

	public:

		GammaRWSampler(double a, double b, double x0, double sd) :
			RWSampler(sd, x0, 1, .44, "x", LOG) {
			alpha = a;
			beta = b;
			am1 = alpha - 1;
			logc = a * std::log(b) - std::lgamma(a);
		}

};

// [[Rcpp::export]]
Rcpp::List test_rw_sampler(arma::vec params, double init, int nSamples) {

	GammaRWSampler cs = GammaRWSampler(params[0], params[1], init, params[2]);

	mcstat2::GibbsSampler sampler = mcstat2::GibbsSampler();
	sampler.addSampler(cs);
	sampler.run(nSamples);

	Rcpp::Rcout << "Acceptance rate: " << cs.getAcceptanceRate() << std::endl <<
		"sd: " << cs.getSd() << std::endl;

	return sampler.getSamples();
}
