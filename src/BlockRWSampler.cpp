#include "BlockRWSampler.h"
#include "transformations.h"


// helper function to create vector with n repeated double values
std::vector<double> repeatVals(double v, int n) {
	std::vector<double> r(n, v);
	return r;
}

mcstat2::BlockRWSampler::BlockRWSampler(
	std::vector<std::string> names,
	std::vector<ProposalType> t_propTypes,
	std::vector<double> t_sd,
	std::vector<double> init_vals,
	std::vector<double> t_C,
	double t_alpha,
	std::vector<double> t_L,
	std::vector<double> t_U
) : BlockSampler(mcstat2::repeatTypes(REAL, names.size()), names) {

	propTypes = t_propTypes;
	sd = t_sd;
	current = init_vals;
	C = t_C;
	alpha = t_alpha;

	// derived values
	dim = names.size();

	// "expand" default settings for logit transformations
	if(t_L.size()==0) { L = repeatVals(0, dim); }
	else { L = t_L; }

	if(t_U.size()==0) { U = repeatVals(1, dim); }
	else { U = t_U; }
}

double mcstat2::BlockRWSampler::getAcceptanceRate() { return accept; }

double mcstat2::BlockRWSampler::getSd(int i) { return sd.at(i); }

arma::vec mcstat2::BlockRWSampler::returnSamples(int i) {
	arma::vec s = {current[i]};
	return s;
}

int mcstat2::BlockRWSampler::getSize(int i) { return 1; }

void mcstat2::BlockRWSampler::drawSample() {

	// build proposal
	std::vector<double> x(current);
	for(int i=0; i<dim; i++) {
		switch (propTypes[i]) {
			case NORMAL:
				x[i] += R::rnorm(0, sd[i]);
				break;
			case LOG:
				x[i] = mcstat2::logProposal(current[i], sd[i]);
				break;
			case LOGIT:
				x[i] = mcstat2::logitProposal(current[i], L[i], U[i], sd[i]);
				break;
		}
	}

	// accept/reject
	
	double logR = logAcceptProb(x, current);
	bool accepted = log(R::runif(0,1)) <= std::min(logR, 0.0);

	// adapt, then update

	accept += ((accepted ? 1.0 : 0.0) - accept) / (double) (++nSamples);
	adapt(x, current);

	if(accepted) {
		update();
		current = x;
	}

}

double mcstat2::BlockRWSampler::logAcceptProb(
	const std::vector<double>& x, const std::vector<double>& x0) {

	// compute likelihood ratio
	double logR = logR_posterior(x, x0);

	// add jacobians
	for(int i=0; i<dim; i++) {
		switch (propTypes[i]) {
			case NORMAL:
				break;
			case LOG:
				logR += mcstat2::loglogJacobian(x0[i]) -
					mcstat2::loglogJacobian(x[i]);
				break;
			case LOGIT:
				logR += mcstat2::loglogitJacobian(x0[i]) -
					mcstat2::loglogitJacobian(x[i]);
				break;
		}
	}

	return logR < 0 ? logR : 0;
}

void mcstat2::BlockRWSampler::adapt(
	const std::vector<double>& x,
	const std::vector<double>& x0) {

	std::vector<double> xpartial(x0);

	double commonScale = (double) (nSamples);

	for(int i=0; i<dim; i++) {
		xpartial[i] = x[i];

		sd[i] *= std::exp(
			C[i] / commonScale * (std::exp(logAcceptProb(xpartial, x0)) - alpha)
		);

		xpartial[i] = x0[i];
	}

}


//
// Rcpp exports
//

// [[Rcpp::depends(RcppArmadillo)]]

class JointGaussianSampler : public mcstat2::BlockRWSampler {

	// Construct a Block RW sampler for approximating integrals of a multivariate
	// Gaussian distribution, where two dimensions are highly correlated, and the
	// third is independent.

	public:

		JointGaussianSampler(
			double rho, double x3var,
			std::vector<double> t_sd,
			std::vector<double> init_vals,
			std::vector<double> t_C) : BlockRWSampler({"x1", "x2", "x3"},
				{NORMAL, NORMAL, NORMAL}, t_sd, init_vals, t_C, .23) {

			// initialize covariance matrix
			sigma = arma::eye<arma::mat>(3,3);
			sigma.at(0,1) = rho;
			sigma.at(1,0) = rho;
			sigma.at(2,2) = x3var;

			// compute integration constant
			double val, sign;
			log_det(val, sign, sigma);
			lCst = -.5 * (val + 3 * std::log(6.2831853072));
		}

	private:

		double lCst;
		arma::mat sigma;

		double ll(const std::vector<double>& t_x) {
			arma::vec x = arma::vec(t_x);
			arma::vec y = arma::solve(sigma, x);
			return lCst - .5 * dot(x,y);
		}

		double logR_posterior(
			const std::vector<double>& x,
			const std::vector<double>& x0
		) { return ll(x) - ll(x0); }

};


// [[Rcpp::export]]
Rcpp::List test_blockrw_sampler(std::vector<double> params,
	std::vector<double> sd, std::vector<double> init, std::vector<double> C,
	int nSamples) {

	JointGaussianSampler gs =
		JointGaussianSampler(params[0], params[1], sd, init, C);

		mcstat2::GibbsSampler sampler = mcstat2::GibbsSampler();
		sampler.addSampler(gs);
		sampler.run(nSamples);

		Rcpp::Rcout << "Acceptance rate: " << gs.getAcceptanceRate() << std::endl
			<< "sds: " << gs.getSd(0) << ", " << gs.getSd(1) << ", " <<
			gs.getSd(2) << std::endl;

	return sampler.getSamples();
}
