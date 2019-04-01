#ifndef _BLOCKRWSAMPLER_H
#define _BLOCKRWSAMPLER_H

#include <RcppArmadillo.h>
#include "GibbsSampler.h"

namespace mcstat2 {

	//
	// multivariate adaptive random walk sampler
	//

	class BlockRWSampler : public BlockSampler {
		public:

			enum ProposalType { NORMAL, LOG, LOGIT };

			BlockRWSampler(
				std::vector<std::string> names,
				std::vector<ProposalType> t_propTypes,
				std::vector<double> t_sd,
				std::vector<double> init_vals,
				std::vector<double> t_C,
				double t_alpha,
				std::vector<double> t_L = std::vector<double>(),
				std::vector<double> t_U = std::vector<double>()
			);

			double getAcceptanceRate();
			double getSd(int i);

			// implement parent class methods
			void drawSample();
			vec returnSamples(int i);
			int getSize(int i);

		private:

			int nSamples, dim;
			double accept, alpha;
			std::vector<double> sd, L, U, C, current;
			std::vector<ProposalType> propTypes;

			void adapt(const std::vector<double>& x,
							   const std::vector<double>& x0);

			// compute log-acceptance probability for the pair (x,x0)
			double logAcceptProb(const std::vector<double>& x,
													 const std::vector<double>& x0);

		protected:

			// return ll(x) + pi(x) - ( ll(x0) + pi(x0) ), i.e., the log Metropolis
			// ratio for the proposed value x relative to the current value x0;
			// jacobians for sampling transformations are added within sampling
			// function
			virtual double logR_posterior(const std::vector<double>& x,
																		const std::vector<double>& x0) = 0;

			// run this code if the proposal is accepted
			virtual void update() {};

	};

}

#endif
