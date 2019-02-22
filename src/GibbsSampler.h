#ifndef _GIBBSSAMPLER_H
#define _GIBBSSAMPLER_H

#include <RcppArmadillo.h>

namespace mcstat2 {

	using namespace Rcpp;
	using namespace arma;

	class MCMCCheckpoint {

	private:
		int it, thin, checkPointIt, nSamples;
		std::clock_t lap, start;

	public:

		MCMCCheckpoint(int nSamples, int thin);

		void reset();
		void run();
		void finish();
	};


	// abstract class for a sampler that can be used with mcstat::GibbsSampler
	class BlockSampler {
		public:
			enum SamplerType { REAL, VECTOR };

			BlockSampler() { };

			// return type and name of ith sampler in block
			SamplerType getType(int i);
			std::string getName(int i);

			// draw samples for entire block
			virtual void drawSample() = 0;

			// extract the output from the ith sampler in the block
			virtual vec returnSamples(int i) = 0;

			// returns 1 or the size of sampled vectors for the ith sampler in block;
			// by default this is done by running the sample function, but inherited
			// classes can overload this
			virtual int getSize(int i) = 0;

			// for printing acceptance rates and tuning parameters, etc. to Rcout
			virtual void printStats(int i) {};

			// returns the number of samplers in block
			int getNumSamplers();

		protected:

			// for initializing MCMC output containers
			std::vector<SamplerType> types;

			// for labeling the samples when sent back to R as a List object
			std::vector<std::string> names;

			// backwards compatibility for simpler creation of Sampler classes
			virtual void refreshTypesAndNames() {};
	};


	// abstract class for a sampler that can be used with mcstat::GibbsSampler
	class Sampler : public BlockSampler  {

	public:

		Sampler() { }

		SamplerType getType();
		std::string getName();

		/* sampler should:
			1. update all relevant data parameters and output lists
			2. autotune itself, if necessary
			3. return the value of the sampled parameter
		 */
		virtual vec sample() = 0;

		// for printing acceptance rates and tuning parameters, etc. to Rcout
		virtual void printStats() {};

		// returns 1 or the size of sampled vectors; by default this is done by
		// running the sample function, but inherited classes can overload this
		virtual int getSize();

	protected:

		// for initializing MCMC output containers
		SamplerType type;
		// for labeling the samples when sent back to R as a List object
		std::string name;

	private:

		// for saving sampler output
		arma::vec tsample;

		// implement sampling methods inherited from BlockSampler
		void drawSample();
		vec returnSamples(int i);
		void printStats(int i);
		int getSize(int i);
		void refreshTypesAndNames();

	};

	// runs a collection of mcstat::Sampler objects in a Gibbs sampling scheme
	class GibbsSampler {

	private:
		std::vector<BlockSampler*> samplers;
		std::vector<mat> samples;
		int thin;

	public:

		GibbsSampler() { thin = 1; };

		void addSampler(BlockSampler &);
		void setThinning(int);

		// loops over the samplers, calling their sample functions
		void run(int);

		List getSamples();
	};
}

#endif
