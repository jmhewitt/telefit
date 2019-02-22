#include "GibbsSampler.h"

using namespace Rcpp;

mcstat2::MCMCCheckpoint::MCMCCheckpoint(int t_nSamples, int t_thin) {
	nSamples = t_nSamples;
	thin = t_thin;
	checkPointIt = (int) nSamples * 0.1;
	it = 0;
	start = std::clock();
	lap = start;
}

void mcstat2::MCMCCheckpoint::reset() {
	it = 0;
	start = std::clock();
	lap = start;
}

void mcstat2::MCMCCheckpoint::run() {
	if(it++ % checkPointIt == 0) {
		// compute time since last checkpoint
		std::clock_t tmp_clock = lap;
		lap = std::clock();
		double duration = ( lap - tmp_clock ) / (double) CLOCKS_PER_SEC;

		// compute percent complete
		double pctComplete = (double) it / (double) nSamples * 100.0;
		if(it==1) {
			pctComplete /= (double) thin;
		}

		// compute remaining time
		double total = (lap - start) / (double) CLOCKS_PER_SEC;
		double remaining = (100.0 - pctComplete) * (total / pctComplete) / 60.0;

		// output information
		Rcout << round(pctComplete) << "% complete" << " (" <<
		floor(duration * 10.0) / 10.0 << " seconds; " <<
		floor(remaining * 10.0) / 10.0 << " minutes remaining)" <<
		endl;
	}
}

void mcstat2::MCMCCheckpoint::finish() {
	// compute final time
	lap = std::clock();
	double duration = (lap - start) / (double) CLOCKS_PER_SEC;

	// compute samples per second
	double sampSec = it / duration;

	// output information
	Rcout << endl << "Total time (min): " << floor(duration / 60.0 * 10.0) / 10.0 <<
	endl << "Samples per second: " << floor(sampSec * 10.0) / 10.0 << endl;
}

void mcstat2::GibbsSampler::setThinning(int t_thin) {
	thin = t_thin < 1 ? 1 : t_thin;
}

void mcstat2::GibbsSampler::addSampler(BlockSampler &s) {
	samplers.push_back(&s);
}

void mcstat2::GibbsSampler::run(int nSamples) {

	GetRNGstate();

	// use some Rsugar to wake up R's random number generator on crossbow
	rgamma(1, 2.0, 1.0);

	// initialize trackers
	mcstat2::MCMCCheckpoint checkpoint = mcstat2::MCMCCheckpoint(nSamples, thin);

	// initialize samples
	for(auto block = samplers.begin(); block != samplers.end(); ++block) {
		int blockSize = (*block)->getNumSamplers();
		for(int i=0; i<blockSize; i++) {
			switch((*block)->getType(i)) {
				case BlockSampler::REAL:
					samples.push_back(mat(nSamples, 1, fill::zeros));
					break;
				case BlockSampler::VECTOR:
					samples.push_back(mat(nSamples, (*block)->getSize(i), fill::zeros));
					break;
			}
		}
	}

	int it;
	int totSamples = nSamples * thin;
	int saved = 0;
	std::string step;
	try{
		// reset timers
		checkpoint.reset();

		// gibbs iterations
		for(it=0; it < totSamples; it++) {

			checkUserInterrupt();

			// gibbs steps: iterate over sampling blocks
			auto sample = samples.begin();
			for(auto block = samplers.begin(); block != samplers.end(); ++block) {

				// prep for error handling: identify first step in block
				step = (*block)->getName(0);

				// sample block
				(*block)->drawSample();

				// extract all samples from block
				int blockSize = (*block)->getNumSamplers();
				for(int i=0; i<blockSize; i++) {

					// store ith block element
					vec s = (*block)->returnSamples(i);

					if(it % thin == 0)
						sample->row(saved) = s.t();

					// iterate storage pointer
					++sample;
				}

			}

			// update output index
			if(it % thin == 0)
				saved++;

			// checkpoint behaviors
			if(it % thin == 0)
				checkpoint.run();

		}

	} catch(...) {
		Rcout << "An error occured while sampling " << step <<
		" in iteration " << it << " for sample " << saved << endl;

		// TODO: dump state by cycling through samplers to get their state
	}


	// print final sampler stats
	Rcout << std::setfill('-') << std::setw(80) << "-" << endl;

	// timings
	checkpoint.finish();

	// print sampler stats
	for(auto block = samplers.begin(); block != samplers.end(); ++block) {
		int blockSize = (*block)->getNumSamplers();
		for(int i=0; i<blockSize; i++)
			(*block)->printStats(i);
	}


	PutRNGstate();
}

List mcstat2::GibbsSampler::getSamples() {

	// initialize output containers for parameter names and samples
	int n = 0;
	for(auto block = samplers.begin(); block != samplers.end(); ++block) {
		n += (*block)->getNumSamplers();
	}
	List res(n);
	CharacterVector parameterNames(n);

	// populate output containers
	int i=0;
	std::vector<mat>::iterator sample = samples.begin();
	for(auto block = samplers.begin(); block != samplers.end(); ++block) {
		int blockSize = (*block)->getNumSamplers();
		for(int j=0; j<blockSize; j++) {
			// prepare sample label
			parameterNames[i] = (*block)->getName(j);
			// export samples
			res[i++] = (*sample);
			// iterate samples
			++sample;
		}
	}

	// label samples
	res.names() = parameterNames;

	return res;
}

mcstat2::Sampler::SamplerType mcstat2::Sampler::getType() {
	return type;
}

std::string mcstat2::Sampler::getName() {
	return name;
}

int mcstat2::Sampler::getSize() {
	return type == REAL ? 1 : sample().size();
}

mcstat2::BlockSampler::SamplerType mcstat2::BlockSampler::getType(int i) {
	return types.at(i);
}

std::string mcstat2::BlockSampler::getName(int i) {
	return names.at(i);
}


void mcstat2::Sampler::drawSample() {
	tsample = sample();
}

arma::vec mcstat2::Sampler::returnSamples(int i) {
	return tsample;
}

void mcstat2::Sampler::printStats(int i) {
	printStats();
}

int mcstat2::Sampler::getSize(int i) {
	return getSize();
}

int mcstat2::BlockSampler::getNumSamplers() {
	if(types.size() == 0) refreshTypesAndNames();
	return types.size();
}

void mcstat2::Sampler::refreshTypesAndNames() {
	types.push_back(type);
	names.push_back(name);
}


//
// Rcpp exports
//

// [[Rcpp::depends(RcppArmadillo)]]

class CountSampler : public mcstat2::Sampler {
	private:
		int i;
	public:
		CountSampler(std::string nom, int val) {
			name = nom; type = REAL; i = val;
		}
		arma::vec sample() { arma::vec s = {(double) i++}; return s;  }
};

// [[Rcpp::export]]
List test_gibbs_sampler(arma::vec inits, int nSamples) {

	CountSampler cs1 = CountSampler("a", inits[0]);
	CountSampler cs2 = CountSampler("b", inits[1]);
	CountSampler cs3 = CountSampler("c", inits[2]);

	mcstat2::GibbsSampler sampler = mcstat2::GibbsSampler();
	sampler.addSampler(cs1);
	sampler.addSampler(cs2);
	sampler.addSampler(cs3);

	sampler.run(nSamples);

	return sampler.getSamples();
}
