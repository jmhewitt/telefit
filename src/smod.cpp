#include <RcppArmadillo.h>
#include "mcstat.h"
#include "covs.h"


using namespace Rcpp;
using namespace arma;

//
// spatial model without teleconnection effects
//

class SModel {
	
	
	private:
	
	// "constants"
	double one;
	mat Int;
	
	// dimensions
	int p, ns, nt, N, ns_minusOne;
	
	// data
	mat X, Dy;
	vec Y;
	
	// priors
	double nu_y, ay, by, ary, bry, aeps, beps, sigmasq_y_shape;
	mat Lambda;
	mat LambdaInv;
	
	// tuning parameters
	double rho_y_sd, eps_sd, rho_y_accept, eps_accept;
	
	// nested class holds model parameters and dependent structures
	class Params {
		
		public:
		
		// parameters
		vec beta;
		double sigmasq_y, rho_y, sigmasq_eps;
		
		// derived terms
		vec resid;
		mat Sigma, SigmaInv;
		double logdet_Sigma;
		
		// constructor zero-initializes parameters
		Params(int p=1, int ns=1, int nt=1, int N=1) {
			beta = vec(p, fill::zeros);
			resid = vec(N, fill::zeros);
			
			Sigma = mat(ns, ns, fill::zeros);
			SigmaInv = mat(ns, ns, fill::zeros);
		}
	};
	
	// sampler state
	Params gibbsCur;
	Params gibbsProposed;
	
	
	//
	// sampler functions
	//
	
	double logR_sigmasq_eps() {
		
		//
		// minimum proposed state update for computing ratios
		//
		
		maternCov(gibbsProposed.Sigma, Dy, gibbsCur.sigmasq_y,
				  gibbsCur.rho_y, nu_y,
				  gibbsCur.sigmasq_y * gibbsProposed.sigmasq_eps);
		gibbsProposed.SigmaInv = inv_sympd(gibbsProposed.Sigma);
		log_det(gibbsProposed.logdet_Sigma, one, gibbsProposed.Sigma);
		
		
		//
		// compute and return log MH ratio
		//
		
		vec logexp = gibbsCur.resid.t() * mcstat::dgemkmm(Int,
				gibbsProposed.SigmaInv - gibbsCur.SigmaInv, gibbsCur.resid);
		
		return ( nt * ( gibbsCur.logdet_Sigma - gibbsProposed.logdet_Sigma ) -
				logexp.at(0) ) / 2.0;
	}
	void sigmasq_eps_update() {
		gibbsCur.sigmasq_eps = gibbsProposed.sigmasq_eps;
		gibbsCur.Sigma = gibbsProposed.Sigma;
		gibbsCur.SigmaInv = gibbsProposed.SigmaInv;
		gibbsCur.logdet_Sigma = gibbsProposed.logdet_Sigma;
	}

	
	double logR_rho_y() {
		
		//
		// minimum proposed state update for computing ratios
		//
		
		maternCov(gibbsProposed.Sigma, Dy, gibbsCur.sigmasq_y,
				  gibbsProposed.rho_y, nu_y,
				  gibbsCur.sigmasq_y * gibbsCur.sigmasq_eps);
		gibbsProposed.SigmaInv = inv_sympd(gibbsProposed.Sigma);
		log_det(gibbsProposed.logdet_Sigma, one, gibbsProposed.Sigma);
		
		
		//
		// compute and return log MH ratio
		//
		
		vec logexp = gibbsCur.resid.t() * mcstat::dgemkmm(Int,
				gibbsProposed.SigmaInv - gibbsCur.SigmaInv, gibbsCur.resid);
		
		return ( nt * ( gibbsCur.logdet_Sigma - gibbsProposed.logdet_Sigma ) -
				logexp.at(0) ) / 2.0;
	}
	void rho_y_update() {
		gibbsCur.rho_y = gibbsProposed.rho_y;
		gibbsCur.Sigma = gibbsProposed.Sigma;
		gibbsCur.SigmaInv = gibbsProposed.SigmaInv;
		gibbsCur.logdet_Sigma = gibbsProposed.logdet_Sigma;
	}
	
	
	
	
	
	
	public:
	
	// sampler output
	mat beta_samples;
	vec sigmasq_y_samples, rho_y_samples, sigmasq_eps_samples, ll_samples;
	
	
	SModel(int _p, int _ns, int _nt) {
		p = _p;		// number of mean parameters
		ns = _ns;	// number of local locations
		nt = _nt;	// number of time points
		
		// compute additional dimension-based constants
		N = ns * nt;
		ns_minusOne = ns - 1;
		
		// initialize state
		gibbsCur = Params(p, ns, nt, N);
		gibbsProposed = Params(p, ns, nt, N);
		
		// initialize "constants"
		one = 1.0;
		Int = eye(nt, nt);
	}
	
	
	void setData(const mat & _X, const vec & _Y, const mat & _Dy) {
		X = _X;		// spatio-temporal design matrix (ns*nt x p)
		Y = _Y;		// spatio-temporal response (ns*nt x 1)
		Dy = _Dy;	// local distance matrix
	}
	
	
	void setPriors(double _nu_y, double _ay, double _by, double _ary,
				   double _bry, double _aeps, double _beps,
				   const mat & _Lambda) {
		nu_y = _nu_y;		// fixed smoothness parameters for covariance matrices
		ay = _ay;			// prior parameters for sigmasq_y
		by = _by;
		ary = _ary;			// prior parameters for rho_y
		bry = _bry;
		aeps = _aeps;		// prior parameters for nugget
		beps = _beps;
		Lambda = _Lambda;	// prior covariance for mean parameters
		
		
		//
		// compute constants associated with priors
		//
		
		// conjugate sampling beta
		LambdaInv = inv_sympd(Lambda);
		
		// conjugate sampling sigmasq_y
		sigmasq_y_shape = ((double) N / 2.0) + ay;
	}
	
	
	void tuneSampler(double _rho_y_sd, double _eps_sd) {
		rho_y_sd = _rho_y_sd;			// set variance for logit(rho_y)
		eps_sd = _eps_sd;				// set variance for log(nugget)
	}
	
	
	
	void fit(int maxIt, bool returnll, Function errDump, double C, double RWrate) {
		
		// use some Rsugar to wake up R's random number generator on crossbow
		rgamma(1, 2.0, 1.0);
		
		// initialize output
		beta_samples = mat(maxIt, p, fill::zeros);
		sigmasq_y_samples = vec(maxIt, fill::zeros);
		rho_y_samples = vec(maxIt, fill::zeros);
		sigmasq_eps_samples = vec(maxIt, fill::zeros);
		if(returnll)
			ll_samples = vec(maxIt, fill::zeros);
		
		// initialize parameters
		gibbsCur.beta = mcstat::mvrnorm(Lambda);
		gibbsCur.sigmasq_y = 1.0 / R::rgamma(ay, 1.0 / by);
		gibbsCur.rho_y = R::runif(ary, bry);
		gibbsCur.sigmasq_eps = 1.0 / R::rgamma(ay, 1.0 / by);
		
		// initialize temp
		double logR = 0.0;
		double adaptSize;
		bool accept;
		mat CkSigmaInvX;
		vec logexp;
		
		// initialize dependencies
		maternCov(gibbsCur.Sigma, Dy, gibbsCur.sigmasq_y, gibbsCur.rho_y,
				  nu_y, gibbsCur.sigmasq_y * gibbsCur.sigmasq_eps);
		gibbsCur.SigmaInv = inv_sympd(gibbsCur.Sigma);

		// initialize tracking variables
		int checkpointIt = (int) maxIt * 0.1;
		std::clock_t start, lap, tmp_clock;
		start = std::clock();
		lap = std::clock();
		double duration, total, pctComplete, remaining;
		rho_y_accept = eps_accept = 1.0;
		
		// gibbs iterations
		char step;
		int it;
		try{
		for(it=0; it<maxIt; it++) {
			
			checkUserInterrupt();
			
			//
			// conjugate beta
			//
			
			step='B';
			
			// update posterior
			CkSigmaInvX = mcstat::dgemkmm(Int, gibbsCur.SigmaInv, X);
			mat postBetaCov = inv_sympd( LambdaInv + X.t() * CkSigmaInvX );
			vec postBetaMean = postBetaCov * CkSigmaInvX.t() * Y;
			
			// sample and save
			gibbsCur.beta = mcstat::mvrnorm(postBetaMean, postBetaCov);
			beta_samples.row(it) = gibbsCur.beta.t();
			
			// recompute residuals
			for(int i=0; i<nt; i++) {
				int rStart = i * ns;
				int rEnd = rStart + ns_minusOne;
				gibbsCur.resid.subvec(rStart, rEnd) = Y.rows(rStart, rEnd) -
					X.rows(rStart, rEnd) * gibbsCur.beta;
			}
			
			
			//
			// conjugate sigmasq_y
			//
			
			step='S';
			
			// update posterior
			vec qform = gibbsCur.resid.t() * mcstat::dgemkmm(Int ,
					gibbsCur.SigmaInv * gibbsCur.sigmasq_y, gibbsCur.resid);
			double sigmasq_y_scale = by + qform.at(0) / 2.0;
			
			// sample and save
			gibbsCur.sigmasq_y = 1.0 / R::rgamma(sigmasq_y_shape, 1.0 / sigmasq_y_scale);
			sigmasq_y_samples.at(it) = gibbsCur.sigmasq_y;
			
			// recompute local covariances
			maternCov(gibbsCur.Sigma, Dy, gibbsCur.sigmasq_y, gibbsCur.rho_y,
					  nu_y, gibbsCur.sigmasq_y * gibbsCur.sigmasq_eps);
			gibbsCur.SigmaInv = inv_sympd(gibbsCur.Sigma);
			log_det(gibbsCur.logdet_Sigma, one, gibbsCur.Sigma);
			
			
			//
			// RW rho_y
			//
			
			step='R';
			
			// MH sample and update
			gibbsProposed.rho_y = mcstat::logitProposal(gibbsCur.rho_y, ary,
														bry, rho_y_sd);
			logR = logR_rho_y() + mcstat::loglogitJacobian(gibbsCur.rho_y) -
				mcstat::loglogitJacobian(gibbsProposed.rho_y);
			accept = log(R::runif(0,1)) <= std::min(logR, 0.0);
			if(accept)
				rho_y_update();
			rho_y_accept += ((accept ? 1.0 : 0.0) - rho_y_accept) / (double) (it + 1);
			
			// save parameter
			rho_y_samples.at(it) = gibbsCur.rho_y;
			
			
			//
			// RW sigmasq_eps
			//
			
			 step = 'e';
			
			// MH sample and update
			gibbsProposed.sigmasq_eps = mcstat::logProposal(gibbsCur.sigmasq_eps,
															eps_sd);
			logR = logR_sigmasq_eps() +
				mcstat::loglogJacobian(gibbsCur.sigmasq_eps) -
				mcstat::loglogJacobian(gibbsProposed.sigmasq_eps) +
				mcstat::logdinvgamma_unscaled(gibbsProposed.sigmasq_eps, aeps, beps) -
				mcstat::logdinvgamma_unscaled(gibbsCur.sigmasq_eps, aeps, beps);
			accept = log(R::runif(0,1)) <= std::min(logR, 0.0);
			if(accept)
				sigmasq_eps_update();
			eps_accept += ((accept ? 1.0 : 0.0) - eps_accept) / (double) (it + 1);
			
			// save parameter
			sigmasq_eps_samples.at(it) = gibbsCur.sigmasq_eps;
			 
			
			
			//
			// ll
			//
			
			step = 'l';
			
			if(returnll) {
				logexp = gibbsCur.resid.t() * mcstat::dgemkmm(Int ,
							gibbsCur.SigmaInv, gibbsCur.resid);
				
				ll_samples.at(it) = ( - nt * gibbsCur.logdet_Sigma +
									    logexp.at(0) ) / 2.0;
			}
			
			
			
			//
			// adjust tuning (following a method presented in Andrieu and Thoms (2008))
			//
			
			adaptSize = C / sqrt( (double) (it + 1));
			rho_y_sd *= exp( adaptSize * (rho_y_accept - RWrate) );
			eps_sd *= exp( adaptSize * (eps_accept - RWrate) );
			
			
			
			//
			// checkpoint behaviors
			//
			
			if(it % checkpointIt == 0) {
				// compute time since last checkpoint
				tmp_clock = lap;
				lap = std::clock();
				duration = ( lap - tmp_clock ) / (double) CLOCKS_PER_SEC;
				
				// compute percent complete
				pctComplete = (double) it / (double) maxIt * 100.0;
				
				// compute remaining time
				total = ( lap - start ) / (double) CLOCKS_PER_SEC;
				remaining = ( 100.0 - pctComplete ) * ( total / pctComplete ) / 60.0;
				
				// output information
				Rcout << round(pctComplete) << "% complete" << " (" <<
					floor(duration * 10.0) / 10.0 << " seconds; " <<
					floor(remaining * 10.0) / 10.0 << " minutes remaining)" <<
					endl;
			}
			
		}
		} catch(...) {
			
			Rcout << "An error occured while sampling '" << step << "'"
				<< " in iteration " << it << endl;
			
			//
			// dump state
			//
			
			List current = List::create(
				_["beta"] = gibbsCur.beta,
				_["sigmasq_y"] = gibbsCur.sigmasq_y,
				_["rho_y"] = gibbsCur.rho_y,
				_["sigmasq_eps"] = gibbsCur.sigmasq_eps,
				_["Sigma"] = gibbsCur.Sigma
			);
			
			List proposed = List::create(
				_["beta"] = gibbsProposed.beta,
				_["sigmasq_y"] = gibbsProposed.sigmasq_y,
				_["rho_y"] = gibbsProposed.rho_y,
				_["sigmasq_eps"] = gibbsProposed.sigmasq_eps,
				_["Sigma"] = gibbsProposed.Sigma
			);
			
			List samples = List::create(
				_["beta"] = beta_samples,
				_["sigmasq_y"] = sigmasq_y_samples,
				_["rho_y"] = rho_y_samples,
				_["sigmasq_eps"] = sigmasq_eps_samples,
				_["ll"] = ll_samples
			);
			
			List tuning = List::create(
				_["rho_y_sd"] = rho_y_sd,
				_["eps_sd"] = eps_sd
			);
			
			List errStuff = List::create(
				_["samples"] = samples,
				_["proposed"] = proposed,
				_["current"] = current,
			    _["tuning"] = tuning
			);
			
			errDump(errStuff);
			
			Rcout << "  Sampler state saved" << endl;
		}
		
		// print final stats
		for(int i=0; i<80; i++)
			Rcout << "-";
		Rcout << endl << "Acceptance rates:" << endl;
		Rcout << "  sigmasq_eps: " << round(eps_accept * 100.0) << "%" << endl;
		Rcout << "  rho_y: " << round(rho_y_accept * 100.0) << "%" << endl;
	}
	
	
};


//
// R exports
//

RcppExport SEXP _sfit( SEXP p, SEXP ns, SEXP nt, SEXP X,
					    SEXP Y, SEXP Dy, SEXP nu_y,
					    SEXP ay, SEXP by, SEXP ary, SEXP bry,
					    SEXP aeps, SEXP beps, SEXP Lambda,
					    SEXP rho_y_sd, SEXP eps_sd,
					    SEXP maxIt, SEXP returnll,
					    SEXP errDump, SEXP C, SEXP RWrate ) {

	using namespace Rcpp;
	
	// instantiate sampler
	SModel smod = SModel(as<int>(p), as<int>(ns), as<int>(nt));
	
	
	//
	// configure sampler
	//
	
	smod.setData(as<mat>(X), as<vec>(Y), as<mat>(Dy));
	
	smod.setPriors(as<double>(nu_y), as<double>(ay), as<double>(by),
					as<double>(ary), as<double>(bry),
					as<double>(aeps), as<double>(beps), as<mat>(Lambda));
	
	smod.tuneSampler(as<double>(rho_y_sd), as<double>(eps_sd));
	
	
	// run sampler
	smod.fit(as<int>(maxIt), as<bool>(returnll), as<Function>(errDump),
			  as<double>(C), as<double>(RWrate));
	
	return List::create(
						_["beta"] = smod.beta_samples,
						_["sigmasq_y"] = smod.sigmasq_y_samples,
						_["rho_y"] = smod.rho_y_samples,
						_["sigmasq_eps"] = smod.sigmasq_eps_samples,
						_["ll"] = smod.ll_samples
						);

}
