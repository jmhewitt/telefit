#include <RcppArmadillo.h>
#include "mcstat.h"
#include "covs.h"


using namespace Rcpp;
using namespace arma;

//
// spatial teleconnection model
//

class STModel {
	
	
	private:
	
	// "constants"
	double one;
	
	// dimensions
	int p, r, ns, nt, N, ns_minusOne;
	
	// data
	mat X, Z, Dy, Dz;
	vec Y;
	
	// priors
	double nu_y, nu_r, ay, by, ar, br, ary, bry, arr, brr, aeps, beps,
		   sigmasq_y_shape;
	mat Lambda;
	mat LambdaInv;
	
	// tuning parameters
	double rho_y_sd, rho_r_sd, eps_sd, sigmasq_r_sd, rho_y_accept, rho_r_accept,
		   eps_accept, sigmasq_r_accept;
	
	// nested class holds model parameters and dependent structures
	class Params {
		
		public:
		
		// parameters
		vec beta;
		double sigmasq_y, rho_y, rho_r, sigmasq_r, sigmasq_eps;
		
		// derived terms
		vec resid;
		mat R, C, CInv, Sigma, SigmaInv;
		double logdet_Sigma, logdet_C;
		
		// constructor zero-initializes parameters
		Params(int p=1, int r=1, int ns=1, int nt=1, int N=1) {
			beta = vec(p, fill::zeros);
			resid = vec(N, fill::zeros);
			
			R = mat(r, r, fill::zeros);
			CInv = mat(r, r, fill::zeros);
			C = mat(r, r, fill::zeros);
			
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
		
		vec logexp = gibbsCur.resid.t() * mcstat::dgemkmm(gibbsCur.C,
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
	
	double logR_sigmasq_r() {
		//
		// minimum proposed state update for computing ratios
		//
		maternCov(gibbsProposed.R, Dz, gibbsProposed.sigmasq_r, gibbsCur.rho_r,
				  nu_r, 0);
		//poweredExpCov(gibbsProposed.R, Dz, gibbsProposed.sigmasq_r, gibbsCur.rho_r,
		//		  nu_r, 0);
		gibbsProposed.CInv = Z.t() * gibbsProposed.R * Z;
		gibbsProposed.CInv.diag() += 1;
		gibbsProposed.C = inv_sympd(gibbsProposed.CInv);
		log_det(gibbsProposed.logdet_C, one, gibbsProposed.C);
		
		
		//
		// compute and return log MH ratio
		//
		
		vec logexp = gibbsCur.resid.t() *
			mcstat::dgemkmm(gibbsProposed.C - gibbsCur.C, gibbsCur.SigmaInv,
							gibbsCur.resid);
		
		return ( ns * ( gibbsProposed.logdet_C - gibbsCur.logdet_C ) -
				logexp.at(0) ) / 2.0;
	}
	void sigmasq_r_update() {
		gibbsCur.sigmasq_r = gibbsProposed.sigmasq_r;
		gibbsCur.C = gibbsProposed.C;
		gibbsCur.CInv = gibbsProposed.CInv;
		gibbsCur.logdet_C = gibbsProposed.logdet_C;
	}
	
	double logR_rho_r() {
		//
		// minimum proposed state update for computing ratios
		//
		maternCov(gibbsProposed.R, Dz, gibbsCur.sigmasq_r, gibbsProposed.rho_r,
				  nu_r, 0);
		//poweredExpCov(gibbsProposed.R, Dz, gibbsCur.sigmasq_r, gibbsProposed.rho_r,
		//		  nu_r, 0);
		gibbsProposed.CInv = Z.t() * gibbsProposed.R * Z;
		gibbsProposed.CInv.diag() += 1;
		gibbsProposed.C = inv_sympd(gibbsProposed.CInv);
		log_det(gibbsProposed.logdet_C, one, gibbsProposed.C);
		
		
		//
		// compute and return log MH ratio
		//
		
		vec logexp = gibbsCur.resid.t() *
			mcstat::dgemkmm(gibbsProposed.C - gibbsCur.C, gibbsCur.SigmaInv,
							gibbsCur.resid);
		
		return ( ns * ( gibbsProposed.logdet_C - gibbsCur.logdet_C ) -
				 logexp.at(0) ) / 2.0;
	}
	void rho_r_update() {
		gibbsCur.rho_r = gibbsProposed.rho_r;
		gibbsCur.C = gibbsProposed.C;
		gibbsCur.CInv = gibbsProposed.CInv;
		gibbsCur.logdet_C = gibbsProposed.logdet_C;
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
		
		vec logexp = gibbsCur.resid.t() * mcstat::dgemkmm(gibbsCur.C,
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
	vec sigmasq_y_samples, rho_y_samples, rho_r_samples, sigmasq_r_samples,
		sigmasq_eps_samples, ll_samples;
	
	// secondary sampler output
	mat alpha_samples;
	
	
	STModel(int _p, int _r, int _ns, int _nt) {
		p = _p;		// number of mean parameters
		r = _r;		// number of remote locations
		ns = _ns;	// number of local locations
		nt = _nt;	// number of time points
		
		// compute additional dimension-based constants
		N = ns * nt;
		ns_minusOne = ns - 1;
		
		// initialize state
		gibbsCur = Params(p, r, ns, nt, N);
		gibbsProposed = Params(p, r, ns, nt, N);
		
		// initialize "constants"
		one = 1.0;
	}
	
	
	void setData(const mat & _X, const mat & _Z, const vec & _Y,
				 const mat & _Dy, const mat & _Dz) {
		X = _X;		// spatio-temporal design matrix (ns*nt x p)
		Z = _Z;		// spatio-temporal remote covariate matrix (r x nt)
		Y = _Y;		// spatio-temporal response (ns*nt x 1)
		Dy = _Dy;	// local distance matrix
		Dz = _Dz;	// remote distance matrix
	}
	
	
	void setPriors(double _nu_y, double _nu_r, double _ay, double _by,
				   double _ar, double _br, double _ary, double _bry,
				   double _arr, double _brr, double _aeps, double _beps,
				   const mat & _Lambda) {
		nu_y = _nu_y;		// fixed smoothness parameters for covariance matrices
		nu_r = _nu_r;
		ay = _ay;			// prior parameters for sigmasq_y
		by = _by;
		ar = _ar;			// prior parameters for sigmasq_r
		br = _br;
		ary = _ary;			// prior parameters for rho_y
		bry = _bry;
		arr = _arr;			// prior parameters for rho_r
		brr = _brr;
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
	
	
	void tuneSampler(double _rho_y_sd, double _rho_r_sd, double _eps_sd,
					 double _sigmasq_r_sd) {
		rho_y_sd = _rho_y_sd;			// set variance for logit(rho_y)
		rho_r_sd = _rho_r_sd;			// set variance for logit(rho_r)
		eps_sd = _eps_sd;				// set variance for log(nugget)
		sigmasq_r_sd = _sigmasq_r_sd;	// set variance for log(sigmasq_r)
	}
	
	
	
	void fit(int maxIt, bool returnll, Function errDump, double C, double RWrate) {
		
		// use some Rsugar to wake up R's random number generator on crossbow
		rgamma(1, 2.0, 1.0);
		
		// initialize output
		beta_samples = mat(maxIt, p, fill::zeros);
		sigmasq_y_samples = vec(maxIt, fill::zeros);
		rho_y_samples = vec(maxIt, fill::zeros);
		rho_r_samples = vec(maxIt, fill::zeros);
		sigmasq_r_samples = vec(maxIt, fill::zeros);
		sigmasq_eps_samples = vec(maxIt, fill::zeros);
		if(returnll)
			ll_samples = vec(maxIt, fill::zeros);
		
		// initialize parameters
		gibbsCur.beta = mcstat::mvrnorm(Lambda);
		gibbsCur.sigmasq_y = 1.0 / R::rgamma(ay, 1.0 / by);
		gibbsCur.rho_y = R::runif(ary, bry);
		gibbsCur.rho_r = R::runif(arr, brr);
		gibbsCur.sigmasq_r = 1.0 / R::rgamma(ar, 1.0 / br);
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
		maternCov(gibbsCur.R, Dz, gibbsCur.sigmasq_r, gibbsCur.rho_r, nu_r, 0);
		//poweredExpCov(gibbsCur.R, Dz, gibbsCur.sigmasq_r, gibbsCur.rho_r, nu_r, 0);
		gibbsCur.CInv = Z.t() * gibbsCur.R * Z;
		gibbsCur.CInv.diag() += 1;
		gibbsCur.C = inv_sympd(gibbsCur.CInv);
		log_det(gibbsCur.logdet_C, one, gibbsCur.C);

		// initialize tracking variables
		int checkpointIt = (int) maxIt * 0.1;
		std::clock_t start, lap, tmp_clock;
		start = std::clock();
		lap = std::clock();
		double duration, total, pctComplete, remaining;
		rho_y_accept = rho_r_accept = sigmasq_r_accept = eps_accept = 1.0;
		
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
			CkSigmaInvX = mcstat::dgemkmm(gibbsCur.C, gibbsCur.SigmaInv, X);
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
			vec qform = gibbsCur.resid.t() * mcstat::dgemkmm(gibbsCur.C ,
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
			// RW rho_r
			//
			
			step='r';
			
			// MH sample and update
			gibbsProposed.rho_r = mcstat::logitProposal(gibbsCur.rho_r, arr,
														brr, rho_r_sd);
			logR = logR_rho_r() + mcstat::loglogitJacobian(gibbsCur.rho_r) -
				mcstat::loglogitJacobian(gibbsProposed.rho_r);
			accept = log(R::runif(0,1)) <= std::min(logR, 0.0);
			if(accept)
				rho_r_update();
			rho_r_accept += ((accept ? 1.0 : 0.0) - rho_r_accept) / (double) (it + 1);
			
			// save parameter
			rho_r_samples.at(it) = gibbsCur.rho_r;
			
			
			//
			// RW sigmasq_r
			//
			
			step = 'Z';
			
			// MH sample and update
			gibbsProposed.sigmasq_r = mcstat::logProposal(gibbsCur.sigmasq_r,
														  sigmasq_r_sd);
			logR = logR_sigmasq_r() + mcstat::loglogJacobian(gibbsCur.sigmasq_r) -
				mcstat::loglogJacobian(gibbsProposed.sigmasq_r) +
				mcstat::logdinvgamma_unscaled(gibbsProposed.sigmasq_r, ar, br) -
				mcstat::logdinvgamma_unscaled(gibbsCur.sigmasq_r, ar, br);
			accept = log(R::runif(0,1)) <= std::min(logR, 0.0);
			if(accept)
				sigmasq_r_update();
			sigmasq_r_accept += ((accept ? 1.0 : 0.0) - sigmasq_r_accept) / (double) (it + 1);
			
			// save parameter
			sigmasq_r_samples.at(it) = gibbsCur.sigmasq_r;
			
			
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
				logexp = gibbsCur.resid.t() * mcstat::dgemkmm(gibbsCur.C ,
							gibbsCur.SigmaInv, gibbsCur.resid);
				
				ll_samples.at(it) = ( - nt * gibbsCur.logdet_Sigma +
									    ns * gibbsCur.logdet_C -
									    logexp.at(0) ) / 2.0;
			}
			
			
			
			//
			// adjust tuning
			//
			
			adaptSize = C / sqrt( (double) (it + 1));
			rho_y_sd *= exp( adaptSize * (rho_y_accept - RWrate) );
			rho_r_sd *= exp( adaptSize * (rho_r_accept - RWrate) );
			eps_sd *= exp( adaptSize * (eps_accept - RWrate) );
			sigmasq_r_sd *= exp( adaptSize * (sigmasq_r_accept - RWrate) );
			
			
			
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
				_["sigmasq_r"] = gibbsCur.sigmasq_r,
				_["rho_y"] = gibbsCur.rho_y,
				_["rho_r"] = gibbsCur.rho_r,
				_["sigmasq_eps"] = gibbsCur.sigmasq_eps,
				_["R"] = gibbsCur.R,
				_["CInv"] = gibbsCur.CInv,
				_["Sigma"] = gibbsCur.Sigma
			);
			
			List proposed = List::create(
				_["beta"] = gibbsProposed.beta,
				_["sigmasq_y"] = gibbsProposed.sigmasq_y,
				_["sigmasq_r"] = gibbsProposed.sigmasq_r,
				_["rho_y"] = gibbsProposed.rho_y,
				_["rho_r"] = gibbsProposed.rho_r,
				_["sigmasq_eps"] = gibbsProposed.sigmasq_eps,
				_["R"] = gibbsProposed.R,
				_["CInv"] = gibbsProposed.CInv,
				_["Sigma"] = gibbsProposed.Sigma
			);
			
			List samples = List::create(
				_["beta"] = beta_samples,
				_["sigmasq_y"] = sigmasq_y_samples,
				_["sigmasq_r"] = sigmasq_r_samples,
				_["rho_y"] = rho_y_samples,
				_["rho_r"] = rho_r_samples,
				_["sigmasq_eps"] = sigmasq_eps_samples,
				_["ll"] = ll_samples
			);
			
			List tuning = List::create(
				_["rho_y_sd"] = rho_y_sd,
				_["rho_r_sd"] = rho_r_sd,
				_["eps_sd"] = eps_sd,
				_["sigmasq_r_sd"] = sigmasq_r_sd
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
		Rcout << "  sigmasq_r: " << round(sigmasq_r_accept * 100.0) << "%" << endl;
		Rcout << "  sigmasq_eps: " << round(eps_accept * 100.0) << "%" << endl;
		Rcout << "  rho_y: " << round(rho_y_accept * 100.0) << "%" << endl;
		Rcout << "  rho_r: " << round(rho_r_accept * 100.0) << "%" << endl;
	}
	
	void sampleAlphas(int burn) {
		/* Assuming that the sampler output variables have data, this will use
		   those data to composition sample alphas.
		 */
		
		
		// initialize constants
		mat ZZT = Z * Z.t();
		int nsr = ns * r;
		
		// initialize structures
		mat Sigma = mat(ns, ns, fill::zeros);
		mat R = mat(r, r, fill::zeros);
		mat RZZInv = mat(r, r, fill::zeros);
		vec muAlpha = vec(nsr, fill::zeros);
		vec alpha = vec(nsr, fill::zeros);
		
		// initialize output
		int maxIt = beta_samples.n_rows;
		int nsamples = maxIt - burn;
		alpha_samples = mat(nsamples, nsr, fill::zeros);
		
		// initialize tracking variables
		int checkpointIt = (int) maxIt * 0.1;
		std::clock_t start, lap, tmp_clock;
		start = std::clock();
		lap = std::clock();
		double duration, total, pctComplete, remaining;
		
		
		// generate composition samples
		for(int it=burn; it<maxIt; it++) {
			
			//
			// covariance structures
			//
			
			maternCov( Sigma, Dy, sigmasq_y_samples.at(it), rho_y_samples.at(it),
					   nu_y, sigmasq_y_samples.at(it) * sigmasq_eps_samples.at(it) );
			
			maternCov( R, Dz, sigmasq_r_samples.at(it), rho_r_samples.at(it),
					   nu_r, 0.0 );
			//poweredExpCov( R, Dz, sigmasq_r_samples.at(it), rho_r_samples.at(it),
			//		  nu_r, 0.0 );
			
			RZZInv = inv_sympd( inv_sympd(R) + ZZT );
			
			//
			// mean
			//
			
			muAlpha.zeros();
			
			for(int i=0; i<nt; i++) {
				int rStart = i * ns;
				int rEnd = rStart + ns_minusOne;
				muAlpha += kron(
					Y.rows(rStart, rEnd) - X.rows(rStart, rEnd) *
													 beta_samples.row(it).t(),
					RZZInv * Z.col(i)
				);
			}
			
			
			alpha = muAlpha + vectorise( chol(RZZInv, "lower") *
										 randn<mat>(r,ns) *
										 chol(Sigma, "upper") );
			
			alpha_samples.row(it-burn) = alpha.t();
			
			
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
	}
	
	
};


//
// R exports
//

RcppExport SEXP _stfit( SEXP p, SEXP r, SEXP ns, SEXP nt, SEXP X, SEXP Z,
					    SEXP Y, SEXP Dy, SEXP Dz, SEXP nu_y, SEXP nu_r,
					    SEXP ay, SEXP by, SEXP ar, SEXP br, SEXP ary, SEXP bry,
					    SEXP arr, SEXP brr, SEXP aeps, SEXP beps, SEXP Lambda,
					    SEXP rho_y_sd, SEXP rho_r_sd, SEXP eps_sd,
					    SEXP sigmasq_r_sd, SEXP maxIt, SEXP returnll,
					    SEXP errDump, SEXP C, SEXP RWrate ) {

	using namespace Rcpp;
	
	// instantiate sampler
	STModel stmod = STModel(as<int>(p), as<int>(r), as<int>(ns), as<int>(nt));
	
	
	//
	// configure sampler
	//
	
	stmod.setData(as<mat>(X), as<mat>(Z), as<vec>(Y), as<mat>(Dy), as<mat>(Dz));
	
	stmod.setPriors(as<double>(nu_y), as<double>(nu_r), as<double>(ay),
					as<double>(by), as<double>(ar), as<double>(br), as<double>(ary),
					as<double>(bry), as<double>(arr), as<double>(brr),
					as<double>(aeps), as<double>(beps), as<mat>(Lambda));
	
	stmod.tuneSampler(as<double>(rho_y_sd), as<double>(rho_r_sd),
					  as<double>(eps_sd), as<double>(sigmasq_r_sd));
	
	
	// run sampler
	stmod.fit(as<int>(maxIt), as<bool>(returnll), as<Function>(errDump),
			  as<double>(C), as<double>(RWrate));
	
	return List::create(
						_["beta"] = stmod.beta_samples,
						_["sigmasq_y"] = stmod.sigmasq_y_samples,
						_["sigmasq_r"] = stmod.sigmasq_r_samples,
						_["rho_y"] = stmod.rho_y_samples,
						_["rho_r"] = stmod.rho_r_samples,
						_["sigmasq_eps"] = stmod.sigmasq_eps_samples,
						_["ll"] = stmod.ll_samples
						);

}


RcppExport SEXP _compositionAlpha( SEXP p, SEXP r, SEXP ns, SEXP nt,
								   SEXP X, SEXP Z, SEXP Y, SEXP Dy, SEXP Dz,
								   SEXP nu_y, SEXP nu_r, SEXP beta,
								   SEXP sigmasq_y, SEXP rho_y, SEXP rho_r,
								   SEXP sigmasq_r, SEXP sigmasq_eps, SEXP burn,
								   SEXP summaryOnly) {
	
	using namespace Rcpp;
	
	int r_ = as<int>(r);
	int ns_ = as<int>(ns);
	int burn_ = as<int>(burn);
	
	// instantiate sampler
	STModel stmod = STModel(as<int>(p), r_, ns_, as<int>(nt));
	
	//
	// configure sampler
	//
	
	stmod.setData(as<mat>(X), as<mat>(Z), as<vec>(Y), as<mat>(Dy), as<mat>(Dz));
	
	stmod.setPriors(as<double>(nu_y), as<double>(nu_r), 0.0,
					0.0, 0.0, 0.0, 0.0,
					0.0, 0.0, 0.0,
					0.0, 0.0, mat(1,1,fill::eye) );
	
	//
	// add posterior samples
	//
	
	stmod.beta_samples = as<mat>(beta);
	stmod.sigmasq_y_samples = as<vec>(sigmasq_y);
	stmod.rho_y_samples = as<vec>(rho_y);
	stmod.rho_r_samples = as<vec>(rho_r);
	stmod.sigmasq_r_samples = as<vec>(sigmasq_r);
	stmod.sigmasq_eps_samples = as<vec>(sigmasq_eps);
	
	
	// run composition sampler
	stmod.sampleAlphas(burn_);
	
	
	//
	// compute summary objects
	//
	
	// posterior mean and sd for each alpha
	mat est = mean(stmod.alpha_samples);
	mat sd = stddev(stmod.alpha_samples, 1);
	
	// posterior variance for the remote field of alphas for each local point
	cube ZVar = cube(r_, r_, ns_, fill::zeros);
	for(int i=0; i<ns_; i++) {
		int rStart = i*r_;
		int rEnd = rStart + r_ - 1;
		ZVar.slice(i) = cov( stmod.alpha_samples.cols(rStart, rEnd), 1 );
	}
	
	// posterior covariance between beta and alpha
	mat covBetaAlpha = cov( stmod.beta_samples.rows(burn_, stmod.beta_samples.n_rows-1),
						    stmod.alpha_samples, 1 );
	
	// mean of the betas used required for merging covariances
	mat betaMean = mean(stmod.beta_samples.rows(burn_, stmod.beta_samples.n_rows-1));

	// build and return results
	List res;
	if(as<bool>(summaryOnly)==true) {
		res = List::create(
			_["est"] = est,
			_["sd"] = sd,
			_["ZVar"] = ZVar,
			_["covBetaAlpha"] = covBetaAlpha,
			_["beta"] = betaMean,
			_["nsamples"] = stmod.alpha_samples.n_rows
		);
	} else {
		res = List::create(
			_["est"] = est,
			_["sd"] = sd,
			_["ZVar"] = ZVar,
			_["covBetaAlpha"] = covBetaAlpha,
			_["beta"] = betaMean,
			_["samples"] = stmod.alpha_samples,
			_["nsamples"] = stmod.alpha_samples.n_rows
		);
		
	}
	return wrap( res );
}