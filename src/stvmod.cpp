#include <RcppArmadillo.h>
#include "mcstat.h"
#include "covs.h"


using namespace Rcpp;
using namespace arma;

//
// spatial teleconnection model with spatially varying local coefficients
//

// NOTES:
//	1) This class assumes that R will pass in the local covariates in a dense
//		matrix.  The class function, STVModel.setData will convert this to a
//		sparse, stacked block diagonal matrix.

class STVModel {
	
	private:
	
	// "constants"
	double one;
	
	// dimensions
	int p, r, ns, nt, N, nsp, ns_minusOne, nsnu_t, nt0, N0;
	
	// data
	SpMat<double> X, Xnew;
	mat Z, Dy, Dz, Znew;
	vec Y;
	
	// priors
	double nu_y, nu_r, ay, by, ar, br, ary, bry, arr, brr, aeps, beps,
		   sigmasq_y_shape, nu_t;
	mat Psi;
	
	// tuning parameters
	double rho_y_sd, rho_r_sd, eps_sd, sigmasq_r_sd, rho_y_accept, rho_r_accept,
		   eps_accept, sigmasq_r_accept;
	
	// nested class holds model parameters and dependent structures
	class Params {
		
		public:
		
		// parameters
		vec beta;
		double sigmasq_y, rho_y, rho_r, sigmasq_r, sigmasq_eps;
		mat T;
		
		// derived terms
		vec resid;
		mat R, C, CInv, Sigma, SigmaInv, TInv;
		double logdet_Sigma, logdet_C;
		
		// constructor zero-initializes parameters
		Params(int p=1, int r=1, int ns=1, int nt=1, int N=1, int nsp=1) {
			beta = vec(nsp, fill::zeros);
			resid = vec(N, fill::zeros);
			
			R = mat(r, r, fill::zeros);
			CInv = mat(r, r, fill::zeros);
			C = mat(r, r, fill::zeros);
			
			Sigma = mat(ns, ns, fill::zeros);
			SigmaInv = mat(ns, ns, fill::zeros);
			
			T = mat(p, p, fill::zeros);
			TInv = mat(p, p, fill::zeros);
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
	cube T_samples;
	vec sigmasq_y_samples, rho_y_samples, rho_r_samples, sigmasq_r_samples,
		sigmasq_eps_samples, ll_samples;
	
	// composition sample output
	mat alpha_samples;
	cube fcst_samples, local_samples, remote_samples;
	
	
	
	//
	// class functions
	//
	
	STVModel(int _p, int _r, int _ns, int _nt) {
		p = _p;		// number of mean parameters
		r = _r;		// number of remote locations
		ns = _ns;	// number of local locations
		nt = _nt;	// number of time points
		
		// compute additional dimension-based constants
		N = ns * nt;
		nsp = ns * p;
		ns_minusOne = ns - 1;
		
		// initialize state
		gibbsCur = Params(p, r, ns, nt, N, nsp);
		gibbsProposed = Params(p, r, ns, nt, N, nsp);
		
		// initialize "constants"
		one = 1.0;
	}
	
	
	void setData(const mat & _X, const mat & _Z, const vec & _Y,
				 const mat & _Dy, const mat & _Dz) {
		Z = _Z;		// spatio-temporal remote covariate matrix (r x nt)
		Y = _Y;		// spatio-temporal response (ns*nt x 1)
		Dy = _Dy;	// local distance matrix
		Dz = _Dz;	// remote distance matrix
		
		// spatially varying coefficient design matrix (ns*nt x ns*p)
		X = sp_mat(repmat(linspace<uvec>(0, nsp-1, nsp), nt, 1),
				   linspace<uvec>(0, N*p, N+1),
				   vectorise(_X, 1),
				   nsp, N).t();
	}
	
	void setCompositionData(const mat & _Xnew, const mat & _Znew) {
		
		//
		// set data needed for computing forecasts
		//	 (recall that there may be a different number of prediction timepoints)
		//
		
		// extract number of forecast timepoints, and set dimension constants
		nt0 = _Znew.n_cols;
		N0 = ns * nt0;
		
		Znew = _Znew;	// spatio-temporal remote covariate matrix (r x nt0)
		
		// spatially varying coefficient design matrix (ns*nt0 x ns*p)
		Xnew = sp_mat(repmat(linspace<uvec>(0, nsp-1, nsp), nt0, 1),
				   linspace<uvec>(0, N0*p, N0+1),
				   vectorise(_Xnew, 1),
				   nsp, N0).t();
	}
	
	
	void setPriors(double _nu_y, double _nu_r, double _ay, double _by,
				   double _ar, double _br, double _ary, double _bry,
				   double _arr, double _brr, double _aeps, double _beps,
				   const mat & _Psi, double _nu_t) {
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
		Psi = _Psi;			// prior parameters for mean process
		nu_t = _nu_t;
		
		
		//
		// compute constants associated with priors
		//
		
		// conjugate sampling sigmasq_y
		sigmasq_y_shape = ((double) N / 2.0) + ay;
		
		// conjugate sampling T
		nsnu_t = ns + nu_t;
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
		beta_samples = mat(maxIt, nsp, fill::zeros);
		T_samples = cube(p, p, maxIt, fill::zeros);
		sigmasq_y_samples = vec(maxIt, fill::zeros);
		rho_y_samples = vec(maxIt, fill::zeros);
		rho_r_samples = vec(maxIt, fill::zeros);
		sigmasq_r_samples = vec(maxIt, fill::zeros);
		sigmasq_eps_samples = vec(maxIt, fill::zeros);
		if(returnll)
			ll_samples = vec(maxIt, fill::zeros);
		
		// initialize parameters
		gibbsCur.T = mcstat::rinvwishart(Psi, nu_t);
		gibbsCur.TInv = inv_sympd(gibbsCur.T);
		gibbsCur.beta = repmat(mcstat::mvrnorm(gibbsCur.T), ns, 1);
		gibbsCur.sigmasq_y = 1.0 / R::rgamma(ay, 1.0 / by);
		gibbsCur.rho_y = R::runif(ary, bry);
		gibbsCur.rho_r = R::runif(arr, brr);
		gibbsCur.sigmasq_r = 1.0 / R::rgamma(ar, 1.0 / br);
		gibbsCur.sigmasq_eps = 1.0 / R::rgamma(ay, 1.0 / by);
		
		// initialize temp
		double logR = 0.0;
		double adaptSize;
		bool accept;
		mat CkSigmaInvX, postBetaCov, conjT;
		vec logexp, postBetaMean, bs;
		conjT = mat(p, p, fill::zeros);
		
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
			CkSigmaInvX = mcstat::dsemkmm(gibbsCur.C, gibbsCur.SigmaInv, X); //update for sparse X?
			postBetaCov = inv_sympd( kron(gibbsCur.SigmaInv, gibbsCur.TInv) + X.t() * CkSigmaInvX );
			postBetaMean = postBetaCov * CkSigmaInvX.t() * Y;
			
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
			// conjugate T
			//
			
			step='T';
			
			// update posterior
			
			// reset matrix
			conjT.zeros();
			// compute cross terms
			for(int i=0; i<ns; i++) {
				bs = gibbsCur.beta.rows(i*p, (i+1)*p-1);
				for(int j=0; j<i; j++) {
					conjT = conjT + gibbsCur.SigmaInv.at(i,j) *
						bs * gibbsCur.beta.rows(j*p, (j+1)*p-1).t();
				}
			}
			// add symmetry
			conjT = conjT + conjT.t();
			// add diagonal entries
			for(int i=0; i<ns; i++) {
				bs = gibbsCur.beta.rows(i*p, (i+1)*p-1);
				conjT = conjT + gibbsCur.SigmaInv.at(i,i) * bs * bs.t();
			}
			
			// sample and save
			gibbsCur.T = mcstat::rinvwishart(conjT + Psi, nsnu_t);
			gibbsCur.TInv = inv_sympd(gibbsCur.T);
			T_samples.slice(it) = gibbsCur.T;
			
			
			
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
			// adjust tuning (following a method presented in Andrieu and Thoms (2008))
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
				_["Sigma"] = gibbsCur.Sigma,
				_["T"] = gibbsCur.T
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
				_["Sigma"] = gibbsProposed.Sigma,
				_["T"] = gibbsProposed.T
			);
			
			List samples = List::create(
				_["beta"] = beta_samples,
				_["sigmasq_y"] = sigmasq_y_samples,
				_["sigmasq_r"] = sigmasq_r_samples,
				_["rho_y"] = rho_y_samples,
				_["rho_r"] = rho_r_samples,
				_["sigmasq_eps"] = sigmasq_eps_samples,
				_["T"] = T_samples,
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
	
	
	void compositionSample(int burn, bool alphas, bool forecast) {
		/* Assuming that the sampler output variables have data, this will use
		 those data to composition sample alphas and then use these to forecast
		 the response at new locations.
		 */
		
		// initialize constants
		mat ZZT = Z * Z.t();
		int nsr = ns * r;
		
		// initialize structures
		mat Sigma = mat(ns, ns, fill::zeros);
		mat Sigma_cholU = mat(ns, ns, fill::zeros);
		mat R = mat(r, r, fill::zeros);
		mat RZZInv = mat(r, r, fill::zeros);
		vec muAlpha = vec(nsr, fill::zeros);
		vec alpha;
		mat local, remote, noise;
		if(alphas) {
			alpha = vec(nsr, fill::zeros);
		}
		if(forecast) {
			local = mat(ns, nt0, fill::zeros);
			remote = mat(ns, nt0, fill::zeros);
			noise = mat(ns, nt0, fill::zeros);
		}
		
		// TODO: allow forecasts to be made at new or subsampled locations
		
		// initialize output
		int maxIt = beta_samples.n_rows;
		int nsamples = maxIt - burn;
		if(alphas) {
			alpha_samples = mat(nsamples, nsr, fill::zeros);
		}
		if(forecast) {
			fcst_samples = cube(ns, nt0, nsamples, fill::zeros);
			local_samples = cube(ns, nt0, nsamples, fill::zeros);
			remote_samples = cube(ns, nt0, nsamples, fill::zeros);
		}
		
		// initialize tracking variables
		int checkpointIt = (int) maxIt * 0.1;
		std::clock_t start, lap, tmp_clock;
		start = std::clock();
		lap = std::clock();
		double duration, total, pctComplete, remaining;
		
		
		// generate composition samples
		for(int it=burn; it<maxIt; it++) {
			
			checkUserInterrupt();
			
			//
			// covariance structures
			//
			
			maternCov( Sigma, Dy, sigmasq_y_samples.at(it), rho_y_samples.at(it),
					  nu_y, sigmasq_y_samples.at(it) * sigmasq_eps_samples.at(it) );
			
			Sigma_cholU = chol(Sigma, "upper");
			
			maternCov( R, Dz, sigmasq_r_samples.at(it), rho_r_samples.at(it),
					  nu_r, 0.0 );
			
			RZZInv = inv_sympd( inv_sympd(R) + ZZT );
			
			//
			// posterior mean for teleconnection effects
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
			
			
			// sample teleconnection effects
			alpha = muAlpha + vectorise( chol(RZZInv, "lower") *
										 randn<mat>(r,ns) *
										 Sigma_cholU );
			
			// save teleconnection effects
			if(alphas)
				alpha_samples.row(it-burn) = alpha.t();
			
			
			// sample and save posterior predictive samples
			if(forecast) {
				
				// compute remote effects
				remote = reshape(alpha, r, ns).t() * Znew;
				
				// compute local effects
				local = reshape(Xnew * beta_samples.row(it).t(), ns, nt0);
				
				// sample spatially correlated noise, independent across time
				noise = Sigma_cholU.t() * randn<mat>(ns, nt0);
				
				// save forecast objects
				local_samples.slice(it-burn) = local;
				remote_samples.slice(it-burn) = remote;
				fcst_samples.slice(it-burn) = local + remote + noise;
			}
			
			
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

RcppExport SEXP _stvfit( SEXP p, SEXP r, SEXP ns, SEXP nt, SEXP X, SEXP Z,
					     SEXP Y, SEXP Dy, SEXP Dz, SEXP nu_y, SEXP nu_r,
					     SEXP ay, SEXP by, SEXP ar, SEXP br, SEXP ary, SEXP bry,
					     SEXP arr, SEXP brr, SEXP aeps, SEXP beps, SEXP Psi,
					     SEXP rho_y_sd, SEXP rho_r_sd, SEXP eps_sd,
					     SEXP sigmasq_r_sd, SEXP maxIt, SEXP returnll,
					     SEXP errDump, SEXP C, SEXP RWrate, SEXP nu_t ) {

	using namespace Rcpp;
	
	// instantiate sampler
	STVModel stvmod = STVModel(as<int>(p), as<int>(r), as<int>(ns), as<int>(nt));
	
	
	//
	// configure sampler
	//
	
	stvmod.setData(as<mat>(X), as<mat>(Z), as<vec>(Y), as<mat>(Dy), as<mat>(Dz));
	
	stvmod.setPriors(as<double>(nu_y), as<double>(nu_r), as<double>(ay),
					as<double>(by), as<double>(ar), as<double>(br), as<double>(ary),
					as<double>(bry), as<double>(arr), as<double>(brr),
					as<double>(aeps), as<double>(beps), as<mat>(Psi),
					 as<double>(nu_t));
	
	stvmod.tuneSampler(as<double>(rho_y_sd), as<double>(rho_r_sd),
					  as<double>(eps_sd), as<double>(sigmasq_r_sd));
	
	
	// run sampler
	stvmod.fit(as<int>(maxIt), as<bool>(returnll), as<Function>(errDump),
			  as<double>(C), as<double>(RWrate));
	
	return List::create(
						_["beta"] = stvmod.beta_samples,
						_["sigmasq_y"] = stvmod.sigmasq_y_samples,
						_["sigmasq_r"] = stvmod.sigmasq_r_samples,
						_["rho_y"] = stvmod.rho_y_samples,
						_["rho_r"] = stvmod.rho_r_samples,
						_["sigmasq_eps"] = stvmod.sigmasq_eps_samples,
						_["T"] = stvmod.T_samples,
						_["ll"] = stvmod.ll_samples
						);
	 
	
	return List::create();

}


RcppExport SEXP _stvcomposition( SEXP p, SEXP r, SEXP ns, SEXP nt,
								 SEXP X, SEXP Z, SEXP Y, SEXP Dy, SEXP Dz,
								 SEXP nu_y, SEXP nu_r, SEXP beta,
								 SEXP sigmasq_y, SEXP rho_y, SEXP rho_r,
								 SEXP sigmasq_r, SEXP sigmasq_eps,
								 SEXP burn, SEXP summaryOnly, SEXP Xnew,
								 SEXP Znew, SEXP alphas, SEXP forecast ) {
	
	using namespace Rcpp;
	
	int r_ = as<int>(r);
	int ns_ = as<int>(ns);
	int burn_ = as<int>(burn);
	bool alphas_ = as<bool>(alphas);
	bool forecast_ = as<bool>(forecast);
	
	// instantiate sampler
	STVModel stvmod = STVModel(as<int>(p), r_, ns_, as<int>(nt));
	
	//
	// configure sampler
	//
	
	stvmod.setData(as<mat>(X), as<mat>(Z), as<vec>(Y), as<mat>(Dy), as<mat>(Dz));
	
	stvmod.setPriors(as<double>(nu_y), as<double>(nu_r), 0.0,
					0.0, 0.0, 0.0, 0.0,
					0.0, 0.0, 0.0,
					0.0, 0.0, mat(1,1,fill::eye), 0.0 );
	
	if(forecast_)
		stvmod.setCompositionData(as<mat>(Xnew), as<mat>(Znew));
	
	//
	// add posterior samples
	//
	
	stvmod.beta_samples = as<mat>(beta);
	stvmod.sigmasq_y_samples = as<vec>(sigmasq_y);
	stvmod.rho_y_samples = as<vec>(rho_y);
	stvmod.rho_r_samples = as<vec>(rho_r);
	stvmod.sigmasq_r_samples = as<vec>(sigmasq_r);
	stvmod.sigmasq_eps_samples = as<vec>(sigmasq_eps);
	
	
	// run composition sampler
	stvmod.compositionSample(burn_, alphas_, forecast_);
	
	// post-process teleconnection effects
	mat est, sd;
	if(alphas_) {
		
		//
		// compute summary objects
		//
		
		// posterior mean and sd for each alpha
		est = mean(stvmod.alpha_samples);
		sd = stddev(stvmod.alpha_samples, 1);
		
	}
	

	//
	// build and return results
	//
	
	// compile teleconnection results
	List alpha_res;
	if(alphas_) {
		if(as<bool>(summaryOnly)==true) {
			alpha_res = List::create(
							_["est"] = est,
							_["sd"] = sd,
							_["nsamples"] = stvmod.alpha_samples.n_rows
						);
		} else {
			alpha_res = List::create(
							_["est"] = est,
							_["sd"] = sd,
							_["samples"] = stvmod.alpha_samples,
							_["nsamples"] = stvmod.alpha_samples.n_rows
						);
		}
	}
	
	// compile forecast results
	List forecast_res;
	if(forecast_) {
		forecast_res = List::create(
						_["forecasts"] = stvmod.fcst_samples,
						_["local"] = stvmod.local_samples,
						_["remote"] = stvmod.remote_samples
		);
	}
	
	// compile appropriate results
	List res;
	if(alphas_ & forecast_) {
		res = List::create(
				_["forecast"] = forecast_res,
				_["alpha"] = alpha_res
		);
	} else if(alphas_) {
		res = List::create(
				_["alpha"] = alpha_res
		);
	} else {
		res = List::create(
				_["forecast"] = forecast_res
		);
	}
	
	return wrap( res );
}