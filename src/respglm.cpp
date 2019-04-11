/*
	Implementation of Gibbs sampler for RESP-GLM parameters
*/

#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include "covs.h"
#include "BlockRWSampler.h"
#include "glm_gmrf.h"
#include "distributions.h"

using namespace Rcpp;
using namespace arma;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::SparseMatrix;
using Eigen::SimplicialLLT;

// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<double> EigenSpMat;

struct Data {
	VectorXd *Y;
	MatrixXd *X;
};

struct Params {
	double nu, // Matern smoothness for ocean
				 rho, // Matern range for ocean
				 sigmasq, // Matern sill for ocean
				 kappa; // correlation parameter for GMRF structure
	VectorXd beta, // regression coefficients
					 eta0; // loaded teleconnection effects
};

struct Consts {
	mat *dknots,  // matrix of distances between ocean knots (k x k)
			*dzknots; // matrix of distances from ocean locations to knots (nz x k)
	MatrixXd *W, *A;  // EOF patterns and scores (nz x neofs), (neofs x t)
	int nz, nEOFs, nknots, // number of ocean locations, EOF patterns, knots
			n, nt, p; // number of areal units, timepoints, covariates
	SparseMatrix<double> *Q; // GMRF structure for areal units
};

struct Priors {
	double rho_L, rho_U, kappa_L, kappa_U; // bounds for spatial cov. priors
};

struct Scratch {
	mat R, c; // spatial structures for ocean covariance matrices
	MatrixXd psi; // spatial structure for ocean covariance matrices
	EigenSpMat SigmaInv; // precision matrix for areal units
	double lly, lleta0, llbeta, lltheta; // log-likelihood components
	VectorXd eta0mu; // mean for gaussian approx. to full conditional posterior
};

struct Config {
	Data data; Params params; Consts consts; Priors priors; Scratch scratch;
};


// Deprecated
// EigenSpMat buildEta0Prec(const EigenSpMat& SigmaInv, const MatrixXd& A,
// 	const MatrixXd& psi) {
// 	/* compute the prior precision matrix for eta0
//
// 		Parameters:
// 			prec - (output) prior precision matrix
// 			SigmaInv - precision matrix for areal units (n x n)
// 			A - EOF scores, each column should contain all scores for one timepoint
// 				(k x t)
// 			psi - prior covariance for transformed teleconnection effects
// 	*/
//
// 	int t = A.cols();
// 	return kroneckerProduct(
// 		(A.transpose() * psi * A + MatrixXd::Identity(t,t)).inverse(),
// 		SigmaInv
// 	);
// }

EigenSpMat buildEta0Prec(const EigenSpMat& SigmaInv, LLT<MatrixXd>& localLLT) {
	/* compute the prior precision matrix for eta0

		Parameters:
			prec - (output) prior precision matrix
			SigmaInv - precision matrix for areal units (n x n)
			A - EOF scores, each column should contain all scores for one timepoint
				(k x t)
			psi - prior covariance for transformed teleconnection effects
	*/

	int t = localLLT.matrixL().rows();
	MatrixXd localPrec = localLLT.solve(MatrixXd::Identity(t,t));

	return kroneckerProduct(localPrec, SigmaInv);
}


LLT<MatrixXd> decomposeLocalCov(const MatrixXd& A, mat& R, mat& c,
	const mat& d, const mat& dz, const MatrixXd& W, double scale, double range,
	double smoothness) {
	/* Compute LLT decomposition of A^t Psi A + I

		 Parameters:
		 	 A - EOF scores, each column should contain all scores for one timepoint
				   (k x t)
			 R - pre-allocated storage for prior covariance matrix for teleconnection
 					 effects at knot locations (k x k)
	 		 c - pre-allocated storage for prior cross-covariances (nz x k)
	 		 d - matrix of distances between knots (k x k)
	 		 dz - matrix of distances from ocean locations to knots (nz x k)
	 		 W - matrix of EOF patterns (nz x neofs)
	 		 scale - Matern sill parameter
	 		 range - Matern range parameter
	 		 smoothness - Matern smoothness parameter
	*/

	// build and decompose prior covariance for teleconnection effects at knots
	maternCov(R, d, scale, range, smoothness, 0);
	Map<MatrixXd> rmat(R.memptr(), R.n_rows, R.n_cols);
	LLT<MatrixXd> lltR(rmat);

	// build and compute mapping onto EOF patterns
	maternCov(c, dz, scale, range, smoothness, 0);

	// compute y, a decomposition of psi
	Map<MatrixXd> cmat(c.memptr(), c.n_rows, c.n_cols);
	MatrixXd y = lltR.matrixL().solve(cmat.transpose() * W);

	int t = A.cols();
	LLT<MatrixXd> localLLT(
		A.transpose() * y.transpose() * y * A + MatrixXd::Identity(t,t)
	);

	return localLLT;
}

// Deprecated
MatrixXd buildPsi(mat& R, mat& c, const mat& d, const mat& dz,
	const MatrixXd& W, double scale, double range, double smoothness) {
	/* compute entries in prior covariance matrix psi for transformed
		 teleconnection effects. also updates R and c as side effects.

		 Parameters:
		   psi - (output) prior covariance matrix
			 R - pre-allocated storage for prior covariance matrix for teleconnection
			 			effects at knot locations (k x k)
			 c - pre-allocated storage for prior cross-covariances (nz x k)
			 d - matrix of distances between knots (k x k)
			 dz - matrix of distances from ocean locations to knots (nz x k)
			 W - matrix of EOF patterns (nz x neofs)
			 scale - Matern sill parameter
			 range - Matern range parameter
			 smoothness - Matern smoothness parameter
	*/

	// build and decompose prior covariance for teleconnection effects at knots
	maternCov(R, d, scale, range, smoothness, 0);
	Map<MatrixXd> rmat(R.memptr(), R.n_rows, R.n_cols);
	LLT<MatrixXd> lltR(rmat);

	// build and compute mapping onto EOF patterns
	maternCov(c, dz, scale, range, smoothness, 0);

	Map<MatrixXd> cmat(c.memptr(), c.n_rows, c.n_cols);
	MatrixXd y = lltR.matrixL().solve(cmat.transpose() * W);

	return y.transpose() * y;
}


// TODO: write function that computes likelihoods

class TeleSampler : public mcstat2::BlockRWSampler {

	public:

		/* One-block sampler to update covariance parameters for teleconnection
			 effects.  Will also update the teleconnection effects to help facilitate
			 larger step sizes.  Teleconnection effects can be saved via additional
			 "extraction" samplers. */

		TeleSampler(Config& t_cfg, std::vector<double> sds,
			std::vector<double> init, std::vector<double> C, double alpha,
			std::vector<double> priors_L, std::vector<double> priors_U) :
		BlockRWSampler({"kappa", "sigmasq", "rho"}, {LOGIT, LOG, LOGIT},
			sds, init, C, alpha, priors_L, priors_U) {
				cfg = &t_cfg;
				prop = Config(t_cfg);
				R = mat(cfg->consts.nknots, cfg->consts.nknots, fill::zeros);
				c = mat(cfg->consts.nz, cfg->consts.nknots, fill::zeros);
				eta0mu = VectorXd(cfg->consts.n * cfg->consts.nt);
		}

	private:

		Config *cfg, prop;
		mat R, c;
		VectorXd eta0mu;

		double logR_posterior(const std::vector<double>& x,
													const std::vector<double>& x0) {

//prop.scratch.lltheta


			// pre-compute some of the proposal ll components
			// prop.scratch.lleta0


			// return likelihood ratio!
			return 0;
		}

		void preAcceptProb(const std::vector<double>& x) {

			// compute covariance matrices

			SparseMatrix<double> SigmaInv = *prop.consts.Q * x[0];

			// MatrixXd psi = buildPsi(R, c, *(prop.consts.dknots),
			// 	*(prop.consts.dzknots), *(prop.consts.W), x[1], x[2], cfg->params.nu);

			LLT<MatrixXd> localLLT = decomposeLocalCov(*(prop.consts.A), R, c,
				*(prop.consts.dknots), *(prop.consts.dzknots), *(prop.consts.W),
				x[1], x[2], cfg->params.nu);

			//EigenSpMat Qprior = buildEta0Prec(SigmaInv, *(prop.consts.A), psi);
			EigenSpMat Qprior = buildEta0Prec(SigmaInv, localLLT);

			// gaussian approximation to full conditional for teleconnection effects
			SimplicialLLT<EigenSpMat> eta0CovL;
			gaussian_approx_eta0(cfg->params.eta0.data(), Qprior, 1, eta0mu.data(),
				eta0CovL, cfg->params.beta.data(), prop.data.Y->data(),
				prop.data.X->data(), prop.consts.n, prop.consts.nt, prop.consts.p,
				mcstat2::glm::glmfamily::poisson);

			// compute transition probabilities

			// propose new teleconnection effects
			prop.params.eta0 = eta0mu + mcstat2::mvrnorm_spchol(eta0CovL);

			// compute likelihood
			prop.scratch.lly = mcstat2::glm::ll(prop.data.Y->data(),
				prop.params.eta0.data(), cfg->params.beta.data(), prop.data.X->data(),
				prop.consts.n, prop.consts.nt, prop.consts.p,
				mcstat2::glm::glmfamily::poisson);

			// TODO: figure out how to deal with the adaptation likelihoods. they
			// might just be approximations since we won't consider updating the eta0
			// each time as well since this involves some additional randomness

		}

		void update() {
			// save updated model structures

			// compute updates to the additional model structures
		}
};


//
// Rcpp exports
//

// TODO: include sparse matrix Q


// [[Rcpp::export]]
List respglm_fit(arma::mat& dknots, arma::mat& dzknots, Eigen::MatrixXd& W,
	int nSamples, List priors, std::vector<double>& inits,
	std::vector<double>& sds, std::vector<double>& C,
	Eigen::SparseMatrix<double>& Q, Eigen::VectorXd inits_eta0,
	Eigen::VectorXd inits_beta, Eigen::VectorXd& Y, Eigen::MatrixXd& X,
  Eigen::MatrixXd& A) {
		/* Sample covariance parameters and latent effects for RESP GLM model.

			 Parameters:
				 dknots - matrix of distances between ocean knots
				 dzknots - matrix of distances from ocean knots to ocean locations
				 W - matrix of EOF patterns
				 nSamples - number of posterior samples to draw
				 priors - Rcpp::List containing parameterizations of priors; see
				 	R documentation to get details of the list's structure
				inits - vector containing initial values for  (kappa, sigmasq, rho)
				sds - initial proposal sd's for random walk samplers
				C - scaling factors for adapting RW proposal sd's
				Q - GMRF structure for areal units
				inits_eta0 - vector of initial values of loaded teleconnection effects
				inits_beta - vector of initial values for regression coefficients
				Y - vector of all observations
				X - matrix of covariates
				A - matrix of eof scores, each col. should have all scores for one
					timepoint (neofs x t)
		*/

	/*    - Fixed effect samplers
						* Metropolis step where the proposal distribution is the Gaussian
							approximation proposal
				- Extraction samplers
					  * BlockSampler objects that extract the teleconnection effects and
							likelihood from the model; this is a deterministic step
	*/

	// extract data

	Config cfg = Config();

	cfg.data.Y = &Y;
	cfg.data.X = &X;

	cfg.consts.dknots = &dknots;
	cfg.consts.dzknots = &dzknots;
	cfg.consts.W = &W;
	cfg.consts.nz = W.rows();
	cfg.consts.nEOFs = W.cols();
	cfg.consts.nknots = dknots.n_rows;
	cfg.consts.Q = &Q;
	cfg.consts.n = Q.rows();
	cfg.consts.nt = Y.size() / cfg.consts.n;
	cfg.consts.p = X.cols();
	cfg.consts.A = &A;

	cfg.priors.rho_L = as<double>(priors["rho_L"]);
	cfg.priors.rho_U = as<double>(priors["rho_U"]);
	cfg.priors.kappa_L = as<double>(priors["kappa_L"]);
	cfg.priors.kappa_U = as<double>(priors["kappa_U"]);

	cfg.params.nu = as<double>(priors["nu"]);
	cfg.params.kappa = inits[0];
	cfg.params.rho = inits[1];
	cfg.params.sigmasq = inits[2];
	cfg.params.beta = inits_beta;
	cfg.params.eta0 = inits_eta0;

	// initialize scratch

	/*
	MatrixXd psi(cfg.consts.nEOFs, cfg.consts.nEOFs);
	mat R(cfg.consts.nknots, cfg.consts.nknots, fill::zeros);
	mat c(cfg.consts.nz, cfg.consts.nknots, fill::zeros);
	SparseMatrix<double> SigmaInv = *cfg.consts.Q * cfg.params.kappa;
	VectorXd eta0mu(cfg.consts.n * cfg.consts.nt);
	SimplicialLLT<EigenSpMat> eta0CovL;

	cfg.scratch.psi = psi;
	cfg.scratch.R = R;
	cfg.scratch.c = c;
	cfg.scratch.SigmaInv = SigmaInv;
	cfg.scratch.eta0mu = eta0mu;
	*/

	// buildPsi(psi, R, c, dknots, dzknots, W, cfg.params.sigmasq, cfg.params.rho,
	// 	cfg.params.nu);

	// initialize samplers

	TeleSampler ts = TeleSampler(cfg, sds, inits, C, .23,
		{cfg.priors.kappa_L, 0, cfg.priors.rho_L},
		{cfg.priors.kappa_U, 1, cfg.priors.rho_U});

	mcstat2::GibbsSampler sampler = mcstat2::GibbsSampler();
	sampler.addSampler(ts);
	sampler.run(nSamples);

	return sampler.getSamples();
}
