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
			n, nt, p, // number of areal units, timepoints, covariates
			df; // rank deficiency for GMRF structure Q
	SparseMatrix<double> *Q; // GMRF structure for areal units
};

struct Priors {
	double rho_L, rho_U, kappa_a, kappa_b; // bounds for spatial cov. priors
	double s_a, s_b; // prior parameters for sigmasq
};

struct Scratch {
	SparseMatrix<double> eta0Prec; // prior prec. matrix for eta0
	double lly, lleta0, llbeta, lltheta; // log-likelihood components
};

struct Config {
	Data data; Params params; Consts consts; Priors priors; Scratch scratch;
};


void buildEta0Prec(const EigenSpMat& SigmaInv,
	const LLT<MatrixXd>& localLLT, EigenSpMat& prec) {
	/* compute the prior precision matrix for eta0

		Parameters:
			SigmaInv - precision matrix for areal units (n x n)
			localLLT - LLT decomposition of A^t Psi A + I (eta0 covariance)
			prec - (output) prior precision matrix for eta0
	*/
	int t = localLLT.matrixL().rows();
	MatrixXd localPrec = localLLT.solve(MatrixXd::Identity(t,t));
	prec = kroneckerProduct(localPrec, SigmaInv);
}


void localCovLLT(const MatrixXd& A, mat& R, mat& c, const mat& d, const mat& dz,
	const MatrixXd& W, double scale, double range, double smoothness,
	LLT<MatrixXd>& llt) {
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

	// compute and decompose local covariance for eta0
	int t = A.cols();
	llt.compute(
		A.transpose() * y.transpose() * y * A + MatrixXd::Identity(t,t)
	);
}

class StateSampler : public mcstat2::BlockSampler {

	public:

		/* Extract likelihood and teleconnection effects, which are updated
			throughout sampling */

		StateSampler(Config& t_cfg) :
			BlockSampler({VECTOR, REAL}, {"eta0", "ll"}) {
				cfg = &t_cfg;
			}

		void drawSample() { };

		vec returnSamples(int i) {
			vec s;
			if(i==0) {
				s = vec(cfg->params.eta0.data(), cfg->params.eta0.size());
			} else if(i==1) {
				s = vec(1);
				s.at(0) = cfg->scratch.lly;
			}
			return s;
		}

		int getSize(int i) { return i==0 ? cfg->params.eta0.size() : 1; }

	private:

		Config *cfg;

};

class TeleSampler : public mcstat2::BlockRWSampler {

	public:

		/* One-block sampler to update covariance parameters for teleconnection
			 effects.  Will also update the teleconnection effects to help facilitate
			 larger step sizes.  Teleconnection effects can be saved via additional
			 "extraction" samplers. */

		TeleSampler(Config& t_cfg, std::vector<double> sds,
			std::vector<double> init, std::vector<double> C, double alpha,
			std::vector<double> priors_L, std::vector<double> priors_U, bool adapt) :
		BlockRWSampler({"kappa", "sigmasq", "rho"}, {LOG, LOG, LOGIT},
			sds, init, C, alpha, priors_L, priors_U, adapt) {

				// copy model configuration
				cfg = &t_cfg;
				prop = Config(t_cfg);

				// pre-allocate space for covariance matrix and approx. components
				R = mat(cfg->consts.nknots, cfg->consts.nknots, fill::zeros);
				c = mat(cfg->consts.nz, cfg->consts.nknots, fill::zeros);
				eta0mu = VectorXd(cfg->consts.n * cfg->consts.nt);

				// pre-compute eta0 precision
				SparseMatrix<double> SigmaInv = *(cfg->consts.Q) * cfg->params.kappa;
				localCovLLT(*(cfg->consts.A), R, c, *(cfg->consts.dknots),
					*(cfg->consts.dzknots), *(cfg->consts.W), cfg->params.sigmasq,
					cfg->params.rho, cfg->params.nu, localLLT);
				buildEta0Prec(SigmaInv, localLLT, cfg->scratch.eta0Prec);

				// pre-compute eta0 density
				cfg->scratch.lleta0 = mcstat2::ldigmrfKron(cfg->params.eta0,
					eta0mu, localLLT, cfg->params.kappa, cfg->consts.df,
					cfg->scratch.eta0Prec);
		}

	private:

		Config *cfg, prop;
		mat R, c;
		VectorXd eta0mu;
		LLT<MatrixXd> localLLT;
		double fwdll, backll;

		double logR_posterior(const std::vector<double>& x,
													const std::vector<double>& x0) {

			// TODO: for adaptation purposes, localLLT *should* be updated
			prop.scratch.lleta0 = mcstat2::ldigmrfKron(prop.params.eta0,
				eta0mu, localLLT, x[0], prop.consts.df, prop.scratch.eta0Prec);

			prop.scratch.lltheta =
				mcstat2::ldinvgamma(x[1], prop.priors.s_a, prop.priors.s_b) +
				mcstat2::ldinvgamma(x[0], prop.priors.kappa_a, prop.priors.kappa_b);

			return prop.scratch.lly + prop.scratch.lleta0 + prop.scratch.lltheta +
			  backll -
				(cfg->scratch.lly + cfg->scratch.lleta0 + cfg->scratch.lltheta + fwdll);
		}

		void preAcceptProb(const std::vector<double>& x) {

			// extract proposals

			prop.params.kappa = x[0];
			prop.params.sigmasq = x[1];
			prop.params.rho = x[2];

			// compute proposed covariance matrices for eta0

			SparseMatrix<double> SigmaInv = *prop.consts.Q * prop.params.kappa;

			localCovLLT(*(prop.consts.A), R, c, *(prop.consts.dknots),
				*(prop.consts.dzknots), *(prop.consts.W), prop.params.sigmasq,
				prop.params.rho, cfg->params.nu, localLLT);

			buildEta0Prec(SigmaInv, localLLT, prop.scratch.eta0Prec);

			// gaussian approximation to full conditional for teleconnection effects
			SimplicialLLT<EigenSpMat> eta0CovL;
			gaussian_approx_eta0(cfg->params.eta0.data(), prop.scratch.eta0Prec, 1,
				eta0mu.data(), eta0CovL, cfg->params.beta.data(), prop.data.Y->data(),
				prop.data.X->data(), prop.consts.n, prop.consts.nt, prop.consts.p,
				mcstat2::glm::glmfamily::poisson);

			// propose new teleconnection effects
			prop.params.eta0 = eta0mu + mcstat2::mvrnorm_spchol(eta0CovL);

			// pre-compute likelihood
			prop.scratch.lly = mcstat2::glm::ll(prop.data.Y->data(),
				prop.params.eta0.data(), cfg->params.beta.data(), prop.data.X->data(),
				prop.consts.n, prop.consts.nt, prop.consts.p,
				mcstat2::glm::glmfamily::poisson);

			// compute forward transition probability for eta0
			fwdll = mcstat2::ldmvrnorm_spchol(prop.params.eta0, eta0mu, eta0CovL);

			// compute backward transition probability for eta0
			gaussian_approx_eta0(prop.params.eta0.data(), cfg->scratch.eta0Prec, 1,
				eta0mu.data(), eta0CovL, cfg->params.beta.data(), prop.data.Y->data(),
				prop.data.X->data(), prop.consts.n, prop.consts.nt, prop.consts.p,
				mcstat2::glm::glmfamily::poisson);
			backll = mcstat2::ldmvrnorm_spchol(cfg->params.eta0, eta0mu, eta0CovL);
		}

		void update() {
			cfg->params.kappa = prop.params.kappa;
			cfg->params.sigmasq = prop.params.sigmasq;
			cfg->params.rho = prop.params.rho;
			cfg->params.eta0 = prop.params.eta0;

			cfg->scratch.lly = prop.scratch.lly;
			cfg->scratch.lleta0 = prop.scratch.lleta0;
			cfg->scratch.lltheta = prop.scratch.lltheta;
			cfg->scratch.eta0Prec = prop.scratch.eta0Prec;
		}
};


//
// Rcpp exports
//

// [[Rcpp::export]]
List respglm_fit(arma::mat& dknots, arma::mat& dzknots, Eigen::MatrixXd& W,
	int nSamples, List priors, std::vector<double>& inits,
	std::vector<double>& sds, std::vector<double>& C,
	Eigen::SparseMatrix<double>& Q, Eigen::VectorXd inits_eta0,
	Eigen::VectorXd inits_beta, Eigen::VectorXd& Y, Eigen::MatrixXd& X,
  Eigen::MatrixXd& A, int df) {
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
				df - rank deficiency in Q
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
	cfg.consts.df = df;

	cfg.priors.rho_L = priors["rho_L"];
	cfg.priors.rho_U = priors["rho_U"];
	cfg.priors.kappa_a = priors["kappa_a"];
	cfg.priors.kappa_b = priors["kappa_b"];
	cfg.priors.s_a = priors["sigmasq_a"];
	cfg.priors.s_b = priors["sigmasq_b"];

	cfg.params.nu = as<double>(priors["nu"]);
	cfg.params.kappa = inits[0];
	cfg.params.rho = inits[1];
	cfg.params.sigmasq = inits[2];
	cfg.params.beta = inits_beta;
	cfg.params.eta0 = inits_eta0;

	// initialize scratch

	cfg.scratch.lly = - 1e6;
	cfg.scratch.lleta0 = - 1e6;
	cfg.scratch.llbeta = - 1e6;
	cfg.scratch.lltheta = - 1e6;

	// initialize samplers

	bool adapt = *(std::max_element(C.begin(), C.end())) != 0 ? true : false;
	TeleSampler ts = TeleSampler(cfg, sds, inits, C, .4,
		{0, 0, cfg.priors.rho_L}, {0, 1, cfg.priors.rho_U}, adapt);

  StateSampler ss = StateSampler(cfg);

	mcstat2::GibbsSampler sampler = mcstat2::GibbsSampler();
	sampler.addSampler(ts);
	sampler.addSampler(ss);
	sampler.run(nSamples);

	Rcpp::Rcout << "kappa: " << ts.getSd(0) << " sigmasq: " << ts.getSd(1) <<
		" rho: " << ts.getSd(2) << std::endl;

	return sampler.getSamples();
}
