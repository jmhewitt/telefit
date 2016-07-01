#ifndef _telefit_COVS_H
#define _telefit_COVS_H


#include <RcppArmadillo.h>

using namespace arma;


void maternCov( mat & cov, const mat & d, double scale, double range,
				double smoothness, double nugget );

void maternArray( vec & d, double scale, double range,
			      double smoothness, double nugget );

#endif