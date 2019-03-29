#ifndef _TRANSFORMATIONS_H
#define _TRANSFORMATIONS_H

#include <RcppArmadillo.h>


namespace mcstat2 {

	// transformation functions
	double logit(double x);
	double invlogit(double x);

	// jacobians for transformation functions
	double loglogJacobian(double x);
	double loglogitJacobian(double x);

	// proposal functions for RW samplers
	double logitProposal(double x, double min_x, double max_x, double sd);
	double logProposal(double x, double sd);

}


#endif
