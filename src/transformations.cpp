#include "transformations.h"

double mcstat2::logit(double x) { return log( x / (1.0 - x) ); }

double mcstat2::invlogit(double x) {
	double expX = exp(x);
	return std::isinf(expX)!=0 ? 1 : expX / (1.0 + expX);
}

double mcstat2::logitProposal(double x, double min_x, double max_x, double sd) {
	double w = max_x - min_x;
	return mcstat2::invlogit(
		mcstat2::logit((x - min_x)/w) + R::rnorm(0, sd) ) * w + min_x;
}

double mcstat2::logProposal(double x, double sd) {
	return exp( log(x) + R::rnorm(0, sd) );
}

double mcstat2::loglogJacobian(double x) { return -log( std::abs(x) ); }

double mcstat2::loglogitJacobian(double x) {
	return -log( std::abs(x*(1.0-x)) );
}
