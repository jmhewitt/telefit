#include "stpnotation.h"


List CompositionSamples::toSummarizedList() {
	
	// post process teleconnection knots
	
	mat alpha_knots_est, alpha_knots_sd;
	List alpha_knots_sum;
	
	if(!localOnly) {
		alpha_knots_est = mean(alpha_knots);
		alpha_knots_sd = stddev(alpha_knots, 1);
		
		alpha_knots_sum = List::create(
									   _["est"] = alpha_knots_est,
									   _["sd"] = alpha_knots_sd,
									   _["nSamples"] = alpha_knots.n_rows
									   );
	}
	
	// post process full teleconnection field
	
	mat alpha_est, alpha_sd;
	
	List alpha_sum;
	
	if(!localOnly) {
		if(return_full_alpha) {
			alpha_est = mean(alpha);
			alpha_sd = stddev(alpha, 1);
		
			alpha_sum = List::create(
									 _["est"] = alpha_est,
									 _["sd"] = alpha_sd,
									 _["nSamples"] = alpha.n_rows
									 );
		}
	}
	
	// package forecast results
	
	List forecast_sum;
 
	if(localOnly) {
		forecast_sum = List::create(
									_["forecast"] = forecast
									);
	} else {
		forecast_sum = List::create(
									_["forecast"] = forecast,
									_["local"] = local,
									_["remote"] = remote
									);
	}
	
	
	// return results
	List ret;
	
	if(return_full_alpha & return_forecast) {
		ret = List::create(
			_["alpha_knots"] = alpha_knots_sum,
			_["alpha"] = alpha_sum,
			_["forecast"] = forecast_sum
		);
	} else if(return_full_alpha) {
		ret = List::create(
			_["alpha_knots"] = alpha_knots_sum,
			_["alpha"] = alpha_sum
		);
	} else if(return_forecast & (!localOnly)) {
		ret = List::create(
			_["alpha_knots"] = alpha_knots_sum,
			_["forecast"] = forecast_sum
		);
	} else if(return_forecast & localOnly) {
		ret = List::create(
							_["forecast"] = forecast_sum
							);
	}
	
	return ret;
}
