// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <algorithm>

using namespace Rcpp;
using std::endl;


// [[Rcpp::export]]
DataFrame online_fallback_faster(NumericVector pval,
	NumericVector gammai,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();
	NumericVector alphai(N);
	LogicalVector R(N);

	alphai[0] = alpha * gammai[0];
	R[0] = (pval[0] <= alphai[0]);

	Progress p(N, display_progress);

	for(int i = 1; i < N; i++) {
		p.increment();
		alphai[i] = alpha * gammai[i] + R[i-1] * alphai[i-1];
		R[i] = (pval[i] <= alphai[i]);
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}
