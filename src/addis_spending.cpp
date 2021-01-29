// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>#include <algorithm>

using namespace Rcpp;
using std::endl;

// [[Rcpp::export]]
DataFrame addis_spending_faster(NumericVector pval,
	NumericVector gammai = NumericVector(0),
	double alpha = 0.05,
	double lambda = 0.25,
	double tau = 0.5,
	bool display_progress = true) {
	
	int N = pval.size();
	NumericVector alphai(N);
	LogicalVector R(N);

	alphai[0] = alpha * (tau - lambda) * gammai[0];
	R[0] = (pval[0] <= alphai[0]);

	int selectsum = (pval[0] <= tau);
	int candsum = (pval[0] <= lambda);

	Progress p(N, display_progress);

	for (int i = 1; i < N; i++) {
		p.increment();
		alphai[i] = alpha * (tau - lambda) * gammai[selectsum - candsum];
		R[i] = (pval[i] <= alphai[i]);
		selectsum = selectsum + (pval[i] <= tau);
		candsum = candsum + (pval[i] <= lambda);
	}
	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}

// [[Rcpp::export]]
DataFrame addis_spending_dep_faster(NumericVector pval,
	IntegerVector L,
	NumericVector gammai = NumericVector(0),
	double alpha = 0.05,
	double lambda = 0.25,
	double tau = 0.5,
	bool display_progress = true) {

	int N = pval.size();
	NumericVector alphai(N);
	LogicalVector R(N);
	LogicalVector select(N);
	LogicalVector cand(N);

	alphai[0] = alpha * (tau - lambda) * gammai[0];
	R[0] = (pval[0] <= alphai[0]);
	select[0] = (pval[0] <= tau);
	cand[0] = (pval[0] <= lambda);

	Progress p(N, display_progress);

	for (int i = 1; i < N; i++) {
		p.increment();
		int selectsum = 0;
		int candsum = 0;
		int maxL = std::max(0, i - L[i]);
		if (maxL > 0) {
			for (int j = 0; j <= maxL; j++) {
				if (select[j])
					selectsum++;
				if (cand[j])
					candsum++;
			}
		}

		alphai[i] = alpha * (tau - lambda) * gammai[1 + std::min(L[i]-1, i-1) + selectsum - candsum];
		R[i] = (pval[i] <= alphai[i]);
		select[i] = (pval[i] <= tau);
		cand[i] = (pval[i] <= lambda);
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}