// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <algorithm>

using namespace Rcpp;
using std::endl;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
DataFrame londstar_async_faster(NumericVector pval,
	IntegerVector E,
	NumericVector betai,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	alphai(0) = betai(0);
	R(0) = (pval(0) <= alphai(0));

	Progress p(N * N, display_progress);

	for (int i = 1; i < N; i++) {
		int Dsum = 0;
		for (int j = 0; j <= i-1; j++) {
			p.increment();
			if (R(j) && (E(j)-1 <= i-1))
				Dsum++;

		}
		int D = std::max(Dsum, 1);
		alphai(i) = betai(i) * D;
		R(i) = (pval(i) <= alphai(i));
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}

// [[Rcpp::export]]
DataFrame londstar_dep_faster(NumericVector pval,
	IntegerVector L,
	NumericVector betai,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	alphai(0) = betai(0);
	R(0) = (pval(0) <= alphai(0));

	Progress p(N * N, display_progress);

	for (int i = 1; i < N; i++) {
		int Dsum = 0;
		for (int j = 0; j <= i-1; j++) {
			p.increment();
			if (R(j) && (j < i - L(i)))
				Dsum++;
		}
		int D = std::max(Dsum, 1);
		alphai(i) = betai(i) * D;
		R(i) = (pval(i) <= alphai(i));
	}

	return DataFrame::create(_["pval"] = pval,
		_["lag"] = L,
		_["alphai"] = alphai,
		_["R"] = R);

}

// [[Rcpp::export]]
List londstar_batch_faster(NumericVector pval,
	IntegerVector batch,
	IntegerVector batchsum,
	NumericVector betai,
	double alpha = 0.05,
	bool display_progress = true) {

	int B = batch.size();

	NumericMatrix alphai(B, max(batch));
	LogicalMatrix R(B, max(batch));

	for (int i = 0; i < batch(0); i++) {
		alphai(0,i) = betai(i);
		R(0,i) = (pval(i) <= alphai(0,i));
	}

	int mysum = 0;
	for (int a = 1; a < batch.size(); a++) {
		mysum += batch(a);
	}

	Progress p(mysum, display_progress);

	for (int b = 1; b < B; b++) {
		int Dsum = 0;
		for (int j = 0; j <= b-1; j++) {

			for (int k = 0; k <= R.ncol()-1; k++) {

				if(R(j,k))
					Dsum++;
			}
		}
		int D = std::max(Dsum, 1);
		for (int x = 0; x < batch(b); x++) {
			p.increment();
			alphai(b,x) = betai(batchsum(b-1) + x) * D;
			R(b,x) = (pval(batchsum(b-1) + x) <= alphai(b,x));
		}
	}

	return List::create(
		_["alphai"] = alphai,
		_["R"] = R);
}
