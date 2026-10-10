// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// In each version, r holds one entry per rejection whose decision is known,
// in order: the index of the step at which it became known (step i-1 for
// test i), or for the batch version the batch it was made in. The number of
// known rejections never decreases, so r is extended rather than rebuilt.

// [[Rcpp::export]]
DataFrame lordstar_async_faster(NumericVector pval,
	IntegerVector E,
	NumericVector gammai,
	double w0 = 0.005,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	std::vector<int> r;
	// newR[t]: number of rejections whose decision is first known at test t,
	// i.e. rejections j with j < t and E[j] <= t
	std::vector<int> newR(N + 1, 0);
	alphai[0] = gammai[0] * w0;
	R[0] = (pval[0] <= alphai[0]);
	if (R[0]) {
		newR[std::min(N, std::max(1, (int)E[0]))]++;
	}

	Progress p(N, display_progress);

	int cond = 0;
	for (int i = 1; i < N; i++) {
		p.increment();

		int prev = cond;
		cond += newR[i];
		for (int y = prev; y < cond; y++) {
			r.push_back(i - 1);
		}

		int K = r.size();

		if (K <= 1) {

			if (K > 0) {

				alphai[i] = gammai[i] * w0 + (alpha - w0) * gammai[i-r[0]-1];

			} else {

				alphai[i] = gammai[i] * w0;
			}

		} else {

			double gammaisum = 0;

			for (int g = 1; g < K; g++) {
				gammaisum += gammai[i-r[g]-1];
			}

			alphai[i] = gammai[i] * w0 + (alpha - w0) * gammai[i-r[0]-1] +
			  alpha * gammaisum;
		}

		R[i] = (pval[i] <= alphai[i]);
		if (R[i]) {
			newR[std::min(N, std::max(i + 1, (int)E[i]))]++;
		}
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);

}

// [[Rcpp::export]]
DataFrame lordstar_dep_faster(NumericVector pval,
	IntegerVector L,
	NumericVector gammai,
	double w0 = 0.005,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	std::vector<int> r;
	// cumR[k]: number of rejections among the first k tests
	std::vector<int> cumR(N + 1, 0);
	alphai[0] = gammai[0] * w0;
	R[0] = (pval[0] <= alphai[0]);
	cumR[1] = R[0];

	Progress p(N, display_progress);

	// The lags satisfy L[i+1] <= L[i] + 1, so the number of known rejections,
	// cond, never decreases
	int cond = 0;
	for (int i = 1; i < N; i++) {
		p.increment();

		int prev = cond;
		cond = (i - L[i] >= 1) ? cumR[i - L[i]] : 0;
		for (int y = prev; y < cond; y++) {
			r.push_back(i - 1);
		}

		int K = r.size();

		if (K <= 1) {

			if (K > 0) {

				alphai[i] = gammai[i] * w0 + (alpha - w0) * gammai[i-r[0]-1];

			} else {

				alphai[i] = gammai[i] * w0;

			}

		} else {

			double gammaisum = 0;

			for (int g = 1; g < K; g++) {
				gammaisum += gammai[i-r[g]-1];
			}

			alphai[i] = gammai[i] * w0 + (alpha - w0) * gammai[i-r[0]-1] +
			  alpha * gammaisum;
		}

		R[i] = (pval[i] <= alphai[i]);
		cumR[i + 1] = cumR[i] + R[i];
	}

	return DataFrame::create(_["pval"] = pval,
		_["lag"] = L,
		_["alphai"] = alphai,
		_["R"] = R);

}

// [[Rcpp::export]]
List lordstar_batch_faster(NumericVector pval,
	IntegerVector batch,
	IntegerVector batchsum,
	NumericVector gammai,
	double w0 = 0.005,
	double alpha = 0.05,
	bool display_progress = true) {

	int B = batch.size();

	NumericMatrix alphai(B, max(batch));
	LogicalMatrix R(B, max(batch));
	std::vector<int> r;

	int mysum = 0;
	for (int a = 1; a < batch.size(); a++) {
		mysum += batch[a];
	}

	Progress p(mysum, display_progress);

	for (int i = 0; i < batch[0]; i++) {
		alphai(0,i) = gammai[i] * w0;
		R(0,i) = (pval[i] <= alphai(0,i));
	}

	for (int b = 1; b < B; b++) {

		// add the rejections of the previous batch
		for (int x = 0; x < batch[b-1]; x++) {
			if (R(b-1,x)) {
				r.push_back(b - 1);
			}
		}

		int K = r.size();

		for (int x = 0; x < batch[b]; x++) {
			p.increment();

			if (K <= 1) {
				if (K > 0) {
					alphai(b,x) = gammai[batchsum[b-1] + x] * w0 + (alpha - w0) *
					gammai[batchsum[b-1] + x - batchsum[r[0]]];

				} else {
					alphai(b,x) = gammai[batchsum[b-1] + x] * w0;
				}
			} else {
				double gammaisum = 0;
				for (int g = 1; g < K; g++) {
					gammaisum += gammai[batchsum[b-1] + x - batchsum[r[g]]];
				}
				alphai(b,x) = gammai[batchsum[b-1] + x] * w0 + (alpha - w0) *
				gammai[batchsum[b-1] + x - batchsum[r[0]]] +
				alpha * gammaisum;
			}
			R(b,x) = (pval[batchsum[b-1] + x] <= alphai(b,x));
		}
	}

	return List::create(
		_["alphai"] = alphai,
		_["R"] = R);
}
