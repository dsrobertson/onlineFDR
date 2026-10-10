// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// In each version, r holds one entry per rejection whose decision is known,
// in order: the index of the step at which it became known (step i-1 for
// test i), or for the batch version the batch it was made in. The number of
// known rejections never decreases, so r is extended rather than rebuilt.
// Cjplus[j] is the number of known candidates (p-values <= lambda) after the
// j-th entry of r, found from running counts rather than by scanning.

// [[Rcpp::export]]
DataFrame saffronstar_async_faster(NumericVector pval,
	IntegerVector E,
	NumericVector gammai,
	double w0 = 0.025,
	double lambda = 0.5,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	std::vector<int> r;
	// newR[t]: number of rejections j whose decision is first known at test
	// t, i.e. j < t and E[j] <= t; newCand[t]: the candidates first known at t
	std::vector<int> newR(N + 1, 0);
	std::vector< std::vector<int> > newCand(N + 1);
	// Fenwick tree over test indices, counting the known candidates
	std::vector<int> tree(N + 1, 0);
	alphai(0) = std::min((1-lambda)*gammai(0)*w0, lambda);
	R(0) = (pval(0) <= alphai(0));
	if (R(0)) {
		newR[std::min(N, std::max(1, (int)E(0)))]++;
	}

	Progress p(N, display_progress);

	int cond = 0;
	int candsum = 0;
	for (int i = 1; i < N; i++) {
		p.increment();

		if (pval(i-1) <= lambda) {
			newCand[std::min(N, std::max(i, (int)E(i-1)))].push_back(i-1);
		}
		for (size_t c = 0; c < newCand[i].size(); c++) {
			for (int k = newCand[i][c] + 1; k <= N; k += k & (-k)) {
				tree[k]++;
			}
			candsum++;
		}

		int prev = cond;
		cond += newR[i];
		for (int y = prev; y < cond; y++) {
			r.push_back(i - 1);
		}

		int K = r.size();

		// number of known candidates with index greater than s
		auto after = [&](int s) {
			int n = 0;
			for (int k = s + 1; k > 0; k -= k & (-k)) {
				n += tree[k];
			}
			return candsum - n;
		};

		double alphaitilde;
		if (K > 1) {

			double Cjplussum = 0;
			for (int j = 0; j < K; j++) {
				Cjplussum += gammai(i-r[j]-after(r[j])-1);
			}

			int Cjplus0 = after(r[0]);
			Cjplussum -= gammai(i-r[0]-Cjplus0-1);

			alphaitilde = (1-lambda)*(w0*gammai(i-candsum) + (alpha - w0)*
			gammai(i-r[0]-Cjplus0-1) + alpha*Cjplussum);

		} else if (K == 1) {

			int Cjplus0 = after(r[0]);
			alphaitilde = (1-lambda)*(w0*gammai(i-candsum) + (alpha-w0)*
			gammai(i-r[0]-Cjplus0-1));

		} else {
			alphaitilde = (1-lambda)*w0*gammai(i-candsum);
		}
		alphai(i) = std::min(lambda, alphaitilde);
		if (pval(i) <= alphai(i)) {
			R(i) = 1;
			newR[std::min(N, std::max(i + 1, (int)E(i)))]++;
		}
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}

// [[Rcpp::export]]
DataFrame saffronstar_dep_faster(NumericVector pval,
	IntegerVector L,
	NumericVector gammai,
	double w0 = 0.025,
	double lambda = 0.5,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	std::vector<int> r;
	// cumR[k], cumC[k]: number of rejections and of candidates among the first
	// k tests
	std::vector<int> cumR(N + 1, 0);
	std::vector<int> cumC(N + 1, 0);
	alphai(0) = std::min((1-lambda)*gammai(0)*w0, lambda);
	R(0) = (pval(0) <= alphai(0));
	cumR[1] = R(0);

	Progress p(N, display_progress);

	// The lags satisfy L[i+1] <= L[i] + 1, so the number of known rejections,
	// cond, never decreases
	int cond = 0;
	for (int i = 1; i < N; i++) {
		p.increment();

		cumC[i] = cumC[i-1] + (pval(i-1) <= lambda);

		// the decisions of tests 0, ..., hi-1 are known
		int hi = i - L(i);

		int prev = cond;
		cond = (hi >= 1) ? cumR[hi] : 0;
		for (int y = prev; y < cond; y++) {
			r.push_back(i - 1);
		}

		int candsum = (hi >= 1) ? cumC[hi] : 0;

		int K = r.size();

		// number of known candidates with index greater than s
		auto after = [&](int s) {
			return (hi > s + 1) ? cumC[hi] - cumC[s + 1] : 0;
		};

		double alphaitilde;
		if (K > 1) {

			double Cjplussum = 0;
			for (int j = 0; j < K; j++) {
				Cjplussum += gammai(i-r[j]-after(r[j])-1);
			}

			int Cjplus0 = after(r[0]);
			Cjplussum -= gammai(i-r[0]-Cjplus0-1);

			alphaitilde = (1-lambda)*(w0*gammai(i-candsum) + (alpha-w0)*
			gammai(i-r[0]-Cjplus0-1) + alpha*Cjplussum);

		} else if (K == 1) {

			int Cjplus0 = after(r[0]);
			alphaitilde = (1-lambda)*(w0*gammai(i-candsum) + (alpha-w0)*
			gammai(i-r[0]-Cjplus0-1));

		} else {
			alphaitilde = (1-lambda)*w0*gammai(i-candsum);
		}
		alphai(i) = std::min(lambda, alphaitilde);
		if (pval(i) <= alphai(i)) {
			R(i) = 1;
		}
		cumR[i + 1] = cumR[i] + R(i);
	}

	return DataFrame::create(_["pval"] = pval,
		_["lag"] = L,
		_["alphai"] = alphai,
		_["R"] = R);
}

// [[Rcpp::export]]
List saffronstar_batch_faster(NumericVector pval,
	IntegerVector batch,
	IntegerVector batchsum,
	NumericVector gammai,
	double w0 = 0.025,
	double lambda = 0.5,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();
	int B = batch.size();

	NumericMatrix alphai(B, max(batch));
	LogicalMatrix R(B, max(batch));
	IntegerVector cand(N);
	std::vector<int> r;
	// cumCj[k]: number of candidates in batches 0, ..., k-1
	std::vector<int> cumCj(B + 1, 0);

	int mysum = 0;
	for (int a = 1; a < batch.size(); a++) {
		mysum += batch(a);
	}

	Progress p(mysum, display_progress);

	for (int i = 0; i < batch(0); i++) {
		cand(i) = (pval(i) <= lambda);
		alphai(0,i) = std::min((1-lambda)*gammai(i) * w0, lambda);
		R(0,i) = (pval(i) <= alphai(0,i));
	}

	cumCj[1] = sum(cand);

	for (int b = 1; b < B; b++) {

		// add the rejections of the previous batch
		for (int x = 0; x < batch(b-1); x++) {
			if (R(b-1,x)) {
				r.push_back(b - 1);
			}
		}

		int candsum = cumCj[b];
		int K = r.size();
		double alphaitilde;

		// candidates in the batches after the one of each rejection
		std::vector<int> Cjplus(K);
		for (int j = 0; j < K; j++) {
			Cjplus[j] = (r[j] + 1 <= b - 1) ? cumCj[b] - cumCj[r[j] + 1] : 0;
		}

		for (int x = 0; x < batch(b); x++) {
			cand(batchsum(b-1) + x) = (pval(batchsum(b-1) + x) <= lambda);

			p.increment();

			if (K > 1) {

				double Cjplussum = 0;
				for (int j = 0; j < K; j++) {
					Cjplussum += gammai(batchsum(b-1) + x - batchsum(r[j]) - Cjplus[j]);
				}

				Cjplussum -= gammai(batchsum(b-1) + x - batchsum(r[0]) - Cjplus[0]);

				alphaitilde = (1-lambda)*(w0*gammai(batchsum(b-1) + x - candsum) +
				  (alpha - w0)*gammai(batchsum(b-1) + x - batchsum(r[0]) - Cjplus[0]) +
				  alpha*Cjplussum);

			} else if (K == 1) {

				alphaitilde = (1-lambda)*(w0*gammai(batchsum(b-1) + x - candsum) +
				  (alpha-w0)*gammai(batchsum(b-1) + x - batchsum(r[0]) - Cjplus[0]));

			} else {
				alphaitilde = (1-lambda)*w0*gammai(batchsum(b-1) + x - candsum);
			}
			alphai(b,x) = std::min(lambda, alphaitilde);
			R(b,x) = (pval(batchsum(b-1) + x) <= alphai(b,x));
		}

		int sum = 0;
		for (int z = batchsum(b-1); z < batchsum(b); z++) {
			if (cand(z))
				++sum;
		}

		cumCj[b + 1] = cumCj[b] + sum;
	}

	return List::create(
		_["alphai"] = alphai,
		_["R"] = R);
}
