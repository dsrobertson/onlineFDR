// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <algorithm>
#include <vector>

using namespace Rcpp;
using std::endl;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
DataFrame saffronstar_async_faster(NumericVector pval,
	IntegerVector E,
	NumericVector gammai,
	double w0 = 0.0125,
	double lambda = 0.5,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	IntegerVector cand(N);
	IntegerVector Cjplus(N);
	alphai[0] = std::min(gammai[0]*w0, lambda);
	R[0] = (pval[0] <= alphai[0]);
	int candsum = 0;

	std::vector<bool> r;

	Progress p(N * N, display_progress);
	for (int i = 1; i < N; i++) {

		cand[i-1] = (pval[i-1] <= lambda);

		r.clear();
		for (int j = 0; j <= i-1; j++) {
			p.increment();
			if (R[j] && (E[j]-1 <= i-1))
				r.push_back(j);
		}

		candsum += (int)(cand[i-1] && (E[i-1]-1 <= i-1));

		int K = r.size();

		double alphaitilde;
		if (K > 1) {
			
	    //update Cjplus
			double Cjplussum = 0;
			for (int j = 0; j < K; j++) {

				int from = r[j]+1;
				int to = std::max(i-1, (int)(r[j]+1));
				int sum = 0;

				for (int k = from; k <= to; k++) {
					if (cand[k] && E[k]-1 <= i-1)
						sum++;
				}

				Cjplus[j] = sum;
				Cjplussum += gammai[i-r[j]-Cjplus[j]-1];
			}
			
			Cjplussum -= gammai[i-r[0]-Cjplus[0]-1];

			alphaitilde = w0*gammai[i-candsum] + ((1-lambda)*alpha - w0)*
			gammai[i-r[0]-Cjplus[0]-1] + (1-lambda)*alpha*Cjplussum;
			
		} else if (K == 1) {
			
			int from = r[0]+1;
			int to = std::max(i-1, (int)(r[0]+1));
			Cjplus[0] = 0;
			for (int j = from; j <= to; j++) {
				if (cand[j] && E[j]-1 <= i-1)
					Cjplus[0]++;
			}
			alphaitilde = w0*gammai[i-candsum] + ((1-lambda)*alpha-w0)*
			gammai[i-r[0]-Cjplus[0]-1];
			
		} else {
			alphaitilde = w0*gammai[i-candsum];
		}
		alphai[i] = std::min(lambda, alphaitilde);
		if (pval[i] <= alphai[i]) {
			R[i] = 1;
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
	double w0 = 0.0125,
	double lambda = 0.5,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	IntegerVector cand(N);
	IntegerVector Cjplus(N);
	alphai[0] = std::min(gammai[0]*w0, lambda);
	R[0] = (pval[0] <= alphai[0]);

	std::vector<bool> r;

	Progress p(N * N, display_progress);

	for (int i = 1; i < N; i++) {

		cand[i-1] = (pval[i-1] <= lambda);
		r.clear();
		for (int j = 0; j <= i-1; j++) {
			p.increment();
			if (R[j] && (j <= j - L[i]))
				r.push_back(j);
		}

		int bound = i-1-L[i];
		int candsum = 0;
		for (int m = 0; m <= bound; m++) {
				candsum += cand[m];
		}

		int K = r.size();

		double alphaitilde;
		if (K > 1) {
			
	    //update Cjplus
			double Cjplussum = 0;
			for (int j = 0; j < K; j++) {

				int from = r[j]+1;
				int to = std::max(i-1, (int)(r[j]+1));
				int sum = 0;

				for (int k = from; k <= to; k++) {
					if (cand[k] && k <= i-1-L[i])
						sum++;
				}

				Cjplus[j] = sum;
				Cjplussum += gammai[i-r[j]-Cjplus[j]-1];
			}

			Cjplussum -= gammai[i-r[0]-Cjplus[0]-1];

			alphaitilde = w0*gammai[i-candsum] + ((1-lambda)*alpha-w0)*
			gammai[i-r[0]-Cjplus[0]-1] + (1-lambda)*alpha*Cjplussum;
			
		} else if (K == 1) {
			
			int from = r[0]+1;
			int to = std::max(i-1, (int)(r[0]+1));
			Cjplus[0] = 0;
			for (int j = from; j <= to; j++) {
				if (cand[j] && j <= i-L[i]-1)
					Cjplus[0]++;
			}
			alphaitilde = w0*gammai[i-candsum] + ((1-lambda)*alpha-w0)*
			gammai[i-r[0]-Cjplus[0]-1];
			
		} else {
			alphaitilde = w0*gammai[i-candsum];
		}
		alphai[i] = std::min(lambda, alphaitilde);
		if (pval[i] <= alphai[i]) {
			R[i] = 1;
		}
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
	double w0 = 0.0125,
	double lambda = 0.5,
	double alpha = 0.05,
	bool display_progress = true) {

	int N = pval.size();

	int B = batch.size();
	NumericMatrix alphai(B, max(batch));
	LogicalMatrix R(B, max(batch));
	IntegerVector cand(N);
	IntegerVector Cj(B);

	int mysum = 0;
	for (int a = 1; a < batch.size(); a++) {
		mysum += batch[a];
	}

	Progress p(mysum, display_progress);

	for (int i = 0; i < batch[0]; i++) {
		cand[i] = (pval[i] <= lambda);
		alphai(0,i) = gammai[i] * w0;
		R(0,i) = (pval[i] <= alphai(0,i));
	}

	IntegerVector Cjplus(B);
	Cj[0] = sum(cand);

	for (int b = 1; b < B; b++) {
		NumericVector rcum = cumsum(static_cast<NumericVector>(rowSums(R)));
		int candsum = sum(Cj);

		for (int x = 0; x < batch[b]; x++) {
			cand[batchsum[b-1] + x] = (pval[batchsum[b-1] + x] <= lambda);

			p.increment();

			NumericVector r(0);
			if (max(rcum) > 0) {
				for (int y = 0; y < max(rcum); y++) {
					if (rcum[y] >= y)
						r.push_back(y);
				}
			}

			int K = r.size();
			double alphaitilde;
			if (K > 1) {

	    //update Cjplus
				double Cjplussum = 0;
				for (int j = 0; j < K; j++) {

					int from = r[j]+1;
					int to = b-1;
					int sum = 0;

					for (int k = from; k <= to; k++) {
						if (Cj[k])
							sum++;
					}

					Cjplus[j] = sum * (int)(b-1 >= from);
					Cjplussum += gammai[batchsum[b-1] + x - batchsum[r[j]] - Cjplus[j]];
				}

				Cjplussum -= gammai[batchsum[b-1] + x - batchsum[r[0]] - Cjplus[0]];

				alphaitilde = w0*gammai[batchsum[b-1] + x - candsum] + ((1-lambda)*alpha - w0)*
				gammai[batchsum[b-1] + x - batchsum[r[0]] - Cjplus[0]] + (1-lambda)*alpha*Cjplussum;
				alphai(b,x) = std::min(lambda, alphaitilde);
				R(b,x) = (pval[batchsum[b-1] + x] <= alphai(b,x));

			} else if (K == 1) {

				int from = r[0]+1;
				int to = b-1;
				int sum = 0;
				Cjplus[0] = 0;
				for (int j = from; j <= to; j++) {
					if (Cj[j])
						sum++;
					Cjplus[0] = sum*(int)(b-1 >= from);
				}
				alphaitilde = w0*gammai[batchsum[b-1] + x - candsum] + ((1-lambda)*alpha-w0)*
				gammai[batchsum[b-1] + x - batchsum[r[0]] - Cjplus[0]];
				alphai(b,x) = std::min(lambda, alphaitilde);
				R(b,x) = (pval[batchsum[b-1] + x] <= alphai(b,x));

			} else {
				alphaitilde = w0*gammai[batchsum[b-1] + x - candsum];
				alphai(b,x) = std::min(lambda, alphaitilde);
				R(b,x) = (pval[batchsum[b-1] + x] <= alphai(b,x));
			}
		}
		int from = batchsum[b-1] + 1;
		int to = batchsum[b];
		int sum = 0;
		for (int z = from; z <= to; z++) {
			if(cand[z-1])
				sum++;
		}
		Cj[b] = sum;

	}

	return List::create(
		_["alphai"] = alphai,
		_["R"] = R);
}
