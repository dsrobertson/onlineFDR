#include <algorithm>
#include <Rcpp.h>

using namespace Rcpp;
using std::endl;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
DataFrame londstar_async_faster(NumericVector pval,
	IntegerVector E,
	NumericVector betai,
	double alpha = 0.05) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	alphai[0] = betai[0];
	R[0] = (pval[0] <= alphai[0]);

	for (int i = 1; i < N; i++) {
		int Dsum = 0;
		for (int j = 0; j <= i-1; j++) {
			if (R[j] && (E[j]-1 <= i-1))
				Dsum++;
			// Rcout << "Dsum: " << Dsum << " R: " << R[j] << " E: " << E[j] << " j: " << j << " i-1: " << i-1 << endl;

		}
		int D = std::max(Dsum, 1);
		alphai[i] = betai[i] * D;
		R[i] = (pval[i] <= alphai[i]);
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}

// [[Rcpp::export]]
DataFrame londstar_dep_faster(NumericVector pval,
	IntegerVector L,
	NumericVector betai,
	double alpha = 0.05) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	alphai[0] = betai[0];
	R[0] = (pval[0] <= alphai[0]);

	for (int i = 1; i < N; i++) {
		int Dsum = 0;
		for (int j = 0; j <= i-1; j++) {
			if (R[j] && (j <= j - L[i]))
				Dsum++;
		}
		int D = std::max(Dsum, 1);
		alphai[i] = betai[i] * D;
		R[i] = (pval[i] <= alphai[i]);
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
	double alpha = 0.05) {


	NumericMatrix alphai(batch.size(), max(batch));
	LogicalMatrix R(batch.size(), max(batch));

	for (int i = 0; i < batch[0]; i++) {
		alphai(0,i) = betai[i];
		R(0,i) = (pval[i] <= alphai(0,i));
	}

	for (int b = 1; b < batch.size(); b++) {
		int Dsum = 0;
		for (int j = 0; j <= b-1; j++) {
			for (int k = 0; k <= R.ncol()-1; k++) {
				if(R(j,k))
					Dsum++;
			}
		}
		int D = std::max(Dsum, 1);
		for (int x = 0; x < batch[b]; x++) {
			alphai(b,x) = betai[batchsum[b-1] + x] * D;
			R(b,x) = (pval[batchsum[b-1] + x] <= alphai(b,x));
		}
	}

	return List::create(
		_["alphai"] = alphai,
		_["R"] = R);
}