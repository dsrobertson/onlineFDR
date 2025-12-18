// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <algorithm>

using namespace Rcpp;
using std::endl;

// Debugging function to print contents of a vector.
// Need separate functions for each datatype.
// void printVec(NumericVector vec) {
// 	for (int i = 0; i < vec.size(); i++)
// 		Rcout << vec[i] << " ";
// 	Rcout << endl;
// }
// void printVec(IntegerVector vec) {
// 	for (int i = 0; i < vec.size(); i++)
// 		Rcout << vec[i] << " ";
// 	Rcout << endl;
// }

// [[Rcpp::export]]
DataFrame lond_faster(NumericVector pval,
	NumericVector betai,
	double alpha = 0.05,
	bool original = true,
	bool display_progress = true) {
	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);

	alphai[0] = betai[0];
	R[0] = (pval[0] <= alphai[0]);

	int D = R[0];

	Progress p(N, display_progress);

	if (original == 0){
		for(int i = 1; i < N; i++) {
			p.increment();
			alphai[i] = betai[i]*std::max(D,1);
			if (pval[i] <= alphai[i]) {
				R[i] = 1;
				D++;
			}
		}
	} else {
		//original LOND
		for(int i = 1; i < N; i++) {
			p.increment();
			alphai[i] = betai[i]*(D+1);
			if (pval[i] <= alphai[i]) {
				R[i] = 1;
				D++;
			}
		}
	}
	
	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);

	}
