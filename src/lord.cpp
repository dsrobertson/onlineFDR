#include <Rcpp.h>
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
DataFrame lord_faster(NumericVector pval,
	NumericVector gammai,
	int version,
	double alpha = 0.05,
	double w0 = 0.005,
	double b0 = 0.045,
	double taudiscard = 0.5) {

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R;
	NumericVector W;

	// ++

	if (version == 1) {
		R = LogicalVector(N);
		alphai[0] = gammai[0]*w0;
		R[0] = (pval[0] <= alphai[0]);

		IntegerVector tau(1);
		int K = sum(R);

		for (int i = 1; i < N; i++) {

			if (K <= 1) {

				if (R[i-1])
					tau[0] = i-1;

				double Cjsum = 0;
				for(int j = 0; j < K; j++){
					Cjsum += gammai[ i-tau[j]-1 ];

				}
				alphai[i] = w0*gammai[i] + (alpha-w0)*Cjsum;

			} else {

				if (R[i-1])
					tau.push_back(i-1);

				double Cjsum = 0;
				IntegerVector tau2 = clone(tau);
				tau2.erase(tau2.begin());

				for (int j = 0; j < K-1; j++){
					Cjsum += gammai[ i-tau2[j]-1 ];

				}
				alphai[i] = w0*gammai[i] + (alpha-w0)*gammai[ i-tau[0]-1 ] + alpha*Cjsum;
				
			}

			if (pval[i] <= alphai[i]) {
				R[i] = 1;
				K++;
			}
		}
	}

	// discard

	if (version == 2) {
		R = LogicalVector(N);
		alphai[0] = gammai[0]*w0;
		R[0] = (pval[0] <= alphai[0]);

		LogicalVector selected = (pval <= taudiscard);
		NumericVector S = cumsum(static_cast<NumericVector>(selected));

		IntegerVector kappai(1);
		int K = sum(R);

		for (int i = 1; i < N; i++) {

			double alphaitilde;

			if (K > 1){

				if(R[i-1])
					kappai.push_back(i-1);

				IntegerVector kappaistar(kappai.size());

				int mysum = 0;
				int index = 0;
				int bound = kappai[kappai.size()-1];
				for (int k = 0; k <= bound; k++) {
					mysum += selected[k];
				//this is the sapply workaround
					if (kappai[index] == k){
						kappaistar[index] = mysum;
						index++;
					}
				}

				double Cjsum = 0;
				IntegerVector kappaistar2 = clone(kappaistar);
				kappaistar2.erase(kappaistar2.begin());
				for(int j = 0; j < K-1; j++){
					Cjsum += gammai[ S[i-1]-kappaistar2[j] ];
				}
				
				alphaitilde = w0*gammai[ S[i-1] ] + 
				(taudiscard*alpha - w0)*gammai[ S[i-1]-kappaistar[0] ] +
				taudiscard*alpha*Cjsum;

			} else if (K == 1) {

				if(R[i-1])
					kappai[0] = i-1;

				int kappaistar = 0;
				for (int j = 0; j < kappai[0]; j++)
					kappaistar += selected[j];

				alphaitilde = w0*gammai[ S[i-1] ] +
				(taudiscard*alpha - w0)*gammai[ S[i-1] - kappaistar - 1 ];

			} else {

				alphaitilde = w0*gammai[ S[i-1] ];

			}

			alphai[i] = std::min(taudiscard, alphaitilde);
			if (pval[i] <= alphai[i]) {
				R[i] = 1;
				K++;
			}
		}
	}

	// 3

	if (version == 3) {
		R = LogicalVector(N+1);
		W = NumericVector(N+1);
		R[0] = 1;
		W[0] = w0;
		alphai[0] = gammai[0]*w0;
		double phi = gammai[0]*w0;
		R[1] = (pval[0] <= alphai[0]);
		W[1] = w0-phi+R[1]*b0;

		IntegerVector tau(1);

		for (int i = 1; i < N; i++) {
			if(R[i])
				tau.push_back(i);
			double taumax = max(tau);
			alphai[i] = gammai[ i-taumax ]*W[taumax];
			phi = gammai[ i-taumax ]*W[taumax];

			R[i+1] = (pval[i] <= alphai[i]);
			W[i+1] = W[i] - phi + R[i-1]*b0;
		}

		R.erase(R.begin());
	}

	// dep
	
	if (version == 4) {
		R = LogicalVector(N+1);
		W = NumericVector(N+1);
		R[0] = 1;
		W[0] = w0;
		alphai[0] = gammai[0]*w0;
		double phi = gammai[0]*w0;
		R[1] = (pval[0] <= alphai[0]);
		W[1] = w0-phi+R[1]*b0;

		IntegerVector tau(1);

		for (int i = 1; i < N; i++) {
			if(R[i])
				tau.push_back(i);
			double taumax = max(tau);
			alphai[i] = gammai[i]*W[taumax];
			phi = gammai[i]*W[taumax];

			R[i+1] = (pval[i] <= alphai[i]);
			W[i+1] = W[i] - phi + R[i+1]*b0;
		}

		R.erase(R.begin());
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}