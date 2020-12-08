#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using std::endl;

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
// void printVec(LogicalVector vec) {
// 	for (int i = 0; i < vec.size(); i++)
// 		Rcout << vec[i] << " ";
// 	Rcout << endl;
// }

// [[Rcpp::export]]
DataFrame addis_sync_faster(NumericVector pval,
	NumericVector gammai = NumericVector(0),
	double lambda = 0.5,
	double alpha = 0.05,
	double tau = 0.5,
	double w0 = 0.00625) {
	int N = pval.size();

	if (gammai.size() == 0) {
		gammai = static_cast<NumericVector>(no_init(N));
		for (int i = 0; i < N; i++)
			gammai[i] = 0.4374901658/pow(i+1, 1.6);
	}

	NumericVector alphai(N);
	LogicalVector R(N);
	IntegerVector Cjplus(N);
	IntegerVector cand(N);
	LogicalVector selected = (pval <= tau);
	NumericVector S = cumsum(static_cast<NumericVector>(selected));

	alphai[0] = w0*gammai[0];
	R[0] = (pval[0] <= alphai[0]);

	int K = sum(R);
	int candsum = 0; 
	IntegerVector kappai(1);

	for (int i = 1; i < N; i++) {

		cand[i-1] = (pval[i-1] <= tau*lambda);
		candsum += cand[i-1];

		double alphaitilde;

		if (K > 1) {

			if (R[i-1])
				kappai.push_back(i-1);

			//sapply 
			NumericVector kappaistar(kappai.size());

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

			//update Cjplus
			double Cjplussum = 0;
			for (int j = 0; j < K-1; j++) {
				Cjplus[j] += cand[i-1];
				Cjplussum += gammai[ S[i-1] - kappaistar[j] - Cjplus[j] ];
			}

	    	//update Cjplus again
			Cjplus[K-1] = 0;
			int low = kappai[K-1]+1;
			int high = std::max(i-1, (int)(max(kappai) + 1));
			for (int j = low; j <= high; j++) {
				Cjplus[K-1] += cand[j];
			}

			Cjplussum += gammai[ S[i-1]-kappaistar[K-1]-Cjplus[K-1] ] - 
			gammai[ S[i-1]-kappaistar[0]-Cjplus[0] ];

			alphaitilde = w0*gammai[ S[i-1]-candsum ] + 
			(tau*(1-lambda)*alpha-w0)*gammai[ S[i-1]-kappaistar[0]-Cjplus[0] ] +
			tau*(1-lambda)*alpha*Cjplussum;

		} else if (K == 1) {

			if (R[i-1])
				kappai[0] = i-1;

			int kappaistar = 0;
			for (int j = 0; j <= kappai[0]; j++)
				kappaistar += selected[j];

			Cjplus[0] = 0;
			int low = kappai[0]+1;
			int high = std::max(i-1, kappai[0]+1);
			for (int j = low; j <= high; j++) {
				if (cand[j])
					Cjplus[0]++;
			}

			alphaitilde = w0*gammai[ S[i-1] - candsum  ] + 
			    (tau*(1-lambda)*alpha-w0)*gammai[ S[i-1] - kappaistar - Cjplus[0] ];

		} else {

			alphaitilde = w0*gammai[ S[i-1] - candsum ];

		}

		alphai[i] = std::min(tau*lambda, alphaitilde);
		if (pval[i] <= alphai[i]) {
			R[i] = 1;
			K++;
		}
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}

// [[Rcpp::export]]
DataFrame addis_async_faster(NumericVector pval,
	IntegerVector E,
	NumericVector gammai = NumericVector(0),
	double lambda = 0.5,
	double alpha = 0.05,
	double tau = 0.5,
	double w0 = 0.00625) {

	int N = pval.size();

	if (gammai.size() == 0) {
		gammai = static_cast<NumericVector>(no_init(N));
	// gammai <- 0.4374901658/(seq_len(N)^(1.6))
		for (int i = 0; i < N; i++)
			gammai[i] = 0.4374901658/pow(i+1, 1.6);
	}

	NumericVector alphai(N);
	LogicalVector R(N);
	IntegerVector S(N);
	IntegerVector cand(N);
	IntegerVector Cjplus(N);
	LogicalVector selected = (pval <= tau);

	alphai[0] = w0*gammai[0];
	R[0] = (pval[0] <= alphai[0]);

	int K;
	
	for (int i = 1; i < N; i++) {

		IntegerVector kappai(0);
		// nightmare to code the which statement
		for (int j = 0; j <= i-1; j++) {
			if (R[j] && (E[j]-1 <= i-1))
				kappai.push_back(j);
		}

		K = kappai.size();

		cand[i-1] = (pval[i-1] <= tau*lambda);
		Rcout << "point A" << endl;
		int candsum = 0;
		//C++ trick to loop "seq_len" and use conditional incrementor for "sum"
		for (int j = 0; j <= i-1; j++) {
			if (cand[j] && (E[j]-1 <= i-1))
				candsum++;
		}

		Rcout << "point B" << endl;

		int Ssum = 0;
		for (int j = 0; j <= i-1; j++) {
			if (selected[j] && (E[j]-1 <= i-1))
				Ssum++;
			if (E[j]-1 >= i)
				Ssum++;
		}

		S[i-1] = Ssum;

		double alphaitilde;
		if (K > 1) {

	    //sapply, also nightmare to code 
			NumericVector kappaistar(kappai.size());

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

	    //update Cjplus
			for (int j = 0; j < K; j++) {

				int from = kappai[j]+1;
				int to = std::max(i-1, kappai[j]+1);
				int sum = 0;

				for (int k = from; k <= to; k++) {
					if (cand[k] && E[k]-1 <= i-1)
						sum++;
				}

				Cjplus[j] = sum;
			}

			double Cjplussum = 0;
			//indexing gammai is a nightmare
			for (int j = 0; j < K; j++) {
				Cjplussum += gammai[ S[i-1] - kappaistar[j] - Cjplus[j] ];
			}
			Cjplussum -= gammai[ S[i-1] - kappaistar[0] - Cjplus[0] ];
			
			alphaitilde = w0*gammai[ S[i-1]-candsum ] + 
			(tau*(1-lambda)*alpha-w0)*gammai[ S[i-1]-kappaistar[0]-Cjplus[0] ] +
			tau*(1-lambda)*alpha*Cjplussum;
			
		}  else if (K == 1) {

			int kappaistar = 0;
			for (int j = 0; j <= kappai[0]; j++)
				kappaistar += selected[j];

			int from = kappai[0]+1;
			int to = std::max(i-1, kappai[0]+1);
			Cjplus[0] = 0;
			for (int j = from; j <= to; j++) {
				if (cand[j] && E[j]-1 <= i-1)
					Cjplus[0]++;
			}

			alphaitilde = w0 * gammai[ S[i-1] - candsum  ] + 
			    (tau*(1-lambda)*alpha-w0)*gammai[ S[i-1] - kappaistar - Cjplus[0] ];
			
		} else {

			alphaitilde = w0*gammai[ S[i-1]-candsum ];

		}

		alphai[i] = std::min(tau*lambda, alphaitilde);
		if (pval[i] <= alphai[i]) {
			R[i] = 1;
	    //K++;
		}
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}
