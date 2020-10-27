#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using std::endl;

// Debugging function to print contents of a vector.
// Need separate functions for each datatype.
void printVec(NumericVector vec) {
	for (int i = 0; i < vec.size(); i++)
		Rcout << vec[i] << " ";
	Rcout << endl;
}
void printVec(IntegerVector vec) {
	for (int i = 0; i < vec.size(); i++)
		Rcout << vec[i] << " ";
	Rcout << endl;
}

// [[Rcpp::export]]
DataFrame saffron_faster(NumericVector pval,
	NumericVector gammai = NumericVector(0),
	double lambda = 0.5,
	double alpha = 0.05,
	double w0 = 0.025) {
	
	int N = pval.size();
	
	if (gammai.size() == 0) {
		gammai = static_cast<NumericVector>(no_init(N));
	// gammai <- 0.4374901658/(seq_len(N)^(1.6))
		for (int i = 0; i < N; i++)
			gammai[i] = 0.4374901658/pow(i+1, 1.6);
	}

	NumericVector alphai(N);
	LogicalVector R(N);
	
	alphai[0] = std::min((1-lambda)*gammai[0]*w0, lambda);
	R[0] = (pval[0] <= alphai[0]);

	int candsum = 0;
	IntegerVector Cjplus(N);
	IntegerVector cand(N);
	IntegerVector tau(1);
	
	int K = sum(R);
	
	for(int i = 1; i < N; i++) {
		
		cand[i-1] = (pval[i-1] <= lambda);
		candsum += cand[i-1];

		double alphaitilde;
		if (K > 1) {
			
			if(R[i-1])
				tau.push_back(i-1);
			
	    //update Cjplus
			double Cjplussum = 0;
			for (int j = 0; j < K-1; j++) {
				Cjplus[j] += cand[i-1];
				Cjplussum += gammai[i - tau[j] - Cjplus[j] - 1];
			}
			
	    //update Cjplus again
			Cjplus[K-1] = 0;
			int low = tau[K-1]+1;
			int high = std::max(i-1, (int)(max(tau) + 1));
			for (int j = low; j <= high; j++) {
				Cjplus[K-1] += cand[j];
			}
			
			Cjplussum += gammai[i-tau[K-1]-Cjplus[K-1]-1]-gammai[i-tau[0]-Cjplus[0]-1];

			alphaitilde = (1 - lambda)*(w0*gammai[i-candsum] +
				(alpha - w0)*gammai[i-tau[0]-Cjplus[0]-1] + alpha*Cjplussum);
			
		} else if (K == 1) {
			
			if(R[i-1])
				tau[0] = i-1;
			
			Cjplus[0] = sum(cand[seq(tau[0] + 1, std::max(i-1, (int)(max(tau + 1))))]);
			alphaitilde = (1 - lambda)*(w0*gammai[i-candsum] +
				(alpha-w0)*gammai[i-tau[0]-Cjplus[0]-1]);
			
		} else {
			alphaitilde = (1-lambda)*w0*gammai[i-candsum];
		}
		
		alphai[i] = std::min(lambda, alphaitilde);
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
List saffron_rcpp_full(int N,
	double lambda,
	double w0,
	double alpha,
	int candsum,
	NumericVector gammai,
	IntegerVector Cjplus,
	IntegerVector cand,
	NumericVector pval,
	NumericVector alphai,
	LogicalVector R) {
	
	IntegerVector tau(1);

	for(int i = 1; i < N; i++) {
		
		int K = sum(R);
		cand[i-1] = (pval[i-1] <= lambda);
		candsum += cand[i-1];
		
		if (K > 1) {
			
			if(R[i-1])
				tau.push_back(i-1);
			
	    // seq_len is 1-based
			IntegerVector Kseq = seq_len(K-1) - 1;
			
	    //update Cjplus
			for (int j = 0; j < Kseq.size(); j++) {
				Cjplus[j] += cand[i-1];
			}

			NumericVector gi = gammai[(tau[Kseq] + Cjplus[Kseq]+1 - i) * -1];
			double Cjplussum = sum(gi);

	    //update Cjplus again
			Cjplus[K-1] = sum(cand[seq(tau[K-1] + 1, std::max(i-1, (int)(max(tau) + 1)))]);
			Cjplussum += gammai[i-tau[K-1]-Cjplus[K-1]-1]-gammai[i-tau[0]-Cjplus[0]-1];

			double alphaitilde = (1 - lambda)*(w0*gammai[i-candsum] +
				(alpha - w0)*gammai[i-tau[0]-Cjplus[0]-1] + alpha*Cjplussum);

			alphai[i] = std::min(lambda, alphaitilde);
			R[i] = (pval[i] <= alphai[i]);
			
		} else if (K == 1) {
			
			if(R[i-1])
				tau[0] = i-1;
			
			Cjplus[0] = sum(cand[seq(tau[0] + 1, std::max(i-1, (int)(max(tau + 1))))]);
			double alphaitilde = (1 - lambda)*(w0*gammai[i-candsum] +
				(alpha-w0)*gammai[i-tau[0]-Cjplus[0]-1]);

			alphai[i] = std::min(lambda, alphaitilde);
			R[i] = pval[i] <= alphai[i];
			
		} else {
			double alphaitilde = (1-lambda)*w0*gammai[i-candsum];
			alphai[i] = std::min(lambda, alphaitilde);
			R[i] = pval[i] <= alphai[i];
		}
	}
	return List::create(_["candsum"] = candsum,
		_["Cjplus"] = Cjplus,
		_["cand"] = cand,
		_["tau"] = tau,
		_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}

