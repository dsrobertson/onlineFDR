// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
# include <algorithm>

using namespace Rcpp;
using std::endl;

// [[Rcpp::export]]
DataFrame alphainvesting_faster(NumericVector pval,
	NumericVector gammai = NumericVector(0),
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
	
	alphai[0] = gammai[0]*w0/(1+gammai[0]*w0);
	R[0] = (pval[0] <= alphai[0]);

	int candsum = 0;
	IntegerVector Cjplus(N);
	IntegerVector cand(N);
	IntegerVector tau(1);
	
	int K = sum(R);

	Progress p(N * N,true);
	
	for(int i = 1; i < N; i++) {
		
		cand[i-1] = (pval[i-1] <= alphai[i-1]);
		candsum += cand[i-1];

		double alphaitilde;
		if (K > 1) {
			
			if(R[i-1])
				tau.push_back(i-1);
			
	    //update Cjplus
			double Cjplussum = 0;
			for (int j = 0; j < K-1; j++) {
				p.increment();
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

			alphaitilde = (w0*gammai[i-candsum] +
				(alpha - w0)*gammai[i-tau[0]-Cjplus[0]-1] + alpha*Cjplussum);
			
		} else if (K == 1) {
			
			if(R[i-1])
				tau[0] = i-1;
			
			Cjplus[0] = sum(cand[seq(tau[0] + 1, std::max(i-1, (int)(max(tau + 1))))]);
			alphaitilde = (w0*gammai[i-candsum] +
				(alpha-w0)*gammai[i-tau[0]-Cjplus[0]-1]);
			
		} else {
			alphaitilde = w0*gammai[i-candsum];
		}
		
		alphai[i] = alphaitilde/(1+alphaitilde);
		if (pval[i] <= alphai[i]) {
			R[i] = 1;
			K++;
		}
		
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}

