// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <vector>
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
	NumericVector gammai,
	double lambda = 0.25,
	double alpha = 0.05,
	double tau = 0.5,
	double w0 = 0.025,
	bool display_progress = true) {
	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);
	IntegerVector Cjplus(N);
	IntegerVector cand(N);
	LogicalVector selected = (pval <= tau);
	NumericVector S = cumsum(static_cast<NumericVector>(selected));

	alphai[0] = std::min((tau-lambda)*gammai[0]*w0, lambda);
	R[0] = (pval[0] <= alphai[0]);
  
	int K = sum(R);
	int candsum = 0; 
	IntegerVector kappai(1);

	Progress p(N * N, display_progress);

	for (int i = 1; i < N; i++) {

		cand[i-1] = (pval[i-1] <= lambda);
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
				p.increment();
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

			alphaitilde = (tau - lambda)*(w0*gammai[ S[i-1]-candsum ] + 
			(alpha-w0)*gammai[ S[i-1]-kappaistar[0]-Cjplus[0] ] + alpha*Cjplussum);

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

			alphaitilde = (tau - lambda)*(w0*gammai[ S[i-1] - candsum  ] + 
			    (alpha-w0)*gammai[ S[i-1] - kappaistar - Cjplus[0] ]);

		} else {

			alphaitilde = (tau - lambda)*w0*gammai[ S[i-1] - candsum ];

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
DataFrame addis_async_faster(NumericVector pval,
	IntegerVector E,
	NumericVector gammai,
	double lambda = 0.25,
	double alpha = 0.05,
	double tau = 0.5,
	double w0 = 0.025,
	bool display_progress = false) {

	// ADDIS*_async of Tian and Ramdas (2019), Algorithm 3. Test i (0-based) starts at
	// time t = i + 1 and its outcome is known from time t onwards if E[j] < t, i.e.
	// E[j] <= i. Rejection times kappa_j, kappa_j^* and C_j^+ are all defined by
	// decision times.

	int N = pval.size();

	NumericVector alphai(N);
	LogicalVector R(N);

	// Numbers of started tests that are selected (p <= tau) or candidates
	// (p <= lambda), by decision time; decision times beyond N are never known
	std::vector<int> selDec(N + 2, 0), candDec(N + 2, 0);
	std::vector<int> prefSel(N + 1, 0), prefCand(N + 1, 0);
	std::vector<int> kappa;

	Progress p(N, display_progress);

	for (int i = 0; i < N; i++) {

		if (i > 0) {
			int d = std::min(std::max((int)E[i-1], 1), N + 1);
			selDec[d] += (pval[i-1] <= tau);
			candDec[d] += (pval[i-1] <= lambda);
		}

		// prefix sums over decision times 1, ..., i (outcomes known at time i + 1)
		for (int d = 1; d <= i; d++) {
			prefSel[d] = prefSel[d-1] + selDec[d];
			prefCand[d] = prefCand[d-1] + candDec[d];
		}

		int pending = 0;
		kappa.clear();
		for (int j = 0; j < i; j++) {
			if (E[j] > i)
				pending++;
			else if (R[j])
				kappa.push_back(std::max((int)E[j], 1));
		}
		std::sort(kappa.begin(), kappa.end());

		int S = prefSel[i] + pending;   // S^t
		int C0 = prefCand[i];           // C_0^+

		double wealth = w0 * gammai[S - C0];
		for (std::size_t j = 0; j < kappa.size(); j++) {
			int kstar = prefSel[kappa[j]];               // kappa_j^*
			int Cj = C0 - prefCand[kappa[j]];            // C_j^+
			wealth += (j == 0 ? alpha - w0 : alpha) * gammai[S - kstar - Cj];
		}

		alphai[i] = std::min(lambda, (tau - lambda) * wealth);
		R[i] = (pval[i] <= alphai[i]);
		p.increment();
	}

	return DataFrame::create(_["pval"] = pval,
		_["alphai"] = alphai,
		_["R"] = R);
}
