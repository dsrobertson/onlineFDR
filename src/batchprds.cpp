#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using std::endl;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::NumericVector subset_rcpp(Rcpp::NumericVector v, Rcpp::NumericVector pval, double batch) { 
	Rcpp::NumericVector ind(v.size());
	for (int i = 0; i < v.size(); i++){
		ind[i] = (v[i] == batch);
	}

	return pval[ind > 0];
}

// [[Rcpp::export]]
IntegerVector which_rcpp(NumericVector v, double batch) {

	int n = v.size();
	std::vector< int > res;
	res.push_back(0);

	for(int i = 0; i < n; i++) {
		double x = v[i];
		if(x == batch) {
			res.push_back(i);
		}
	}
	res.erase(res.begin());

	Rcpp::IntegerVector iv( res.begin(), res.end() );
	return iv;
}

// [[Rcpp::export]]
Rcpp::NumericVector Rcpp_sort(Rcpp::NumericVector x, Rcpp::NumericVector y) {
    // Order the elements of x by sorting y
    // First create a vector of indices
    Rcpp::IntegerVector idx = seq_along(x) - 1;
    // Then sort that vector by the values of y
    std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
    // And return x in that order
    return x[idx];
}

// [[Rcpp::export]]
arma::rowvec arma_sub(arma::rowvec x, arma::uvec pos, arma::vec vals) {
    
    x.elem(pos) = vals;  // Subset by element position and set equal to
    return x;
}

// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
   NumericVector y = clone(x);
   std::sort(y.begin(), y.end());
   return y;
}

// [[Rcpp::export]]
Rcpp::List prds_faster(DataFrame d,
	Rcpp::NumericVector gammai,
	double n_batch,
	double alpha = 0.05) {

	Rcpp::LogicalVector R(d.nrows());
	Rcpp::NumericVector alphai(n_batch);
	alphai[0] = gammai[0] * alpha;

	for(int i = 1; i <= n_batch; i++) {
		Rcpp::NumericVector batch_pval = subset_rcpp(d["batch2"], d["pval"], i);
		int n = batch_pval.size();

		Rcpp::NumericVector ordered_pval = stl_sort(batch_pval);
		Rcpp::LogicalVector batchR(n);
		for(int j = 0; j < n; j++){
			batchR[j] = ordered_pval[j] <= (j/n)*alphai[i];
		}

		int max_entry = which_max(batchR);
		// arma_sub(batchR, Rcpp::Range(0, max_entry), Rcpp::rep(1, max_entry));

		Rcpp::NumericVector outR = Rcpp_sort(static_cast<NumericVector>(batchR), batch_pval);
		IntegerVector idx = which_rcpp(d["batch2"], i);
		// arma_sub(R, idx, outR);

		if(i < n_batch) {
			int ntplus = which_rcpp(d["batch2"], i+1).size();
			alphai[i+1] = alpha * gammai[i+1]/ntplus * ( ntplus + Rcpp::sum(R));
		}
	}
	return List::create(_["R"] = R, _["alphai"] = alphai);
}

