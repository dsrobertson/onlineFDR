#include <Rcpp.h>
using namespace Rcpp;
using std::endl;

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

	for(int i = 0; i < n; i++) {
		double x = v[i];
		if(x == batch) {
			res.push_back(i);
		}
	}

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
	Rcpp::List R_list = List::create(R);
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
		//turn batchR into list for subsetting to work
		Rcpp::List batchR_list = List::create(batchR);
		batchR_list_new = batchR[seq_len(max_entry)-1];
		//convert back to vector
		NumericVector batchR_vec = Rcpp::as<std::vector<double>>(batchR_list_new);

		Rcpp::NumericVector outR = Rcpp_sort(batchR_vec, batch_pval);
		IntegerVector idx = which_rcpp(d["batch2"], i);
		R_list[idx]

		if(i < n_batch) {
			int ntplus = which_rcpp(d["batch2"], i+1).size();
			alphai[i+1] = alpha * gammai[i+1]/ntplus * ( ntplus + Rcpp::sum(R));
		}
	}
	return List::create(_["R"] = R, _["alphai"] = alphai);
}

