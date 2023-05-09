#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix gamete_index_combination(int x, int y) {
  int ngam1 = R::choose(x, x / 2);
  int ngam2 = R::choose(y, y / 2);

  NumericMatrix S(ngam2 * ngam1, 2);
  int k = 0;

  for(int j = 0; j < ngam1; j++) {
    for(int m = 0; m < ngam2; m++) {
      S(k, 0) = j;
      S(k, 1) = m;
      k++;
    }
  }

  return S;
}

// Modify the original function
// [[Rcpp::export]]
List myfunc(DataFrame ped_mat, DataFrame ind_pop, IntegerMatrix geno) {
  int n = ped_mat.nrows();
  List L(n);

  //
  CharacterVector indnames = colnames(geno);
  CharacterVector mrknames = colnames(geno);

  for(int j=0; j < indnames.size(); j++){
    Rcpp::Rcout << indnames[j] << " ";
  }
  Rcpp::Rcout << "\n";

  //Converting Data Frame to vector
  //ped_mat
  CharacterVector P1 = ped_mat[0];
  CharacterVector P2 = ped_mat[1];
  IntegerVector pl1 = ped_mat[2];
  IntegerVector pl2 = ped_mat[3];
  //ind_pop
  CharacterVector ind_names = ind_pop[0];
  IntegerVector pop_id = ind_pop[1];
  //For each parental cross combination
  for(int i = 0; i < n; i++) {
    //Get indices for transition matrices
    L[i] = gamete_index_combination(pl1[i], pl2[i]);
    CharacterVector cur_ind_names(0);
    //Get names of individuals for each cross
    for (int j = 0; j < pop_id.size(); ++j) {
      if (pop_id[j] == i) {
        cur_ind_names.push_back(ind_names[j]);
      }
    }
    for(int j=0; j < cur_ind_names.size(); j++){
      Rcpp::Rcout << cur_ind_names[j] << " ";
    }
    Rcpp::Rcout << "\n";
  }
  return L;
}
