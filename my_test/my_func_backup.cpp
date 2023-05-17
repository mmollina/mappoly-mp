#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <string_view>
using namespace Rcpp;

// Helper function to calculate the combinations
IntegerMatrix combn(NumericVector x, int m) {
  int n = x.size();
  if (m > n) return IntegerMatrix(0);
  if (m == n) return IntegerMatrix(x);

  int count = Rf_choose(n, m);
  IntegerMatrix result(m, count);

  std::vector<int> indices(m);
  std::iota(indices.begin(), indices.end(), 0);

  for (int col = 0; col < count; ++col) {
    for (int row = 0; row < m; ++row) {
      result(row, col) = x[indices[row]];
    }

    for (int i = m - 1; i >= 0; --i) {
      if (++indices[i] <= n - m + i) {
        while (++i < m) {
          indices[i] = indices[i - 1] + 1;
        }
        break;
      }
    }
  }

  return result;
}

// Function for rep_each
IntegerVector rep_each(IntegerVector x, int n) {
  int sz = x.size() * n;
  IntegerVector res(sz);
  for (int i = 0; i < x.size(); ++i) {
    std::fill_n(res.begin() + i * n, n, x[i]);
  }
  return res;
}

// Function for rep_len
IntegerVector rep_len(IntegerVector x, int n) {
  IntegerVector res(n);
  std::copy_n(x.begin(), std::min<int>(n, x.size()), res.begin());
  for (int i = x.size(); i < n; i += x.size()) {
    std::copy_n(x.begin(), std::min<int>(n - i, x.size()), res.begin() + i);
  }
  return res;
}

// Helper function for expand_grid
IntegerMatrix expand_grid(IntegerVector v1, IntegerVector v2) {
  int n_v1 = v1.size();
  int n_v2 = v2.size();
  IntegerVector res_v1 = rep_len(v1, n_v1 * n_v2); // Adjust the rep_len() usage
  IntegerVector res_v2 = rep_each(v2, n_v1); // Adjust the rep_each() usage
  return cbind(res_v2, res_v1);
}

// Helper function for which
IntegerVector which(LogicalVector x) {
  IntegerVector idx = seq_len(x.size());
  IntegerVector out;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i]) {
      out.push_back(idx[i]);
    }
  }
  return out;
}

// Helper function for concatenating vectors
IntegerVector concatenate_vectors(IntegerVector x1, IntegerVector x2) {
  int n1 = x1.size();
  int n2 = x2.size();
  IntegerVector result(n1 + n2);
  std::copy(x1.begin(), x1.end(), result.begin());
  std::copy(x2.begin(), x2.end(), result.begin() + n1);
  return result;
}

// Helper function for calculate L and initialize H
List calculate_L_and_initialize_H(int n_fullsib_pop,
                                  int n_mrk,
                                  int n_ind,
                                  NumericVector pl1,
                                  NumericVector pl2) {
  List L(n_fullsib_pop);
  for (int i = 0; i < n_fullsib_pop; i++) {
    int ngam1 = R::choose(pl1[i], pl1[i] / 2);
    int ngam2 = R::choose(pl2[i], pl2[i] / 2);
    IntegerMatrix S = expand_grid(seq_len(ngam2) - 1, seq_len(ngam1) - 1);
    L[i] = S;
  }
  List H(n_mrk); //H[[n_mrk]][[n_ind]]
  for (int k = 0; k < n_mrk; k++) {
    H[k] = List(n_ind);
  }
  return List::create(Named("L") = L, Named("H") = H);
}




// [[Rcpp::export]]
List vs_multiallelic_Rcpp(List PH,
                          List GENO,
                          NumericMatrix unique_pop_mat,
                          NumericMatrix pedigree) {

  int n_fullsib_pop = unique_pop_mat.nrow();
  int n_ind = pedigree.nrow();
  NumericMatrix temp_mat = PH[0];
  int n_mrk = temp_mat.nrow();
  NumericVector pl1 = unique_pop_mat(_,2);
  NumericVector pl2 = unique_pop_mat(_,3);
  List result = calculate_L_and_initialize_H(n_fullsib_pop, n_mrk, n_ind, pl1, pl2);
  List L = result["L"];
  List H = result["H"];
  for(int k = 0; k < n_mrk; k++) { // ***************************************** Markers
    List H_k(n_ind);
    for(int pop_id = 0; pop_id < n_fullsib_pop; pop_id ++){ // ************ Pop id
      IntegerVector ind_id = which(pedigree(_,4) == pop_id + 1) - 1;
      NumericMatrix matrix_PH1 =  PH[unique_pop_mat(pop_id,0) - 1];
      NumericMatrix matrix_PH2 = PH[unique_pop_mat(pop_id,1) - 1];
      NumericVector v1 = matrix_PH1(k, _);
      NumericVector v2 = matrix_PH2(k, _);
      if (any(is_na(v1)).is_true() || any(is_na(v2)).is_true()) {
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind Parents NA
          H_k[ind_id[j]] = L[pop_id];
        } // end individual loop when NA
      } else {
        IntegerMatrix x1 = combn(v1, v1.size() / 2);
        IntegerMatrix x2 = combn(v2, v2.size() / 2);
        IntegerMatrix x(x1.nrow() + x2.nrow(), x1.ncol() * x2.ncol());
        for (int i = 0; i < x1.ncol(); i++) {
          for (int j = 0; j < x2.ncol(); j++) {
            IntegerVector temp = concatenate_vectors(x1(_, i), x2(_, j));
            temp.sort();
            x(_, (i * x2.ncol()) + j) = temp;
          }
        }
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind ALL
          IntegerMatrix matrix_GENO = as<IntegerMatrix>(GENO[ind_id[j]]);
          IntegerVector a = matrix_GENO(k, _);
          if (is_true(any(is_na(a)))) {  // ************************************* Ind NA
            H_k[ind_id[j]] = L[pop_id];
          } else{  // *********************************************************** Ind Visit
            a.sort();
            IntegerVector y;
            for (int i = 0; i < x.ncol(); i++) {
              bool all_equal = true;
              for (int l = 0; l < x.nrow(); l++) {
                if (x(l, i) != a[l]) {
                  all_equal = false;
                  break;
                }
              }
              if (all_equal) {
                y.push_back(i);
              }
            }
            IntegerMatrix L_pop = L[pop_id];
            IntegerMatrix subset_L_pop(y.size(), L_pop.ncol());
            for (int row = 0; row < y.size(); ++row) {
              subset_L_pop(row, _) = L_pop(y[row], _);
            }
            H_k[ind_id[j]] = subset_L_pop;
          }
        }
      }
    } // end population loop
    H[k] = H_k;
  } // end marker loop
  return H;
}

// [[Rcpp::export]]
List vs_biallelic_Rcpp(List PH,
                       IntegerMatrix G,
                       NumericMatrix unique_pop_mat,
                       NumericMatrix pedigree) {
  int n_fullsib_pop = unique_pop_mat.nrow();
  int n_ind = pedigree.nrow();
  NumericMatrix temp_mat = PH[0];
  int n_mrk = temp_mat.nrow();
  NumericVector pl1 = unique_pop_mat(_,2);
  NumericVector pl2 = unique_pop_mat(_,3);
  List result = calculate_L_and_initialize_H(n_fullsib_pop, n_mrk, n_ind, pl1, pl2);
  List L = result["L"];
  List H = result["H"];
  for(int k = 0; k < n_mrk; k++) { // ***************************************** Markers
    List H_k(n_ind);
    for(int pop_id = 0; pop_id < n_fullsib_pop; pop_id ++){ // ************ Pop id
      IntegerVector ind_id = which(pedigree(_,4) == pop_id + 1) - 1;
      NumericMatrix matrix_PH1 =  PH[unique_pop_mat(pop_id,0) - 1];
      NumericMatrix matrix_PH2 = PH[unique_pop_mat(pop_id,1) - 1];
      NumericVector v1 = matrix_PH1(k, _);
      NumericVector v2 = matrix_PH2(k, _);
      if (any(is_na(v1)).is_true() || any(is_na(v2)).is_true()) {
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind Parents NA
          H_k[ind_id[j]] = L[pop_id];
        } // end individual loop when NA
      } else {
        IntegerMatrix x1 = combn(v1, v1.size() / 2);
        IntegerMatrix x2 = combn(v2, v2.size() / 2);
        IntegerVector x(x1.ncol() * x2.ncol());
        for (int i = 0; i < x1.ncol(); i++) {
          for (int j = 0; j < x2.ncol(); j++) {
            x((i * x2.ncol()) + j) = sum(x1(_ ,i)) + sum(x2(_ ,j));
          }
        }
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind ALL
          if (R_IsNA(G(k, ind_id[j]))) {  // ************************************* Ind NA
            H_k[ind_id[j]] = L[pop_id];
          } else{  // *********************************************************** Ind Visit
            IntegerVector y;
            for (int i = 0; i < x.size(); i++) {
              bool all_equal = true;
              if (x(i) == G(k, ind_id[j])) {
                y.push_back(i);
              }
            }
            IntegerMatrix L_pop = L[pop_id];
            IntegerMatrix subset_L_pop(y.size(), L_pop.ncol());
            for (int row = 0; row < y.size(); ++row) {
              subset_L_pop(row, _) = L_pop(y[row], _);
            }
            H_k[ind_id[j]] = subset_L_pop;
          }
        }
      }
    } // end population loop
    H[k] = H_k;
  } // end marker loop
  return H;
}
