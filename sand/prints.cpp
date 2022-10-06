
Rcpp::Rcout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
for(int i=0; i < 20; i++)
  Rcpp::Rcout << state_lengths[i] << " ";
Rcpp::Rcout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";


Rcpp::Rcout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
for(int i=0; i < ploidy.nrow(); i++)
  Rcpp::Rcout << "\t" << ploidy(i,0) << " x " << ploidy(i,1) << ": " << n_states[i] <<"\n";


Rcpp::Rcout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
for(int i=0; i < n_states[0]; i++){
  for(int j=0; j < n_states[0]; j++){
    Rcpp::Rcout <<  R[i*n_states[0] + j] << " ";
  }
  Rcpp::Rcout << "\n";
}



Rcpp::Rcout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
for(int k=0; k < n_states.size(); k++){
  for(int i=0; i < n_states[k]; i++){
    for(int j=0; j < n_states[k]; j++){
      Rcpp::Rcout.precision(2);
      Rcpp::Rcout <<  std::fixed << R[i*n_states[k] + j + accum_states[k]] << " ";
    }
    Rcpp::Rcout << "\n";
  }
  Rcpp::Rcout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
}




Rcpp::Rcout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
for(int l=0; l < n_states.size(); l++){
  for(int i=0; i < rf_cur.size(); i++){
    for(int j=0; j < n_states[l]; j++){
      for(int k=0; k < n_states[l]; k++){
        Rcpp::Rcout.precision(2);
        Rcpp::Rcout <<  std::fixed << exp(T[accum_states[l]*(n_mrk-1) + i*n_states[l]*n_states[l] + j*n_states[l] + k]) << " ";
      }
      Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  }
  Rcpp::Rcout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  Rcpp::Rcout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
}


if(i == 1 &&  k == 0 && it == 0)
  Rcpp::Rcout << "id1 + j: " << id1 + j << "\n";
if(i == 1 &&  k == 0 && it == 0)
  Rcpp::Rcout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~\n";
if(i == 1 &&  k == (n_mrk-1) && it == 0)
  Rcpp::Rcout << "id1 + j: " << id1 + j << "\n";






