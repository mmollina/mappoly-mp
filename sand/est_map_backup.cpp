#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_utils.h"
#include "est_map.h"
#include <thread>
#include <chrono>
#define THRESH 200.0

using namespace std;
using namespace Rcpp;

/* FUNCTION: est_hmm_map
 * -----------------------------------------------------
 */
RcppExport SEXP est_log_hmm_map(SEXP ploidyR,
                                SEXP n_mrkR,
                                SEXP n_indR,
                                SEXP rfR,
                                SEXP crosstypeR,
                                SEXP states_lenR,
                                SEXP cum_states_lenR,
                                SEXP statesR,
                                SEXP emitR,
                                SEXP tolR,
                                SEXP verboseR){
  //convert input to C++ types
  Rcpp::NumericMatrix ploidy(ploidyR); //(n.cross x 2)
  std::vector<int> n_states(ploidy.nrow());
  std::vector<int> accum_states(ploidy.nrow() + 1);
  accum_states[0] = 0;
  for(int i=0; i < ploidy.nrow(); i++)
  {
    n_states[i] = nChoosek(ploidy(i,0), ploidy(i,0)/2) *
      nChoosek(ploidy(i,1), ploidy(i,1)/2) ;
    accum_states[i+1] = accum_states[i] + n_states[i] * n_states[i];
  }
  int n_mrk = Rcpp::as<int>(n_mrkR);
  int n_ind = Rcpp::as<int>(n_indR);
  std::vector<double> rf = Rcpp::as<std::vector<double> >(rfR);
  std::vector<int> crosstype = Rcpp::as<std::vector<int> >(crosstypeR);
  std::vector<int> states_len = Rcpp::as<std::vector<int> >(states_lenR);
  std::vector<int> cum_states_len = Rcpp::as<std::vector<int> >(cum_states_lenR);
  std::vector<int> states = Rcpp::as<std::vector<int> >(statesR);
  std::vector<double> emit = Rcpp::as<std::vector<double> >(emitR);
  double tol = Rcpp::as<double>(tolR);
  int verbose = Rcpp::as<int>(verboseR);
  std::vector<double> alpha(emit.size());
  std::vector<double> beta(emit.size());
  int k, k1, id1, id2, id3, id4, maxit = 1000, flag;
  std::vector<double> term(n_ind), rf_cur(rf.size());
  double loglike = 0;
  // NUMBER OF RECOMBINATIONS
  // i: state in
  // j: state out
  // k: cross type
  // R[i*n_states[k] + j + accum_states[k]]
  std::vector<double> R;
  std::vector<double> Rtemp;
  for(int i=0; i < ploidy.nrow(); i++)
  {
    Rtemp = rec_num(ploidy(i,0), ploidy(i,1));
    R.insert( R.end(), Rtemp.begin(), Rtemp.end() );
  }

  //begin EM algorithm
  for(int it=0; it < maxit; it++)
  {
    //Rcpp::Rcout << ".";
    //Initializing recombination fraction vector for Baum-Welch
    for(int j=0; j<n_mrk-1; j++)
    {
      rf_cur[j] = rf[j];
      rf[j] = 0.0;
    }
    // TRANSITION
    // i: recombination fraction (rf[0]: between first and second)
    // j: state in
    // k: state out
    // l: cross type
    // T[accum_states[l]*(n_mrk-1) + i*n_states[l]*n_states[l] + j*n_states[l] + k]
    std::vector<double> T;
    std::vector<double> Ttemp;
    for(int i=0; i < ploidy.nrow(); i++)
    {
      Ttemp = log_transition(ploidy(i,0), ploidy(i,1), rf_cur);
      T.insert( T.end(), Ttemp.begin(), Ttemp.end() );
    }

    // STATES and EMISSION
    // i: individual
    // j: marker
    // k: states index
    //for(k = 0; (unsigned)k < states_len[j + i*n_mrk]; k++)
    //  alpha[cum_states_len[j + i*n_mrk] + k] = emit[cum_states_len[j + i*n_mrk] + k];

    //Start loop over all individuals
    for(int i=0; i < n_ind; i++)
    {
      R_CheckUserInterrupt();
      for(int j = cum_states_len[0 + i*n_mrk];
          (unsigned)j < cum_states_len[1 + i*n_mrk];
          j++) // fist marker
        alpha[j] = emit[j];
      for(int j = cum_states_len[(n_mrk-1) + i*n_mrk];
          (unsigned)j < cum_states_len[n_mrk + i*n_mrk];
          j++) // last marker
        beta[j] = 0.0;

      //forward
      for(k = 1, k1 = n_mrk-2; k < n_mrk; k++, k1--)
      {
        for(int j1 = 0; (unsigned)j1 < states_len[k + i*n_mrk]; j1++)
        {
          id1 = cum_states_len[k + i*n_mrk];
          id2 = cum_states_len[k - 1 + i*n_mrk];
          alpha[id1 + j1] = alpha[id2] +
            T[accum_states[crosstype[i]]*(n_mrk-1) +
              (k-1)*n_states[crosstype[i]]*n_states[crosstype[i]] +
              states[id2] * n_states[crosstype[i]] + states[id1+j1]];
          for(int j = 1; (unsigned)j < states_len[k - 1 + i*n_mrk];j++) {
            {
              alpha[id1+j1] = addlog(alpha[id1 + j1],
                                     alpha[id2 + j] +
                                    T[accum_states[crosstype[i]]*(n_mrk-1) +
                                    (k-1)*n_states[crosstype[i]]*n_states[crosstype[i]] +
                                    states[id2 + j] * n_states[crosstype[i]] + states[id1+j1]]
              );
            }
            alpha[id1] += emit[id1];
          }
        }
      }


      //backward
      for(k = 1, k1 = n_mrk-2; k < n_mrk; k++, k1--)
      {
        for(int j1 = 0; (unsigned)j1 < states_len[k1 + i*n_mrk]; j1++)
        {
          id3 = cum_states_len[k1 + i*n_mrk];
          id4 = cum_states_len[k1 + 1 + i*n_mrk];
          beta[id3 + j1] = beta[id4] +
            T[accum_states[crosstype[i]]*(n_mrk-1) +
            k1*n_states[crosstype[i]]*n_states[crosstype[i]] +
            states[id3 + j1] * n_states[crosstype[i]] + states[id4]] +
            emit[id4];

          for(int j = 1; (unsigned)j < states_len[k1 - 1 + i*n_mrk];j++) {
            {
              beta[id3+j1] = addlog(beta[id3 + j1],
                                    beta[id4 + j] +
                                      T[accum_states[crosstype[i]]*(n_mrk-1) +
                                      k1*n_states[crosstype[i]]*n_states[crosstype[i]] +
                                      states[id3 + j1] * n_states[crosstype[i]] + states[id4 + j]] +
                                      emit[id4 + j]
              );
            }
          }
        }
      }



      //Updating recombination fraction
      for(int k = 0; k < n_mrk-1; k++)
      {
        int n_gen_k = states_len[k + i*n_mrk];
        int n_gen_k1 = states_len[(k+1) + i*n_mrk];
        vector<vector<double> > gamma(n_gen_k, vector<double>(n_gen_k1));

        double s = 0.0;
        for(int v = 0; (unsigned)v < n_gen_k; v++)
        {
          id1 = cum_states_len[k + i*n_mrk] + v;
          for(int v2 = 0; (unsigned)v2 < n_gen_k1; v2++)
          {
            id2 = cum_states_len[k + 1 + i*n_mrk] + v2;
            gamma[v][v2] = alpha[cum_states_len[k + i*n_mrk] + v] +
                           beta[cum_states_len[(k+1) + i*n_mrk] + v2] +
                           T[accum_states[crosstype[i]]*(n_mrk-1) +
                             k*n_states[crosstype[i]]*n_states[crosstype[i]] +
                             states[id1] * n_states[crosstype[i]] + states[id2]];
            if(v==0 && v2==0) s = gamma[v][v2];
            else s = addlog(s, gamma[v][v2]);
          }
        }
        for(int v = 0; (unsigned)v < n_gen_k; v++)
        {
          id1 = cum_states_len[k + i*n_mrk] + v;
          for(int v2 = 0; (unsigned)v2 < n_gen_k1; v2++)
          {
            id2 = cum_states_len[k + 1 + i*n_mrk] + v2;
            if(s != 0){
              rf[k] +=  R[states[id1] * n_states[crosstype[i]] +
                          states[id2] +
                          accum_states[crosstype[i]]] *
                        exp(gamma[v][v2] - s);
              // if(i == 0)
              //   Rcpp::Rcout << "k: " << k <<
              //                  " v: " << v <<
              //                  " v2: " << v2 <<
              //                    " id1: " << id1 <<
              //                      " id2: " << id2 <<
              //                        "   st_in: " << states[id2] <<
              //                          " ---> st_out: " <<
              //                            states[id1] << " \n";
            }
          }
        }
      }
    }//End loop over all individuals
    // re-scale
    for(int j=0; j<n_mrk-1; j++)
    {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
    }
    //check convergence
      flag=0;
    for(int j=0; j < n_mrk-1; j++)
    {
      if(fabs(rf[j] - rf_cur[j]) > tol*(rf_cur[j]+tol*100.0))
      {
        flag = 1;
        break;
      }
    }
    //Rcpp::Rcout << "\t\n Iter: " << it+1 << "\t";
    for(int j = 0; j < n_mrk-1; j++)
    {
      //Rcpp::Rcout.precision(3);
      //Rcpp::Rcout << std::fixed << rf[j] << " ";
    }
    if(!flag) break;
  }//end of EM algorithm

  //Termination
  for(int i=0; i < n_ind; i++){
    term[i] = alpha[cum_states_len[(n_mrk-1) + i*n_mrk]];
    for(int k = 0; (unsigned)k < states_len[(n_mrk-1) + i*n_mrk]; k++)
      term[i] = addlog(term[i], alpha[cum_states_len[(n_mrk-1) + i*n_mrk] + k]);
    loglike += term[i];
  }
  //if(verbose)
   // Rcpp::Rcout << "\n";

  List z = List::create(wrap(loglike), rf);
  return(z);
}
