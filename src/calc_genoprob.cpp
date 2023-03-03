/*
 MAPpoly: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2022 Marcelo Mollinari

 This file is part of MAPpoly.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

/*
 File: est_hmm_map.cpp

 Description:

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: Oct 06, 2022
 */

#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "combinatorial.h"
#include "hmm_elements.h"
using namespace std;
using namespace Rcpp;

RcppExport SEXP calc_genoprob_multi_fam(SEXP ploidy1_idR,
                                        SEXP ploidy2_idR,
                                        SEXP n_mrkR,
                                        SEXP n_indR,
                                        SEXP haploR,
                                        SEXP emitR,
                                        SEXP rfR,
                                        SEXP verboseR)
{
  //convert input to C++ types
  Rcpp::NumericVector p1(ploidy1_idR);
  Rcpp::NumericVector p2(ploidy2_idR);
  int n_mrk = Rcpp::as<int>(n_mrkR);
  int n_ind = Rcpp::as<int>(n_indR);
  Rcpp::List haplo(haploR);
  Rcpp::List emit(emitR);
  Rcpp::NumericVector rf(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  std::vector<int> pl{2,4,6};
  Rcpp::NumericMatrix Mat(n_mrk, 36);

  //Initializing some variables
  int k, k1, count = 0, v_len = 0;
  long double s;
  int mpi1 = max(p1);
  int mpi2 = max(p2);
  int max_ploidy_id;
  if(mpi1 >= mpi2)
    max_ploidy_id = mpi1;
  else
    max_ploidy_id = mpi2;
  for(int i; i < p1.size(); i++)
    v_len += (pl[p1[i]] + pl[p2[i]]);
  v_len *= n_mrk;
  Rcpp::NumericMatrix pr_h(v_len, 5);
  std::vector<std::vector<std::vector<int>>> g{{{0},{1}},
                                               {{0,1},{0,2},{0,3},
                                                {1,2},{1,3},{2,3}},
                                               {{0,1,2},{0,1,3},{0,1,4},{0,1,5},
                                                {0,2,3},{0,2,4},{0,2,5},{0,3,4},
                                                {0,3,5},{0,4,5},{1,2,3},{1,2,4},
                                                {1,2,5},{1,3,4},{1,3,5},{1,4,5},
                                                {2,3,4},{2,3,5},{2,4,5},{3,4,5}}};

  //Initializing v: states hmm should visit for each marker
  //Initializing e: emission probabilities associated to the states hmm should visit for each marker
  std::vector<std::vector<std::vector<int> > > v;
  std::vector<std::vector<std::vector<double> > > e;
  for(int i=0; i < haplo.size(); i++) // i: number of markers
  {
    Rcpp::List haplo_temp(haplo(i)); //states hmm should visit for marker i
    Rcpp::List emit_temp(emit(i)); //emission probs. for states hmm should visit for marker i
    std::vector<std::vector<int> > v1;
    std::vector<std::vector<double> > e1;
    for(int j=0; j < haplo_temp.size(); j++) //iterate for all j individuals
    {
      Rcpp::NumericMatrix M_temp = haplo_temp(j);
      Rcpp::NumericVector E_temp = emit_temp(j);
      std::vector<int> v2 = Rcpp::as<std::vector<int> >(M_temp);
      std::vector<double> e2 = Rcpp::as<std::vector<double> >(E_temp);
      v1.push_back(v2);
      e1.push_back(e2);
    }
    v.push_back(v1);
    e.push_back(e1);
  }

  //Initializing alpha and beta
  std::vector<std::vector<std::vector<long double> > > alpha(n_ind);
  std::vector<std::vector<std::vector<long double> > > beta(n_ind);
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int i=0; i < n_mrk; i++)
    {
      std::vector<long double> temp3(v[i][ind].size()/2);
      alpha[ind].push_back(temp3);
      beta[ind].push_back(temp3);
    }
  }

  //Initializing transition matrices
  std::vector< std::vector< std::vector< std::vector<double> > > >T;
  for(int j=0; j <= max_ploidy_id; j++)
  {
    std::vector< std::vector< std::vector<double> > > Ttemp;
    for(int i=0; i < n_mrk-1; i++)
    {
      Ttemp.push_back(transition(pl[j], rf[i]));
    }
    T.push_back(Ttemp);
  }

  //Loop over all individuals
  for(int ind=0; ind < n_ind; ind++)
  {
    R_CheckUserInterrupt();
    for(int j=0; (unsigned)j < e[0][ind].size(); j++)
    {
      alpha[ind][0][j] = e[0][ind][j];
    }
    std::fill(beta[ind][n_mrk-1].begin(), beta[ind][n_mrk-1].end(), 1);
    //forward-backward
    for(k=1,k1=n_mrk-2; k < n_mrk; k++, k1--)
    {
      std::vector<long double> temp4 (v[k][ind].size()/2);
      temp4 = forward_emit_highprec(
                           alpha[ind][k-1],
                           v[k-1][ind],
                           v[k][ind],
                           e[k][ind],
                           T[p1[ind]][k-1],
                           T[p2[ind]][k-1]);
      for(int j=0; (unsigned)j < temp4.size(); j++)
      {
        alpha[ind][k][j]=temp4[j];
      }
      std::vector<long double> temp5 (v[k1][ind].size()/2);
      temp5 = backward_emit_highprec(
                            beta[ind][k1+1],
                            v[k1][ind],
                            v[k1+1][ind],
                            e[k1+1][ind],
                            T[p1[ind]][k1],
                            T[p2[ind]][k1]);
      for(int j=0; (unsigned)j < temp5.size(); j++)
      {
        beta[ind][k1][j]=temp5[j];
      }
    }
  } // loop over individuals

  for(int ind=0; ind < n_ind; ind++)
  {
    for(int k=0; k < n_mrk; k++)
    {
      std::vector<long double> w(alpha[ind][k].size());
      s = w[0] = alpha[ind][k][0] * beta[ind][k][0];
      for(int j=1; (unsigned)j < w.size(); j++){
        w[j] = alpha[ind][k][j] * beta[ind][k][j];
        s += w[j];
      }
      for(int j=0; (unsigned)j < w.size(); j++)
        w[j] = w[j]/s;

      std::vector<double> u1(pl[p1[ind]], 0);
      std::vector<double> u2(pl[p2[ind]], 0);

      for(int j=0; (unsigned)j < w.size(); j++)
      {
        for(int l=0; l < pl[p1[ind]]/2; l++){
          u1[g[p1[ind]][v[k][ind][j]][l]] += w[j];
        }
        for(int l=0; l < pl[p2[ind]]/2; l++){
          u2[g[p2[ind]][v[k][ind][j + v[k][ind].size()/2]][l]] += w[j];
        }
      }
      for(int l=0; (unsigned)l < u1.size(); l++){
        pr_h(count,0) =  ind + 1;
        pr_h(count,1) =  k + 1;
        pr_h(count,2) =  1;
        pr_h(count,3) =  l+1;
        pr_h(count,4) =  u1[l];
        count++;
      }
      for(int l=0; (unsigned)l < u2.size(); l++){
        pr_h(count,0) =  ind + 1;
        pr_h(count,1) =  k + 1;
        pr_h(count,2) =  2;
        pr_h(count,3) =  l+1;
        pr_h(count,4) =  u2[l];
        count++;
      }
    }
  }
  return pr_h;
}
