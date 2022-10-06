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
 File: hmm_elements.cpp

 Description:

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: Oct 06, 2022
 */


#include <algorithm>
#include <iostream>
#include <vector>
#include "combinatorial.h"
#include "hmm_elements.h"
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/PrtUtil.h>

using namespace std;
using namespace Rcpp;

/* FUNCTION: prob_k1_given_k_l_m
 This is equation 5 on the paper
 -----------------------------------------------------
 Calculates the genotypic transition probability based
 on l, whichdenotes the number of recombinant bivalents
 between loci k and k + 1 in one parent.
 */
double prob_k1_given_k_l_m(int m, int l, double rf)
{
  return ((pow((1-rf),(m/2-l))*pow(rf,l))/nChoosek(m/2, l));
}

/* FUNCTION: rec_num
 -----------------------------------------------------
 Returns a matrix containing the number of recombination
 events between loci k and k + 1 in one parent given the
 ploidy level  m.
 */
std::vector<std::vector<double> > rec_num(int m)
{
  int g = nChoosek(m, m/2);
  std::vector<std::vector<double> > R(g);
  for(int i = 0; (unsigned)i < R.size(); ++i)
  {
    for(int j = 0; j < g; ++j)
    {
      R[i].push_back(n_rec_given_genk_and_k1(m, i+1, j+1)/(double)m);
      //R[i].push_back(n_rec_given_genk_and_k1(m, i+1, j+1)/((double)m/2.0));
    }
  }
  return(R);
}

/* FUNCTION: transition
 -----------------------------------------------------
 Returns a transition matrix between loci k and k + 1 in
 one parent i.e. Prop(p_{k+1}|p_k), given the ploidy level
 m and the recombination fraction rf.
 */
std::vector<std::vector<double> > transition(int m, double rf)
{
  int g = nChoosek(m, m/2);
  std::vector<std::vector<double> > T(g);
  for(int i = 0; i < g; ++i)
  {
    for(int j = 0; j < g; ++j)
    {
      T[i].push_back(prob_k1_given_k_l_m(m, n_rec_given_genk_and_k1(m, i+1, j+1), rf));
    }
  }
  return(T);
}

/* FUNCTION:  index_func
 * -----------------------------------------------------
 * This function has a similar purpose as the emission function.
 * It return the indices corresponding to
 * states that should be visited given two vectors indicating
 * which homologous contain the allelic variant.
 * The result is a vector twice the size of states that should be
 * visited. The first half indicates the indices in the transition
 * space in parent P and the second half indicates the indices in
 * the transition space in parent Q
 */

std::vector<std::vector<int> > index_func(int m,
                                          std::vector<int>& p,
                                          std::vector<int>& q)
{
  int s, ip=0, iq=0;
  //int g = nChoosek(m, m/2);
  std::vector<std::vector<int> > v1(p.size()+q.size()+2);
  std::vector<std::vector<int> > v2(p.size()+q.size()+2);
  std::vector<int> vp(m), vq(m);
  std::fill(vp.begin(), vp.end()-m/2, true);
  std::fill(vq.begin(), vq.end()-m/2, true);
  do
  {
    iq=0;
    do
    {
      s=0;
      for(int j=0; (unsigned)j < p.size(); j++)
        if(p[j]>=0) s=s+vp[p[j]];
        for(int j=0; (unsigned)j < q.size(); j++)
          if(q[j]>=0) s=s+vq[q[j]];
          v1[s].push_back(ip);
          v2[s].push_back(iq);
          v1[v1.size()-1].push_back(ip);
          v2[v2.size()-1].push_back(iq);
          iq++;
    }
    while (std::prev_permutation(vq.begin(), vq.end()));
    ip++;
  }
  while (std::prev_permutation(vp.begin(), vp.end()));
  for(int i=0; (unsigned)i < v1.size(); i++)
  {
    v1[i].insert(v1[i].end(), v2[i].begin(), v2[i].end());
  }
  return(v1);
}

/*FUNCTION: init_poly
 ------------------------------------------------------------------
 Description: This function computes the probability of a certain
 genotype dG in a F1 population, given the genotypes (dosage) of its
 parents (dP and dQ) and the ploidy level (m). This is also know as
 polysomic segregation. It returns the log of the probability
 */
double init_poly(int m, int dP, int dQ, int dG)
{
  int j, i;
  double seg = 0.0;
  for(i=0, j=dG; i <= dG; i++, j--)
    seg += R::dhyper(i, dP, m-dP, m/2, 0) * R::dhyper(j, dQ, m-dQ, m/2, 0);
  return(seg);
}

/* FUNCTION: forward
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<double> forward(int m,
                            std::vector<double>& fk,
                            std::vector<int>& ik,
                            std::vector<int>& ik1,
                            std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
  }
  return(fk1);
}
/* FUNCTION: backward
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<double> backward(int m,
                             std::vector<double>& fk1,
                             std::vector<int>& ik,
                             std::vector<int>& ik1,
                             std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
  }
  return(fk);
}



/* FUNCTION: forward_emit (with both informative parents)
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<double> forward_emit(int m,
                                 std::vector<double>& fk,
                                 std::vector<int>& ik,
                                 std::vector<int>& ik1,
                                 std::vector<double>& emit,
                                 std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
    fk1[k1] = fk1[k1] * emit[k1];
  }
  return(fk1);
}
/* FUNCTION: backward (with both informative parents)
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<double> backward_emit(int m,
                                  std::vector<double>& fk1,
                                  std::vector<int>& ik,
                                  std::vector<int>& ik1,
                                  std::vector<double>& emit,
                                  std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]] * emit[k1];
    }
  }
  return(fk);
}

/* FUNCTION: forward (with one informative parent)
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<double> forward_emit_one_parent(int m,
                                            std::vector<double>& fk,
                                            std::vector<int>& ik,
                                            std::vector<int>& ik1,
                                            std::vector<double>& emit,
                                            std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size();
  int ngenk1 = ik1.size();
  std::vector<double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T[ik[k]][ik1[k1]];
    }
    fk1[k1] = fk1[k1] * emit[k1];
  }
  return(fk1);
}
/* FUNCTION: backward (with one informative parent)
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<double> backward_emit_one_parent(int m,
                                             std::vector<double>& fk1,
                                             std::vector<int>& ik,
                                             std::vector<int>& ik1,
                                             std::vector<double>& emit,
                                             std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size();
  int ngenk1 = ik1.size();
  std::vector<double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T[ik[k]][ik1[k1]] * emit[k1];
    }
  }
  return(fk);
}
