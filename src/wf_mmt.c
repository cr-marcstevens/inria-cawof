/* src/wf_mmt.c
 * 
 * Copyright (C) 2015, 2016, 2017 Rodolfo Canto Torres 
 *
 * This file is part of CaWoF.
 *
 * CaWoFa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CaWoFa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CaWoF.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file src/wf_mmt.h
 * \brief Implementation of include/wf_mmt.h
 */

#include<wf_mmt.h>

double 
wf_MMT( const gsl_vector *v,void * params)
{
  double k, w, N, p, l, L, T, Proba;
  
  l = gsl_vector_get(v,0);  
  p = gsl_vector_get(v,1);

  k = (( wf_params *) params)->k;
  w = (( wf_params *) params)->w;
  N = (( wf_params *) params)->N;
  
  //We verify the parameters before evaluation
  if ( p<0 || p > w || l<p || l-p> 1-k-w )
    return 1;
 
  L=0.5*binomial(k+l,p/2);  
  T = max(3, L, 2*L-p, 4*L-p-l);
  Proba = binomial(1-k-l,w-p) + binomial(k+l,p) - entropy(w) ;
  if (entropy(w)<1-k)
    return -1*Proba + T;	
  return -1*min(N+Proba,0) + T;	  
}

double 
Optimal_wf_MMT(wf_params *params)
{
  double k, w, N, l, p;  

  int iter=0;
  const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
  double size, wf;

  gsl_multimin_fminimizer *s=NULL;
  gsl_multimin_function F;
  gsl_vector *ss,*x;

  k = params -> k;
  w = params -> w;
 
  //We initialize the parameters for minimization

  p = 0.3 * w;
  l = 0.5 * (1-k-w) + p;
  
  //We choose positon vector for minimization
  x = gsl_vector_alloc(2);
  gsl_vector_set(x,0,l);
  gsl_vector_set(x,1,p);

  //We choose the step size for minimization
  ss = gsl_vector_alloc(2);
  gsl_vector_set(ss, 0, 0.1*l);
  gsl_vector_set(ss, 1, 0.1*p);
  
  //We define function and its parameters 
  F.n = 2;
  F.f = &wf_MMT;
  F.params = params;
  
  //We initialize the mimimizer
  s = gsl_multimin_fminimizer_alloc(T,2);
  gsl_multimin_fminimizer_set(s,&F,x,ss);
  do 
    {
     iter++;
     //We iterate the minimization procedure
     gsl_multimin_fminimizer_iterate(s);
     //We obtain the characteristic size of minimization procedure	
     size=gsl_multimin_fminimizer_size(s);  
     }while (gsl_multimin_test_size(size,1e-10) == GSL_CONTINUE && iter<700);

  //We save the optimal parameters
  params->wf = s->fval;  	  	
  params->l = gsl_vector_get(s->x,0);
  params->p = gsl_vector_get(s->x,1);

  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(ss);
  gsl_vector_free(x); 
 
  return params->wf;
}
