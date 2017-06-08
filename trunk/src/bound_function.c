/* src/bound_function.c
 * 
 * Copyright (C) 2015, 2016 Rodolfo Canto Torres 
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CaWoF.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file src/bound_function.c
 * \brief Implementation of include/bound_function.h
 */

#include <bound_function.h>

double 
pp(double l, wf_params * params)
{
  double k,a;
  k = params->k;
  a = params->a;

  return (k+l)*entropy_inverse(l/(l+k))/a;
}

double dif_pp(double l, void * params)
{
  return pp(l,params) -  (( wf_params *) params)->w ;
}


double
pp_df(double l, wf_params * params)
{
  double k,a;

  k = params->k;
  a = params->a;

  return ( entropy_inverse(l/(l+k)) + 
           k/( (k+l)*entropy_df(entropy_inverse(l/(l+k)),params)))/a;
}

double
reduced_Ba(double l, wf_params *params)
{
  double k,w;

  k = params->k;
  w = params->w;

  return entropy(w)-(1-k-l)*entropy((w-pp(l,params))/(1-k-l))-l;
}

double
reduced_Ba_df(double l,void * params)
{
  double k,w;

  k = (( wf_params *) params)->k;
  w = (( wf_params *) params)->w;

  double p=pp(l,params), dp=pp_df(l, params);

  return entropy((w-p)/(1-k-l)) - 
    entropy_df((w-p)/(1-k-l),params)*(w-p-dp*(1-k-l))/(1-k-l) - 1;
}


double 
l_max(wf_params *params)
{

  int iter=0;
  const gsl_root_fsolver_type *T;
  double l;
  double c,d;

  double w,a;

  w = params->w;
  a = params->a;

  gsl_root_fsolver *s;
  gsl_function F;
	
  c=0.001*entropy(a*w);
  d=entropy(a*w);

  //We define function and its parameters 
  F.function =& dif_pp;
  
  //We transform the value a to a parameter aa
  F.params=params;

  //We initialize the mimimizer
  T=gsl_root_fsolver_brent;
  s= gsl_root_fsolver_alloc(T);

  //We define the search interval  of minimum
  gsl_root_fsolver_set (s,&F,c, d);
  
  do 
    {
      iter++;
      //We iterate the minimization procedure
      gsl_root_fsolver_iterate(s);
      c=gsl_root_fsolver_x_lower(s);
      d=gsl_root_fsolver_x_upper(s);
    } //We verify the size of the new search interval
      while(gsl_min_test_interval( c,d,0,0.001 ) == GSL_CONTINUE && iter<100);

  l =gsl_root_fsolver_root(s);
  gsl_root_fsolver_free(s);

  return l;
}

double
Optimal_reduced_Ba(wf_params * params)
{
  int iter=0;
  const gsl_root_fsolver_type *T;
  double l, c, d, k; 

  gsl_root_fsolver *s;
  gsl_function F;

  //We take the maximum l that we can calculate reduced_Ba	
  k = params->k;
  
  l=l_max(params);
  if ((1-k )< l )
    l= 1-k; 

  //We define function and its parameters 
  F.function = &reduced_Ba_df;

  //We transform the value a to a parameter aa
  F.params = params;

  c=0.0001*l;
  d=0.999*l;

 //We initialize the mimimizer
  T=gsl_root_fsolver_brent;
  s= gsl_root_fsolver_alloc(T);
  //We define the search interval  of minimum
  gsl_root_fsolver_set (s,&F,c, d);

  do 
    {
      iter++;
      //We iterate the minimization procedure
      gsl_root_fsolver_iterate(s);	
      c=gsl_root_fsolver_x_lower(s);
      d=gsl_root_fsolver_x_upper(s);		
    } //We verify the size of the new search interval
      while(gsl_min_test_interval(c,d,0,0.001 ) == GSL_CONTINUE && iter<100);

  l = gsl_root_fsolver_root(s);

  params->l=l;
  params->p = pp(l,params);
  
  gsl_root_fsolver_free(s);

  return reduced_Ba(l, params);

}

double
dif_Optimal_reduced_Ba(double a, void * params)
{
  double wf_wanted;

  wf_wanted =  (( wf_params *) params)->wf;
  (( wf_params *) params)->a=a;

  return Optimal_reduced_Ba(params) - wf_wanted;
}

double
find_coefficient(wf_params * params)
{
  int iter=0; 
  double coefficient;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;

  //We define function and its parameters 
  F.function = dif_Optimal_reduced_Ba;

  //We transform the value wf to a parameter wwff  
  F.params = params;

  //We initialize the mimimizer
  T=gsl_root_fsolver_brent;
  s=gsl_root_fsolver_alloc(T);

  //Define the search interval 
  gsl_root_fsolver_set(s,&F,0.2,0.8);

  do
    {
      iter++;
      //We iterate the minimization procedure
      gsl_root_fsolver_iterate(s);
    } //We verify the size of the search interval
      while (gsl_root_test_interval(gsl_root_fsolver_x_lower(s), 
             gsl_root_fsolver_x_upper(s),0,1e-4) == GSL_CONTINUE && iter<400 );
   
  coefficient = gsl_root_fsolver_root (s);
  gsl_root_fsolver_free(s);	

  return coefficient;		
}
