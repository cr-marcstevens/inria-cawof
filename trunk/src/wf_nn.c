/* src/wf_nn.c
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
 * \file src/wf_nn.c
 * \brief Implementation of include/wf_nn.h
 */

#include<wf_nn.h>

double
wf_NN (double p2, void * paramslpp1) 
{
  double k, w, N, l, p, l2, L, mu, gamma, p1, Proba, y,T; 

  k = (( wf_params *) paramslpp1)->k;
  w = (( wf_params *) paramslpp1)->w;
  N = (( wf_params *) paramslpp1)->N;
  l = (( wf_params *) paramslpp1)->l;
  p = (( wf_params *) paramslpp1)->p;
  p1 = (( wf_params *) paramslpp1)->p1;

  if ( p1/2 > p2  ||  p2 > k+l - p1/2 )
     return 1;

  L = 0.5* binomial ( k + l , p2 ) ;
  l2 =  p1 + binomial(k+ l - p1, p2 - 0.5*p1);
  mu = binomial(k+l, p1) - l;
  gamma = (w - p) / (1- k - l);
  y = ( 1 - gamma) * (1- entropy( ( entropy_inverse (1 - mu/(1-k-l) ) 
       - gamma/2 )/(1-gamma) ) );

  Proba = binomial(1-k-l,w-p) + binomial(k+l,p) - entropy(w) ;
  T = max(5, L, 2*L-l2, 4*L-l-l2, mu, y*(1 - k - l) );
 
  if (entropy(w)<1-k)
    return -1*Proba + T;	
  return -1*min(N+Proba,0) + T;	    


}


double 
Optimal_wf_NN_p(const gsl_vector *v, void * params)
{
 
  wf_params paramslpp1;
  double k,w,N,l, p, p1, p2, mu;
  double a,b,wf,aux;

  k = (( wf_params *) params)->k;
  w = (( wf_params *) params)->w;
  N = (( wf_params *) params)->N;

  l=gsl_vector_get(v,0);
  p=gsl_vector_get(v,1);	

  if ( 0>=p ||  p >= w || p>=k+l || l<=p || w-p >= 1-k-l ) 
    return 1;

  p1 = 0.5*p + (k+l-p)*entropy_inverse((l-p)/(k+l-p) );

  if ( p1 > k+l )
      return 1;

   mu = binomial(k + l, p1) - l;

  if ( mu<=1e-10 || mu>=(1-k-l)*( 1 - entropy( (w-p)/(1-k-l) ) ) ) 
     return 1;
 
  paramslpp1.k = k;
  paramslpp1.w = w;
  paramslpp1.N = N;
  paramslpp1.l = l;
  paramslpp1.p = p;
  paramslpp1.p1 = p1;

  a= p1/2;

  // It will be enough to look for small values, this comes from the the fact 
  // that l2<l

  b = p1/2+0.75*entropy_inverse((l-p1)/(k+l-p1))*(k+l-p1);
  wf = wf_NN(a, &paramslpp1);
  p2 = a;

  for (int i=1;i<400;i++)
  {
    aux=wf_NN(a+i*(b-a)/400,&paramslpp1);
    if(wf>aux)
    {
     p2=a+i*(b-a)/400;
     wf=aux; 
      }
   }
  
  //We save the optimal values of p1 and p2
  (( wf_params *) params)->p1 = p1;
  (( wf_params *) params)->p2 = p2;
  
  return wf;
}

double 
Optimal_wf_NN( wf_params * params)
{
  double k, w, N, l, p;  

  int iter=0;
  const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
  double size, wf, W;

  gsl_multimin_fminimizer *s=NULL;
  gsl_multimin_function F;
  gsl_vector *ss,*x;

  k = (( wf_params *) params)->k;
  w = (( wf_params *) params)->w;
  N = (( wf_params *) params)->N;

  p = 0.15*w;
  l = 0.10*(1-k-w) + p;

  x=gsl_vector_alloc(2);
  gsl_vector_set(x,0,l);
  gsl_vector_set(x,1,p);

  ss=gsl_vector_alloc(2);
  gsl_vector_set(ss,0,0.1*l);
  gsl_vector_set(ss,1,0.1*p);

  F.n=2;
  F.f=&Optimal_wf_NN_p;
  F.params=params;

  s = gsl_multimin_fminimizer_alloc(T,2);
  gsl_multimin_fminimizer_set(s,&F,x,ss);
      
  do 
    {
      iter++;
      //We iterate the minimization procedure
      gsl_multimin_fminimizer_iterate(s);
      //We obtain the characteristic size of minimization procedure	
      size=gsl_multimin_fminimizer_size(s);    
     } while ( gsl_multimin_test_size(size,1e-10) == GSL_CONTINUE && iter<700);
   
  l = gsl_vector_get(s->x,0);
  p = gsl_vector_get(s->x,1);

  iter=0;
      
  gsl_vector_set(x,0,l);
  gsl_vector_set(x,1,p);

  gsl_vector_set(ss,0,0.05*l);
  gsl_vector_set(ss,1,0.05*p);

  s = gsl_multimin_fminimizer_alloc(T,2);
  gsl_multimin_fminimizer_set(s,&F,x,ss);
   
  do 
    {
     iter++;
     //We iterate the minimization procedure
     gsl_multimin_fminimizer_iterate(s);
     //We obtain the characteristic size of minimization procedure	
     size=gsl_multimin_fminimizer_size(s);    
     } while (gsl_multimin_test_size(size,1e-10) == GSL_CONTINUE && iter<700);

  params->l = gsl_vector_get(s->x,0);
  params->p = gsl_vector_get(s->x,1);

      
  wf=s->fval; 

  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(ss);
  gsl_vector_free(x);	   
 
  return wf;
}
