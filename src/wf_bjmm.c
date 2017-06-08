/* src/wf_bjmm.c
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
 * \file src/wf_bjmm.c
 * \brief Implementation of include/wf_bjmm.h
 */

#include<wf_bjmm.h>

double
wf_BJMM ( const gsl_vector *v, void * paramslp)
{
  double k, w, N, l, p, p1, p2, L0, L1, L2, mu1, mu2, Proba, T;

  k = (( wf_params *) paramslp)->k;
  w = (( wf_params *) paramslp)->w;
  N = (( wf_params *) paramslp)->N;
  l = (( wf_params *) paramslp)->l;
  p = (( wf_params *) paramslp)->p;

  p1 = gsl_vector_get(v, 0);
  p2 = gsl_vector_get(v, 1);

  if ( p1 < 0 || p1 < 0.5*p || p1 > k+l-0.5*p || p2 < 0.5*p1 || p2 > k+l-0.5*p1)
    return 1;		
  
  L2 = binomial(k+l,p2);
  L1 = binomial(k+l,p1);
  L0 = binomial(k+l,p);
  mu1 = binomial(p1, 0.5*p) + binomial(k+l-p1, 0.5*p) - binomial(k+l, p1); 
  mu2 = binomial(p2, 0.5*p1) + binomial(k+l-p2, 0.5*p1) - binomial(k+l, p2); 
  
  T = max(4,0.5*L2, L1-mu2-L2, L0 -mu1-mu2-L1, L0-mu1-l);
  Proba = binomial(1-k-l,w-p) + binomial(k+l,p) - entropy(w) ;
 
  if (entropy(w)<1-k)
    return -1*Proba + T;	
  return -1*min(N+Proba,0) + T;	    
  }

double 
Optimal_wf_BJMM_p12(const gsl_vector *vv, void * params)
{
  wf_params paramslp;
  int iterp12=0;
  const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
  double k, w, N, l, p, sizep12, wf;
  gsl_multimin_fminimizer *s=NULL;
  gsl_multimin_function F_p12;
  gsl_vector *ss,*y;

  k = (( wf_params *) params)->k;
  w = (( wf_params *) params)->w;
  N = (( wf_params *) params)->N;

  l=gsl_vector_get(vv,0);
  p=gsl_vector_get(vv,1);	

  //We verify the parametres before evaluation		
  if ( l < 0 || w-p > 1-k-l || p > w || p > k+l || p < 0 )
    return 1;
 
  //We save the parameters l, p for the function wf_MMT
  paramslp.k = k;
  paramslp.w = w;
  paramslp.N = N;
  paramslp.l = l;
  paramslp.p = p;

  //We choose positon vector for minimization
  y = gsl_vector_alloc(2);
  gsl_vector_set(y,0,0.75*p);
  gsl_vector_set(y,1,0.75*gsl_vector_get(y,0));	

  //We choose the step size for minimization
  ss = gsl_vector_alloc(2);
  gsl_vector_set(ss,0,0.2*gsl_vector_get(y,0));
  gsl_vector_set(ss,1,0.2*gsl_vector_get(y,1));

  //We define function and its parameters 	
  F_p12.n = 2;
  F_p12.f = wf_BJMM;
  F_p12.params = &paramslp;

  //We initialize the mimimizer
  s = gsl_multimin_fminimizer_alloc(T,2);
  gsl_multimin_fminimizer_set(s,&F_p12,y,ss);
// printf("wf=%f l =%f p =%f p1=%f p2=%f\n",    s->fval,  paramslp.l,  paramslp.p, gsl_vector_get(s->x,0),  gsl_vector_get(s->x,1) );
  do 
  {

    iterp12++;
   //We iterate the minimization procedure
    gsl_multimin_fminimizer_iterate(s);     
  
   //We obtain the characteristic size of minimization procedure
    sizep12=gsl_multimin_fminimizer_size(s);
   }while(gsl_multimin_test_size(sizep12,1e-12) == GSL_CONTINUE && iterp12<500);
  
  wf = s->fval;
  //We save the optimal values of p1 and p2
  (( wf_params *) params)->p1 = gsl_vector_get(s->x,0);
  (( wf_params *) params)->p2 = gsl_vector_get(s->x,1);
  //printf("wf=%f l =%f p =%f p1=%f p2=%f\n",    s->fval,  paramslp.l,  paramslp.p, gsl_vector_get(s->x,0),  gsl_vector_get(s->x,1) );

  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(ss);
  gsl_vector_free(y);
		


  return wf;
}

double 
Optimal_wf_BJMM(wf_params *params)
{
  double k, w, N, l, p, wf;  

  int iter=0;
  const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
  double size, W;

  gsl_multimin_fminimizer *s=NULL;
  gsl_multimin_function F;
  gsl_vector *ss,*x;

  k = params->k;
//  w = params->w;
  N = params->N;

  F.n=2;
  F.f = &Optimal_wf_BJMM_p12;
 

  W = params->w;
  params->w=entropy_inverse(1-k)/20;
  
  if ((params->w)/2 < W)
     {
      params->l = 0.75*(1-k-W);
      params->p = 0.75*min(k+params->l,W);
         }
  else  
    {
     params->l = 0.75*(1-k-0.5*params->w);
     params->p = 0.75*min(k + params->l, 0.5*(params->w));
     }   
      
  x=gsl_vector_alloc(2);
  ss=gsl_vector_alloc(2);

  do
    {
     iter =0;
     if ((params->w)/2 < W)
       {
        params->w = W;
       }  
     else  
       {
        params->w = (params->w)/2;
       }

     gsl_vector_set(x,0,params->l);
     gsl_vector_set(x,1,params->p);

     gsl_vector_set(ss,0,0.2*params->l);
     gsl_vector_set(ss,1,0.25*params->p);


     F.params = params;
     s = gsl_multimin_fminimizer_alloc(T,2);
     gsl_multimin_fminimizer_set(s,&F,x,ss);

     do 
       {
        iter++;
        //We iterate the minimization procedure
        gsl_multimin_fminimizer_iterate(s);
        //We obtain the characteristic size of minimization procedure	
        size=gsl_multimin_fminimizer_size(s);    
       } while ( gsl_multimin_test_size(size,1e-10) == GSL_CONTINUE 
                       && iter<100);
     params->l = gsl_vector_get(s->x,0);
     params->p = gsl_vector_get(s->x,1);
    
    } while (params->w > W);    

  wf=s->fval;  

  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(ss);
  gsl_vector_free(x);

  return wf;
}
