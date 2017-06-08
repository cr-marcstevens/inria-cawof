/* src/entropy_tools.c
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
 * \file src/entropy_tools.c
 * \brief Implementation of include/entropy_tools.h
 */

#include<entropy_tools.h>

double max(int num, ...)
{

 va_list valist;
 double maximum,aux; 
 int i;

 va_start(valist,num);
 maximum=va_arg(valist,double);

 for(i=1; i<num; i++)
    {
     aux=va_arg(valist, double);
     if (maximum<aux)
        maximum=aux;
    } 
 va_end(valist);   

 return maximum; 
}

double min(double x, double y)
{
 if (x<y)
   return x;
 return y;
}

double 
entropy (double x)
{
  if ((x==0)||(x==1))
    return 0;
  return -(x*log(x) + (1-x)*log(1-x))/log(2);
}

double 
binomial(double x, double y)
{
  if (x==0)
    return 0;
  return x*entropy(y/x);
}

double 
dif_entropy(double x, void * params) 
{
  double a = *((double *) params);
  if ((x == 0) || (x == 1))
    return - a;
  return -(x * log(x) + (1 - x) * log(1 - x)) / log(2) - a;
}

double 
entropy_df(double x, void * params)
{
  return -(log(x/(1-x))/log(2) );
}

void 
entropy_fdf(double x, void * params, double *y, double *dy)
{
  double a = *((double *) params);
  *y = - (x * log(x) + (1 - x) * log(1 - x)) / log(2) - a;
  *dy = log((1 - x) / x) / log(2);
}

double 
entropy_inverse(double y)
{
  //Number of iterations for searching the inverse
  int iter=0; 
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double x,x0=y/2, yy;

  gsl_function_fdf FDF;
  if (y==0)
    return 0;
  else 
    if (y==1)
      return 0.5;

  //We initialize the function for inversing procedure
  FDF.f=&dif_entropy;
  FDF.df=&entropy_df;
  FDF.fdf=&entropy_fdf;

  //We transform the value a to a parameter aa
  yy=y;
  FDF.params=&yy;

  //We choose the Steffenson method for inversing
  T=gsl_root_fdfsolver_steffenson;

  //We initialize the solver with the chose method
  s=gsl_root_fdfsolver_alloc(T);

  /*We initialize the the solver s to use the function and its derivative and 
   the initial guess root x0.*/
  gsl_root_fdfsolver_set(s,&FDF,x0);
    do
      {
        iter++;
	gsl_root_fdfsolver_iterate(s);
          //We actualize the new guess root x0
	x=x0;
	x0=gsl_root_fdfsolver_root(s);   		
	} while (gsl_root_test_delta(x0,x,0,1e-10) == GSL_CONTINUE && iter<100);

  gsl_root_fdfsolver_free(s);	
  return x0;		
}
