/* include/entropy_tools.h
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CaWoF.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file include/entropy_tools.h
 * \brief Auxiliary library for calculation of Work Factor.
 *
 * This Library contains the structure wf_params which comunicate several
 * parameters in the optimization of several work factors. 
 */

#ifndef ENTROPY_TOOLS_H
#define ENTROPY_TOOLS_H

#include<stdlib.h>
#include<math.h>
#include<stdarg.h>

#include<gsl/gsl_errno.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_roots.h>
#include<gsl/gsl_multimin.h>

/**
 * \struct wf_params
 * \brief Parameters of the decoding algorithms.
 *
 * One variable of this type will be use in any function about work factor to 
 * transfer the parameters of decoding algorithms. After the optimization of a
 * function this struct saves the optimal parameters to achieve the Work Factor.
 */

typedef struct {
	double k, w, N, l, p, p1, p2, a, wf;
	} wf_params;

/**
 * \fn min(double x, double y)
 * \brief Minimum of two double values.
 */
double
min(double x, double y);

/**
 * \fn max(int num, ...)
 * \brief Maximum of num double values.
 */
double
max (int num, ...);

/**
 * \fn entropy(double x)
 * \brief This function calculates the (binary) entropy.
 */
double 
entropy (double x);

/**
 * \fn dif_entropy(double x,void * params)
 * \brief The diference entropy(x)-params.
 */
double 
dif_entropy (double x, void * params);

/**
 * \fn binomial(double x, double y)
 * \brief This function approximates Binomial coefficient by entropy function.
 */
double 
binomial(double x, double y);

/**
 * \fn entropy_df(double x,void * params)
 * \brief Derivate of entropy function.
 */
double 
entropy_df(double x, void *params);

/**
 * \fn entropy_fdf(double x, void * params, double *y, double *dy)
 * \brief Obtains the dif_entropy function and its derivate. 
 */	
void 
entropy_fdf(double x, void * params, double *y, double *dy);

/**
 * \fn entropy_inverse(double y)
 * \brief Obtains the inverse of entropy function and return a real in [0,0.5].
 */
double 
entropy_inverse(double y);
#endif  /* ENTROPY_TOOLS_H*/
