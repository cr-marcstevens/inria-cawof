/* include/bound_function.h
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
 * \file include/bound_function.h
 * \brief Library for calculation of Bound Function.
 *
 * This Library calculate theorical Bound Function, its optimal value for a  
 * coefficient a and find the coefficient a for the equation min Ba = wf for 
 * k, w fixed.
 */

#ifndef BOUND_FUNCTION_H
#define BOUND_FUNTION_H

#include"entropy_tools.h"

/**
 * \fn double pp(double l, wf_params * params)
 * \brief Calculates the p such that 2^l=binomial(k+l,ap).
 */
double 
pp(double l, wf_params * params);

/**
 * \fn dif_pp(double l, void * params)
 * \brief the diference pp(l,a)-w.
 */
double 
dif_pp(double l, void * params);

/**
 * \fn pp_df(double l, wf_params * params)
 * \brief Derivate of function pp.
 */
double
pp_df(double l, wf_params * params);

/**
 * \fn double reduced_Ba(double l, wf_params * params)
 * \brief Reduced version of Ba(l,p) 
 *
 * This function will achieve the same minimum value of function Ba(l,p).
 */
double
reduced_Ba(double l, wf_params *params);

/**
 * \fn double reduced_Ba_df(double l, void * params)
 * \brief Derivate of function reduced_Ba.
 */
double
reduced_Ba_df(double l, void * params);

/**
 * \fn double l_max(wf_params *params)
 * \brief Obtains the maximum value of l such Ba is calculable.
 */
double 
l_max(wf_params *params);

/**
 * \fn double Optimal_reduced_Ba(wf_params * params)
 * \brief Obtains the minimum value of  Ba.
 */
double
Optimal_reduced_Ba(wf_params * params);

/**
 * \fn double dif_Optimal_reduced_Ba(double a, void * params)
 * \brief The difference Optimal_reduced_Ba - workfactor.
 */
double
dif_Optimal_reduced_Ba(double a, void * params);

/**
 * \fn double find_coefficient(wf_params * params)
 * \brief Find a such that reduced_Ba_min(a)=wf.
 */
double
find_coefficient(wf_params * params);
#endif /*BOUND_FUNCTION_H*/
