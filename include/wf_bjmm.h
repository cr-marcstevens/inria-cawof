/* include/wf_bjmm.h
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
 * along with CaWoF. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file include/wf_bjmm.h
 * \brief Library for calculation of BJMM's Work Factor.
 */

#ifndef WF_BJMM_H
#define WF_BJMM_H

#include"entropy_tools.h"

/**
 * \fn double wf_bjmm(const gsl_vector *v, void *paramslp)
 * \brief Work Factor of BJMM's algorithm. 
 *
 * Parameters p1 and p2 in gsl_vector v and k, w, N, l ,p in paramslp.
 */
double
wf_BJMM( const gsl_vector *v,void *paramslp);

/**
 * \fn double Optimal_wf_BJMM_p12(const gsl_vector *vv, void *params)
 * \brief Optimal Work Factor of BJMM's algorithm for l and p fixed.
 * 
 * Parameters l and p in vv and k, w, N in params. It saves Optimal Parameters
 * p1 and p2 in paramslp. 
 */
double 
Optimal_wf_BJMM_p12(const gsl_vector *vv, void * params);

/**
 * \fn double Optimal_wf_BJMM(wf_params *params)
 * \brief Optimal Work Factor of BJMM's algorithm.
 * 
 * Parameters k, w, N in params. It saves Optimal Parameters l, p, p1 and p2 
 * in params. 
 */
double
Optimal_wf_BJMM(wf_params *params);
#endif /*WF_BJMM_H*/
