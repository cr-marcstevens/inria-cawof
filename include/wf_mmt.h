/* include/wf_mmt.h
 * 
 * Copyright (C) 2015, 2016,2017 Rodolfo Canto Torres 
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
 * \file include/wf_mmt.h
 * \brief Library for calculation of MMT's Work Factor.
 */

#ifndef WF_MMT_H
#define WF_MMT_H

#include"entropy_tools.h"

/**
 * \fn double wf_MMT(const gsl_vector *v, void *params)
 * \brief Work Factor of MMT's algorithm. 
 *
 * Parameters l and p in gsl_vector v and k, w, N in params
 */
double
wf_MMT( const gsl_vector *v,void *params);

/**
 * \fn double Optimal_wf_MMT(wf_params *params)
 * \brief Optimal Work Factor of MMT's algorithm.
 * 
 * Parameters k, w, N in params. It saves Optimal Parameters l and p in params.
 */
double
Optimal_wf_MMT(wf_params *params);


#endif /*WF_MMT_H*/
