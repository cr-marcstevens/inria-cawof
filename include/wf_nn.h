/* include/wf_nn.h
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
 * \file include/wf_nn.h
 * \brief Library for calculation of Nearest Neighboord's Work Factor.
 */

#ifndef WF_NN_H
#define WF_NN_H

#include"entropy_tools.h"

/**
 * \fn double wf_NN(double p2, void *paramslpp1)
 * \brief Work Factor of Nearest Neighbord's algorithm. 
 *
 * Parameters k, w, N, l ,p, p1 in paramslpp1.
 */
double
wf_NN( double p2, void *paramslpp1);

/**
 * \fn double Optimal_wf_NN_p(const gsl_vector *v, void *params)
 * \brief Optimal Work Factor of Nearest Neighbord's algorithm for l, p fixed.
 * 
 * Parameters l and p in v and k, w, N in params. It saves Optimal Parameters
 * p2 in paramslp ( p1 is chosen from l and p). 
 */

double 
Optimal_wf_NN_p(const gsl_vector *v, void * params);

/**
 * \fn double Optimal_wf_NN(wf_params *params)
 * \brief Optimal Work Factor of Nearest Neighbord's algorithm.
 * 
 * Parameters k, w, N in params. It saves Optimal Parameters l, p, p1 and p2
 * in params.
 */
double
Optimal_wf_NN(wf_params *params);
#endif /*WF_NN_H*/
