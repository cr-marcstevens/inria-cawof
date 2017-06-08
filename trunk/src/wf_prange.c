/* src/wf_prange.c
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
 * \file src/wf_prange.c
 * \brief Implementation of include/wf_prange.h
 */

#include <wf_prange.h>

double
wf_Prange( wf_params* params)
{
  double k,w,N,Proba;

  k = params->k;
  w = params->w;
  N = params->N;
 
  Proba =  binomial(1-k, w) - entropy(w) ;

  if (entropy(w)< 1-k) 
    params->wf = -Proba;
  else
    params->wf = -N-Proba;
  return params->wf;
}


