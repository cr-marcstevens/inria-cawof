/* examples/bound.c
 * 
 * Copyright (C) 2015, 2016 Rodolfo Canto Torres 
 *
 * This file is part of CaWoFa.
 *
 * CaWoFa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3algo of the License, or
 * (at your option) any later version.
 *
 * CaWoFa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CaWoFa.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file examples/bound.c
 * \brief Auxilairy program which calculates the index a of the bound fonction. 
 *
 * This program takes a generic decoding algorithm algo as input and returns
 * a in ]0,1[ such that min B_a = WF(algo) for many code rate k in [0,1] ( and 
 * error rate w as Gilbert Varshamov's bound ).
 */
#include <entropy_tools.h>
#include <bound_function.h>

#include <wf_prange.h>
#include <wf_stern.h>
#include <wf_dumer.h>
#include <wf_mmt.h>
#include <wf_bjmm.h>
#include <wf_nn.h>


#include <stdio.h>

#include <getopt.h>

FILE *streamout;		/* In the case we want to WRITE in a file */
wf_params wfp;			/* Parameters of the code */
double (* wf_algo)(wf_params *);

/**
 * \fn static void usage (int status)
 * \brief Print the help of program.
 */

void
usage (int status)
{
  if (status == EXIT_SUCCESS)
    {
      fprintf (streamout, "Usage: bound [OPTION] algo \n"
	       "Calculates              asymptotic complexity decoding algorithms: \n"
               "0 :                     Prange's Algorithm\n"
               "1 :                     Stern's Algorithm\n"
               "2 :                     Dummer's Algorithm\n"
               "3 :                     MMT's Algorithm\n"
               "4 :                     BJMM's Algorithm\n"
               "5 :                     Nearest Neighbor Algorithm\n"
	       "-o,  --output=FILE      write result to FILE\n"	  
	       "-h,  --help             display this help\n");
    }
  else
    {
      fprintf (stderr, "Try 'bound --help' for more information.\n");
      exit (EXIT_FAILURE);
    }
}

/**
 * \fn int main (int argc, char* argv[])
 * \brief Calculates the index a of the bound fonction for an algorithm.
 */

int
main (int argc, char* argv[])
{

  streamout=stdout;		/*We write in stdout for default */
  double wf;
  int optc;

  /*We declared the set of  options of the program */
  static struct option long_opts[] = 
 	 { {"output", required_argument, NULL, 'o'},
 	   {"help", no_argument, NULL, 'h'},
  	   {NULL, 0, NULL, 0}
	  };
  /*We use the function getopt_long in a loop for reading all the differents
    options of bound program */
  while ((optc =
	  getopt_long (argc, argv, "o:h", long_opts, NULL)) != -1)
    switch (optc)
      {
      case 'o':
	if ((streamout = fopen (optarg, "w")) == NULL)
	  {
	    fprintf (stderr,
		     "bound: error: We could not find or create the file "
		     "'%s'\n", optarg);
	    usage (EXIT_FAILURE);
	  }
	break;
      case 'h':
	usage (EXIT_SUCCESS);
	exit (EXIT_SUCCESS);
	break;
      default:
	usage (EXIT_FAILURE);
      }

  if (optind == argc)	/* In this case, we have read all the optinal arguments 
			   without the argument code rate */
    {
      fprintf (stderr, "bound: error:no algorithm argument\n");
      usage (EXIT_FAILURE);
    }
  //We verify if the code rate is acceptable.
  switch (atoi(argv[optind]))
    {
    case 0:
      wf_algo = &wf_Prange;
      fprintf(streamout,"Prange's Algorithm test:\n" );
      break;
    case 1:
      wf_algo = &Optimal_wf_Stern;
      fprintf(streamout,"Stern's Algorithm test:\n" );
      break;
    case 2:
      wf_algo = &Optimal_wf_Dumer;
      fprintf(streamout,"Dummer's Algorithm test:\n" );
      break;
    case 3:
      wf_algo = &Optimal_wf_MMT;
      fprintf(streamout,"MMT's Algorithm test:\n" );
      break;
    case 4:
      wf_algo = &Optimal_wf_BJMM;
      fprintf(streamout,"BJMM's Algorithm test:\n" );
      break;         		   
    case 5:
      wf_algo = &Optimal_wf_NN;
      fprintf(streamout,"Nearest Neighbors Algorithm test:\n" );
      break;         		   
    default:
      fprintf (stderr, "Algorithm number '%s' is illegal \n", argv[optind]);
      usage (EXIT_FAILURE);
  }

  for (int i=1; i<100; i++)
    {
      wfp.k = i * 0.01;
      wfp.w = entropy_inverse(1-wfp.k)*0.10;
      wfp.wf = (*wf_algo)( & wfp );
      wfp.a=find_coefficient(&wfp);
      wf = Optimal_reduced_Ba(&wfp);
      fprintf(streamout, "%f  %f %f  %f %f\n", wfp.k, wf, wfp.a, wfp.l, wfp.p   );

     }  
  fclose (streamout);
}
