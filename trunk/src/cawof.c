/* src/cawof.c
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
 * \file src/cawof.c
 * \brief Main program of calculation of Work Factor.
 *
 * Calculates of Work Factor of an ISD algorithm in the code rate and error 
 * rate. Default algorithm is BJMM, default code rate is 0.5 and defaut error 
 * rate is s Gilbert-Varshamov's bound. Available algorithms are PRANGE, STERN,
 * DUMER, MMT, BJMM and NN. 
 */

#include "cawof.h"
#include <entropy_tools.h>
#include <wf_prange.h>
#include <wf_stern.h>
#include <wf_dumer.h>
#include <wf_mmt.h>
#include <wf_bjmm.h>
#include <wf_nn.h>

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include <getopt.h>

bool verbose, quiet, format;
wf_params wfp;
char algo[20];
double (* wf_algo)(wf_params *);

/**
 * \fn static void version (void)
 * \brief Gives the version, subversion and revision of program.
 *
 */
void
version (void)
{
  printf ("%s %d.%d.%d.\nThis software calculates the asymptotic complexity for"
          "decoding algorithms.\n",PROG_NAME, PROG_VERSION, PROG_SUBVERSION, 
          PROG_REVISION);
}

/**
 * \fn static void usage (int status)
 * \brief Print the help of program.
 */
void
usage (int status)
{
  if (status == EXIT_SUCCESS)
    {
      printf("Usage: cawofa [OPTION] ... \n"
	     "Calculates of Work Factor of an algorithm in the a code rate. \n"
             "Available algorithms are PRANGE, STERN, DUMER, MMT, BJMM, NN.\n"
             "Default values code rate =0.5, error rate = GV's bound,\n" 
             "algorithm = BJMM. \n\n"
	     "-v,  --verbose          verbose output\n"
	     "-q,  --quiet            quiet output\n"
             "-k,  --dimension        code rate input\n"
             "-a,  --algorithm        algorithm input\n"
	     "-w,  --errors           error rate input\n"	
             "-f,  --format           output format respect to error rate\n"
             "-V,  --version          display version and exit\n"
	     "-h,  --help             display this help\n"              
             );
    }
  else
    {
      fprintf (stderr, "Try 'cawofa --help' for more information.\n");
      exit (EXIT_FAILURE);
    }
}

int
main (int argc, char* argv[])
{

  int optc;
  /*We declared the set of  options of the program */
  static struct option long_opts[] = 
 	 { {"verbose", no_argument, NULL, 'v'},
	   {"quiet", no_argument, NULL, 'q'},
           {"dimension", required_argument, NULL, 'k'},
           {"algorithm", required_argument, NULL, 'a'},
           {"errors", required_argument, NULL, 'w'},  
 	   {"format", no_argument, NULL, 'f'},
           {"version", no_argument, NULL, 'V'},
 	   {"help", no_argument, NULL, 'h'},
  	   {NULL, 0, NULL, 0}
	  };
  /*We use the function getopt_long in a loop for reading all the differents
    options of cawofa program */

  verbose = false;
  quiet = false;
  format = false;
  strcpy(algo, "BJMM");
  wf_algo=&Optimal_wf_BJMM;
  wfp.k = 0.5;
  wfp.w = 0;
  wfp.N = 0;
  wfp.l = -1;
  wfp.p = -1;
  wfp.p1 = -1;
  wfp.p2 = -1;

  while ((optc =
	  getopt_long (argc, argv, "vqk:a:w:fVh", long_opts, NULL)) != -1)
    switch (optc)
      {
      case 'v':
	verbose = true;
	break;
      case 'q':
	quiet = true;
	break;
      case 'k':
        wfp.k = atof (optarg);
        if ( (wfp.k<=0)||(wfp.k>=1) )
          {
            fprintf (stderr, "Code rate '%.10f' is illegal \n", wfp.k);
            usage (EXIT_FAILURE);
           }   
        break;
      case 'a':
        switch (optarg[0])
        {
          case 'P':
            if (strcmp(optarg,"PRANGE ") )
              {
               strcpy(algo, "Prange's");
               wf_algo=&wf_Prange;
              }
              
            else 
              usage (EXIT_FAILURE); 
            break;
          case 'S':
            if (strcmp(optarg,"STERN ") )
              {
               strcpy(algo, "Stern's");
               wf_algo=&Optimal_wf_Stern;
              }
            else 
              usage (EXIT_FAILURE); 
            break;
          case 'D':
            if (strcmp(optarg,"DUMER ") )
              {
               strcpy(algo, "Dumer's");
               wf_algo=&Optimal_wf_Dumer;
              }
            else 
              usage (EXIT_FAILURE); 
            break;
          case 'M':
            if (strcmp(optarg,"MMT ") )
              {
               strcpy(algo, "MMT");
               wf_algo=&Optimal_wf_MMT;
              }
            else 
              usage (EXIT_FAILURE); 
            break;
          case 'B':
            if (strcmp(optarg,"BJMM ") )
              {
               strcpy(algo, "BJMM");
               wf_algo=&Optimal_wf_BJMM;
              }
            else 
              usage (EXIT_FAILURE); 
            break;
          case 'N':
            if (strcmp(optarg,"NN ") )
              {
               strcpy(algo, "Nearest Neighbords ");
               wf_algo=&Optimal_wf_NN;
              }
            else 
              usage (EXIT_FAILURE); 
            break;
          default:     
            usage (EXIT_FAILURE); 
           }
        break;
      case 'w':
	wfp.w = atof (optarg);
        if ( (wfp.w<0) || (wfp.w>1) )
	  {
	    fprintf (stderr, "CaWoFa: error: Error rate must be in ]0,1[ "
		     "'%.10f'\n", wfp.w);
	    usage (EXIT_FAILURE);
	  }
	break;
      case 'N':
	wfp.N = atof (optarg);
        if ( (wfp.w<0))
	  {
	    fprintf (stderr, "CaWoFa: error: Number of solutions no negative "
		     "'%.10f'\n", wfp.N);
	    usage (EXIT_FAILURE);
	  }
	break;
      case 'f':
	format = true;
	break;	
      case 'V':
	version ();
	exit (EXIT_SUCCESS);
	break;
      case 'h':
	usage (EXIT_SUCCESS);
	exit (EXIT_SUCCESS);
	break;
      default:
        fprintf(stderr, "CaWoFa: error: wrong optional argument\n");
	usage (EXIT_FAILURE);
      }


  /*We verify if the code rate is acceptable*/
  if ( wfp.w == 0)
    wfp.w = entropy_inverse(1-wfp.k); 
  else if ( entropy(wfp.w) >1-wfp.k  )
      {
       if (verbose)
          fprintf (stderr, 
                   "CaWoFa: warning: error rate is larger than GV's bound\n");
       if (wfp.N==0)
          wfp.N=entropy(wfp.w) - (1-wfp.k);  
      }  
      
   
  if (wfp.w>(1-wfp.k)/2 ) 
    {
    fprintf (stderr, "CaWoFa: error rate is too large\n");   
    return 0;
       }

  if (!quiet)
    printf ("The work factor of %s algorithm is assymptotically 2^(",algo);
  if (format)
    printf ("%.10f", (*wf_algo)(& wfp)/wfp.w);
  else
    printf ("%.10f", (*wf_algo)(& wfp));      
  
  if ( !quiet )
    if (format)
      printf("w), when the code rate is");     
    else 
      printf("n), when the code rate is");
 
  printf (" %.10f ",wfp.k);
  if (!quiet)
    printf("and error rate is ");
  printf ("%.10f",wfp.w);  
  if (verbose && (wf_algo!=&wf_Prange))
    {
      if (!quiet)
        printf("\nThe optimal parameters are");
      if(wfp.l>=0)
        if (!quiet) 
          printf(" l=%.10f", wfp.l);
        else 
          printf(" %.10f ",wfp.l);       
      if(wfp.p>=0)
        if (!quiet) 
          printf(" p=%.10f", wfp.p);
        else 
          printf(" %.10f ",wfp.p);
      if(wfp.p1>=0)
        if (!quiet) 
          printf(" p1=%.10f", wfp.p1);
        else 
          printf(" %.10f ",wfp.p1);
      if(wfp.p2>=0)
        if (!quiet) 
          printf(" p2=%.10f", wfp.p2);
        else 
          printf(" %.10f ",wfp.p2);
     }   
  printf("\n");       
}
