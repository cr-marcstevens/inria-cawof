/* examples/behaviour.c
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
 * \file examples/behaviour.c
 * \brief Example program, it shows wf/w vs w. 
 *
 * For a code rate k, this program calculates wf/w respects the error rate w
 * ( with values between 0 and Gilbert-Varshamov bound). We can see the wf/w 
 * becomes near to expect value -log_2(1-k) when w approachs 0.
 */

#include <entropy_tools.h>
#include <wf_prange.h>
#include <wf_stern.h>
#include <wf_dumer.h>
#include <wf_mmt.h>
#include <wf_bjmm.h>
#include <wf_nn.h>

#include <stdio.h>
#include <stdbool.h>
#include <string.h>

FILE *streamout;
wf_params wfp;
double (* wf_algo)(wf_params *);
char algo[20];

/**
 * \fn static void usage ( void )
 * \brief Prints the help of program.
 *
 */
void
usage( void )
{
  printf ("Usage: example ALGORITHM CODE_RATE FILE...\n"
          "Saves Work Factor vs Error Rate ( for a constant Code Rate). \n"
          "Available algorithms are PRANGE, STERN, DUMER, MMT, BJMM, NN.\n\n" 
          );   
}  

int
main (int argc, char* argv[])
{

  double wf;


  if (argc != 3) 
    usage();
 
  switch (argv[0][0])
     {
          case 'P':
            if (strcmp(argv[0],"PRANGE ") )
               wf_algo=&wf_Prange;
            else 
              usage ( ); 
            break;
          case 'S':
            if (strcmp(argv[0],"STERN ") )
               wf_algo=&Optimal_wf_Stern;
            else 
              usage ( ); 
            break;
          case 'D':
            if (strcmp(argv[0],"DUMER ") )
               wf_algo=&Optimal_wf_Dumer;
            else 
              usage ( ); 
            break;
          case 'M':
            if (strcmp(argv[0],"MMT ") )
               wf_algo=&Optimal_wf_MMT;
            else 
              usage ( ); 
            break;
          case 'B':
            if (strcmp(argv[0],"BJMM ") )
               wf_algo=&Optimal_wf_BJMM;
            else 
              usage ( ); 
            break;
          case 'N':
            if (strcmp(argv[0],"NN ") )
               wf_algo=&Optimal_wf_NN;
            else 
              usage ( ); 
            break;
          default:     
            usage ( ); 
           }
  
  wfp.k = atof (argv[1]);
  if ( (wfp.k<=0)||(wfp.k>=1) )
    {
      printf ("Code rate '%.10f' is illegal \n", wfp.k);
      usage( );
    }   
  
  if ((streamout = fopen (argv[2], "w")) == NULL)
  {
   printf ("We could not find or create the file '%s'\n", argv[2]);
   usage( );
  } 
 
  for (int i=1; i<=200; i++)
    {
      wfp.w=i*0.005*entropy_inverse(1-wfp.k);
      wf = (*wf_algo)(& wfp);
      fprintf(streamout, "%.11f %.11f %.11f %.11f %.11f %.11f %.11f \n",
              wfp.k, wfp.w, wf/wfp.w, wfp.l, wfp.p, wfp.p1, wfp.p2);
     }  
  fclose (streamout); 
}
