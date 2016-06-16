#include <omp.h>
#include <stdio.h>
//#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define CHUNKSIZE 100
#define num_threads 16
#define threads_max 32
#define num_bins_preproc 10
#define N     1000
#define PI  3.14159621234161928

/*
COMPILE Single.
============================
CC=icc
CFLAGS=-Wall -O3 -openmp

HW3ParallelCorrelationFunc1_0 : HW3ParallelCorrelationFunc1_0.c
   $(CC) -o $@ $< $(CFLAGS)

.PHONY: clean

clean :
   rm -f HW3ParallelCorrelationFunc1_0
============================

RUN
============================
#!/bin/bash

if [[ $# -ne 1 ]] ; then
    echo "usage: ./HW3ParallelCorrelationFunc1_0 num_threads"
    exit 1
fi
export OMP_NUM_THREADS=$1
./HW3ParallelCorrelationFunc1_0
=============================
*/

/*
 * ------------------------------------------------------------------------
 *       ___           _     __  __      _   _
 *      |   \ __ _ _ _| |__ |  \/  |__ _| |_| |_ ___ _ _
 *      | |) / _` | '_| / / | |\/| / _` |  _|  _/ -_) '_|
 *      |___/\__,_|_| |_\_\ |_|  |_\__,_|\__|\__\___|_|
 *
#  .-----..-.  .-. .---.              .-.-.  .---. .-..-. .-..-----.
#  `-' '-'| {  } |/ {-. \     ___     | } }}/ {-. \{ ||  \{ |`-' '-'
#    } {  {  /\  }\ '-} /    {___}    | |-' \ '-} /| }| }\  {  } {
#    `-'  `-'  `-' `---'              `-'    `---' `-'`-' `-'  `-'
#  .----. .---. .---. .---. .----..-.     .--. .-----..-. .---. .-. .-.
#  | }`-'/ {-. \} }}_}} }}_}} |__}} |    / {} \`-' '-'{ |/ {-. \|  \{ |
#  | },-.\ '-} /| } \ | } \ } '__}} '--./  /\  \ } {  | }\ '-} /| }\  {
#  `----' `---' `-'-' `-'-' `----'`----'`-'  `-' `-'  `-' `---' `-' `-'
#  .----..-. .-..-. .-..----..-----..-. .---. .-. .-.
#  } |__}| } { ||  \{ || }`-'`-' '-'{ |/ {-. \|  \  {
#  } '_} \ `-' /| }\  {| },-.  } {  | }\ '-} /| }\  {
#  `--'   `---' `-' `-'`----'  `-'  `-' `---' `-' `-'
#
# This program will make the correlation function given files with [RA, Dec, Z]. 
# It will use the estimator:
#      {   Xi(r) = (N_r/N_d)**2 (DD(r)/RR(r) - 1)   }
#
# The Sample data is separated into 15 logarithmic bins of separation in range 
# 0.1 - 20 h**-1 Mpc ( the first bin at log r = -1 and the last bin is at 
# log r = 1.301 )
#
# The Random Points are redshifts in range 0.02 <= z <= 0.06. Using the SDSS sky
# coverage. 
#
# In 3_3.c adding ability to vary array size and thus read in data on the fly. Malloc.
#
# FILE NAMES:

#      Galaxy Samples:
#      { SDSS_Mr21_rspace.dat, SDSS_Mr20_rspace.dat, SDSS_Mr20_zspace.dat }
#      Random Points:
#      { SDSS_random.dat }

#     DM Samples: Points in Cartesian Coordinates {Xi, Yi, Zi,} Vol=141.3h-1Mpc
#     {DM.dat}
#     Random Points:
#     {DM_random.dat}
#
#     B. Compute the bias function: b(r) = SQRT[ Xi_gal / Xi_DM  ]

------------------------------------------------------------------------
*/

int main (int argc, char **argv)
{
   printf("FIRST line of MAIN\n");
   long int i, j;

   /*
   ---------------------------------------
   READ IN COMMAND LINE ARGUMENTS
   ---------------------------------------
   */
   printf("argc: %d\n",argc);
   if ( argc != 2 ) {
      printf("Usage: bash %s arg1\n", argv[0]);
      printf("where arg1 is the NUM_OF_THREADS (integer <= 32)\n");
      exit(0);
   }

   /*
   ---------------------------------------
   SETS THE NUMBER OF THREADS IF SPECIFIED
   ---------------------------------------
   */
   int k = 1;
   int temp_value = 0;
   int num_threads_set = 0;
   if ( argc > 1 ){
      temp_value = atoi(argv[k]);
      if (temp_value > 0 && temp_value < threads_max){
         num_threads_set = atoi(argv[k]);
      }
   }else {
      num_threads_set = num_threads;
   }

   /*
   ---------------------------------------
       set periodic boundary conditions
   ---------------------------------------
   */

  int periodic_boundary_conditions = 1;  //SET 1 to use PBC.(and analytic random counts) or  SET 0 for !PBC

   /*
   ---------------------------------------
           Initialize Variables
   ---------------------------------------
   */
   int chunk, tid, nthreads;
   //float a[N], b[N], c[N];

   double dec_to_rad = (2. * PI / 360.);
   double deg_to_rad = (2. * PI / 360.);

   //Create the Distance Array.
   //const int num_bins = 15;
   const int num_bins = num_bins_preproc;
   double bin_counts[num_bins_preproc] = {0.0};
   //Comment: The next part was done for convience to prevent renaming everything.
   //The logrNOTSQ_min is the actual logr_min and max values --> r_min, max={.1,20}
   //This trickery is done to avoid having to square_root distance just for an index.
   //double logrNOTSQUARED_min = -1.0;
   //double logrNOTSQUARED_max = 1.3011;
   //MODIFYING ORDER TO MATCH MANODEEPS R-BINS;
   //---------------------
   double r_min = 0.1675;
   double r_max = 5.7885;
   double logrNOTSQUARED_min = log10( r_min );
   double logrNOTSQUARED_max = log10( r_max );
   //double logrNOTSQUARED_max = log10( 4.0618 );
   //---------------------

   //Find the Indexes for log_r if the distance used is the squared distance instead of actual.
   double logr_min = log10( pow( pow(10., logrNOTSQUARED_min), 2.) );
   double logr_max = log10( pow( pow(10., logrNOTSQUARED_max), 2.) );

   long int errorcounts = 0;
   //    r                 0.1,              ...                  , 20
   // logr                 -1,               ...                  ,1.3011
   // index                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
   long int distance_counts[num_bins_preproc] = { 0 };
   long int randdistance_counts[num_bins_preproc] = { 0 };
   double Xi_func[num_bins_preproc] = { 0.0 };
   int dist_index;

   //This should be set below 5495 to take a limited sample.
   //int FILELENGTH = 4000;

   //long int FILELENGTH20r  = 28162;
   //long int FILELENGTH20z  = 28383;
   //long int FILELENGTH21r  = 5495;
   //long int FILELENGTHrand = 42654; 
   //long int FILELENGTHDM     = 50000;

   long int FILELENGTHDM     = 66659;
   //long int FILELENGTHDM     = 10000;
   //long int FILELENGTHDMrand = 100000;

   //long int FILELENGTHDMrand = 250000;
   long int FILELENGTHDMrand = 100000;
   long int N_rand = FILELENGTHDMrand;
   //Filename declaration.
   //string datafilename = 'SDSS_random.dat'
   //static const char random_datafile[] = "SDSS_random.dat";


   //INITIALIZE ARRAYS TO STORE GRID INTEGER VALUES
   //----------------------------------------------
   int box_side_length=1;
   int gridx_LIST[N_rand];
   int gridy_LIST[N_rand];
   int gridz_LIST[N_rand];



   double D, logD, actual_D, logactual_D, actual_Ddist_index;
   double r = 0;
   double x1 = 0;
   double y1 = 0;
   double z1 = 0;

   double rj = 0;
   double x2 = 0;
   double y2 = 0;
   double z2 = 0;


   double delx, dely, delz;
   chunk = CHUNKSIZE;


   /*
   ---------------------------------------
           Grid Data Points
   ---------------------------------------
   */

     //STORE THE GRID INTEGER VALUES FOR EACH POINT
     float grid_binmax = box_side_length;
     float grid_binmin = 0;
     int grid_numbins = 10;

     int grid_x1;
     int grid_y1;
     int grid_z1;
     int grid_x2;
     int grid_y2;
     int grid_z2;
     int gxvar;
     int grid_x_check;
     int grid_y_check;
     int grid_z_check;
     int grid_x_pbc;
     int grid_y_pbc;
     int grid_z_pbc;
     int grid_x_norm;
     int grid_y_norm;
     int grid_z_norm;


    


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /**
      *      ___    _   _  _ ___   ___  __  __
      *     | _ \  /_\ | \| |   \ / _ \|  \/  |
      *     |   / / _ \| .` | |) | (_) | |\/| |
      *     |_|_\/_/ \_\_|\_|___/ \___/|_|  |_|
      */


     // INIITIALIZE FOR RANDOM COUNTS
     //static const char random_datafile[] = "/home/chasonnscr/LargeScaleHW/HW3Codes/DMCorrelationFunctionParallel/DM_random.dat";
     //static const char random_datafile[] = "/scratch/chasonnscr/positdirHPC/PositionXYZdatafile0.dat";
     static const char random_datafile[] = "/scratch/chasonnscr/random/DM_randomResearch.dat";
     //static const char r21_datafile[] = "SDSS_Mr21_rspace.dat";
     // Stores the number of points in random set.


     if (!periodic_boundary_conditions){
       FILE *myrandfile = fopen ( random_datafile, "r" );
       if (myrandfile == NULL) {
          printf("DM_random.dat not opened, trying local filepath...\n");
          static const char random_datafile[] = "./DM_random.dat";
          myrandfile = fopen ( random_datafile, "r" );
          if (myrandfile == NULL) {
            printf("DM_random.dat not opened, exiting...\n");
            exit(0);
          }

       }
       //double * randX_LIST[N_rand] = malloc(N_rand * sizeof(double));
       //double * randY_LIST[N_rand] = malloc(N_rand * sizeof(double));
       //double * randZ_LIST[N_rand] = malloc(N_rand * sizeof(double));
       double randX_LIST[N_rand];
       double randY_LIST[N_rand];
       double randZ_LIST[N_rand];


       //===============
       // READ IN DATA.
       //===============
       for(i = 0; i < (N_rand); ++i){
         fscanf(myrandfile, "%lf", &randX_LIST[i]);
         fscanf(myrandfile, "%lf", &randY_LIST[i]);
         fscanf(myrandfile, "%lf", &randZ_LIST[i]);

         //Find the max value using only X values.
         if (randX_LIST[i]> (float) (box_side_length)) {
              box_side_length = (int) ceil(randX_LIST[i]);
         }

        //printf("BOX SIDE LENGTH: %d \n", box_side_length);

         //printf("i: %ld Data Z: %f \n", i, Z_LIST[i] );
         if (i >= (N_rand-1)){
             printf("Close or exceeded N_data limit. X: %lf \n", randX_LIST[i]);
             //break;
          }
          //++i;
       }

        //CHECK STANDARD GRID_NUMBINS VALUE
       grid_numbins = floor( grid_binmax / r_max );

       printf("Final Grid Bin Resolution: %d Bins......\n", grid_numbins);

       float grid_normalization = grid_numbins/(grid_binmax);
       float grid_binwidth = (grid_binmax/grid_numbins);
       printf("Final Grid Binwidth: %f width [mpc]\n", grid_binwidth);

       fclose(myrandfile);
       printf("Closing File.\n");

       #pragma simd
       for(i = 0; i < (N_rand); ++i){
          //MODEL EQ.
  	      //dist_index=(int)(floor((logD-logr_min)*((num_bins_preproc)/(logr_max-logr_min))));
          gridx_LIST[i]= (int) (floor((randX_LIST[i])*(grid_normalization)));
          gridy_LIST[i]= (int) (floor((randY_LIST[i])*(grid_normalization)));
          gridz_LIST[i]= (int) (floor((randZ_LIST[i])*(grid_normalization)));
          //gridx_LIST[i] = gxvar;
          //printf("GXVAR: %d", gxvar);

       }


       //Begin Nested For Loop to Count Pairs.
       //=========================
       // COMPUTE RAND DATA COUNTS
       //=========================
       printf("Beginning Random Nested Loops...\n");


       #pragma omp parallel shared( randZ_LIST, randY_LIST, randX_LIST, gridx_LIST, gridy_LIST, gridz_LIST, grid_numbins,\
       N_rand,chunk, num_threads_set) private (D, logD, actual_D,\
       logactual_D, x1, y1, z1, x2, y2, z2, delx, dely, delz, dist_index, i, j,grid_x1, grid_y1, grid_z1,\
       grid_x2, grid_y2, grid_z2, grid_x_check,grid_y_check,grid_z_check, grid_x_pbc,grid_y_pbc,grid_z_pbc, grid_x_norm, grid_y_norm,grid_z_norm)
       {
          omp_set_num_threads(num_threads_set);
          long int sum_local_counts[num_bins_preproc];
          memset(sum_local_counts, 0, num_bins_preproc * sizeof(sum_local_counts[0]) );


          #pragma omp for schedule(guided, chunk)
          for(i=0; i< (N_rand-1); ++i){

             x1 = randX_LIST[i];
             y1 = randY_LIST[i];
             z1 = randZ_LIST[i];

             grid_x1 = gridx_LIST[i];
             grid_y1 = gridy_LIST[i];
             grid_z1 = gridz_LIST[i];

             for(j=0; j< (N_rand-1); ++j){
                if (j!=i){
                   x2 = randX_LIST[j];
                   y2 = randY_LIST[j];
                   z2 = randZ_LIST[j];

                   grid_x2 = gridx_LIST[j];
                   grid_y2 = gridy_LIST[j];
                   grid_z2 = gridz_LIST[j];

                   /*
                   ------------------------------------
                   CHECK IF NECESSARY TO COMPUTE DISTANCE
                   ------------------------------------
                   */
                   //CHECKING CONDITIONS OF PERIOD BOX. Mathcing boxes and surrounding boxes 2 in each dimension
                   /*grid_x_check = ((grid_x1==grid_x2 || grid_x1==(grid_x2+1) || grid_x1==(grid_x2-1))||\
                                       (grid_x1==(grid_numbins-1) && grid_x2==0)||(grid_x2==(grid_numbins-1) && grid_x1==0) );

                   grid_y_check = ((grid_y1==grid_y2 || grid_y1==(grid_y2+1) || grid_y1==(grid_y2-1))||\
                                       (grid_y1==(grid_numbins-1) && grid_y2==0)||(grid_y2==(grid_numbins-1) && grid_y1==0) );
                   grid_z_check = ((grid_z1==grid_z2 || grid_z1==(grid_z2+1) || grid_z1==(grid_z2-1))||\
                                       (grid_z1==(grid_numbins-1) && grid_z2==0)||(grid_z2==(grid_numbins-1) && grid_z1==0) );
                  */

                   grid_x_norm = (grid_x1==grid_x2 || grid_x1==(grid_x2+1) || grid_x1==(grid_x2-1));  
                   grid_y_norm = (grid_y1==grid_y2 || grid_y1==(grid_y2+1) || grid_y1==(grid_y2-1)); 
                   grid_z_norm = (grid_z1==grid_z2 || grid_z1==(grid_z2+1) || grid_z1==(grid_z2-1));
                   grid_x_pbc = ((grid_x1==(grid_numbins-1) && grid_x2==0)|| (grid_x2==(grid_numbins-1) && grid_x1==0) );
                   grid_y_pbc = ((grid_y1==(grid_numbins-1) && grid_y2==0)|| (grid_y2==(grid_numbins-1) && grid_y1==0) );
                   grid_z_pbc = ((grid_z1==(grid_numbins-1) && grid_z2==0)|| (grid_z2==(grid_numbins-1) && grid_z1==0) );
                   
                   //CHECK TO VERIFY LOGIC CONDITIONS AND BINNING.
                   /*if (!(grid_x_norm && grid_y_norm && grid_z_norm)){
                      if (grid_x1==0 || grid_y1==0 || grid_z1==0){
                        if (grid_x2==(grid_numbins-1.) || grid_y2==(9) || grid_z2==(9)){
                         printf("grid values: %d %d %d \n", grid_x1, grid_y1, grid_z1);
                         printf("grid values: %d %d %d \n", grid_x2, grid_y2, grid_z2);
                         printf("----------------------------\n");
                        }
                      }
                    }*/
                   
                   grid_x_check = (grid_x_norm || grid_x_pbc);
                   grid_y_check = (grid_y_norm || grid_y_pbc);
                   grid_z_check = (grid_z_norm || grid_z_pbc);
                   //ONLY EXECUTE IF BOX IS APPROPIATE. 
                   if ((grid_x_check && grid_y_check && grid_z_check)){
                    //D = distance_given_2points(*x1, *y1, *z1, *x2, *y2, *z2);
                    //D = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
                    
                    /*if (grid_x_pbc || grid_y_pbc || grid_z_pbc){
                      printf("PBC CHECK IS TRUE!");
                    }*/

                     delx = fabs(x1 - x2);
                     dely = fabs(y1 - y2);
                     delz = fabs(z1 - z2);
                     //CORRECT FOR PERIOD BOUNDARY CONDITIONS
                     if((delx) > (2.*grid_binwidth)){
                       //printf("ENTERED PBC: X");
                       if (x1 > x2){
                         delx = ((x2+grid_binmax)-x1);
                       }
                       else {
                         delx = ((x1+grid_binmax)-x2);
                       }
                     }
                     if ((dely) > (2.*grid_binwidth)){
                       //printf("ENTERED PBC: Y");
                       if (y1 > y2){
                         dely = ((y2+grid_binmax)-y1);
                       }
                       else {
                         dely = ((y1+grid_binmax)-y2);
                       }

                     }
                    if ((delz) > (2.*grid_binwidth)){
                       //printf("ENTERED PBC: Z");
                       if (z1 > z2){
                         delz = ((z2+grid_binmax)-z1);
                       }
                       else {
                         delz = ((z1+grid_binmax)-z2);
                       }

                     }
                     //CALCULATE THE PSUEDO-DISTANCE
                     D = delx*delx + dely*dely + delz*delz;
                     logD = log10(D);
                     //actual_D = sqrt(D);
                     //logactual_D = log10(actual_D);

                     //CALCULATE THE PSUDO-DISTANCE INTO THE APPROPIATE BIN. 
                     if (logD < logr_max){
                       dist_index = (int) (floor((logD - logr_min)*((num_bins_preproc)/(logr_max - logr_min))));
                       //actual_Ddist_index = (int) (floor((logactual_D - logrNOTSQUARED_min)*((num_bins_preproc)/(logrNOTSQUARED_max - logrNOTSQUARED_min))));

                       //ADD COUNT TO APPROPIATE BIN.
                       if (dist_index >= 0 && dist_index < num_bins_preproc){
                          //Increment the appropiate bin.
                          //Non-vectorized: randdistance_counts[dist_index] += 1;
                          sum_local_counts[dist_index] += 1;
                       }
                     //end if statment for logr_max
                     }
                   //end grid check if statement
                   }
                }
             //end inner for loop
             }
          //end outer for loop
          }
          //CRITICAL SECTION TO SUM COUNTS ACROSS THREADS. 
          #pragma omp critical
          {
             //Sum up over the local counts on each thread to get total distance count.
             for(i=0 ; i < num_bins_preproc; ++i)
             {
                randdistance_counts[i] += sum_local_counts[i];
             }

          }

       //END PRAGMA PARALLEL BLOCK
       }

       printf("\n*");
       printf("\n   *");
       printf("\n     *");
       printf("\n       *");
       printf("\n      *");
       printf("\n     *");
       printf("\n    *");
       printf("\n   *        *");
       printf("\n   *");
       printf("\n     *");
       printf("\n      *");
       printf("\n       **");
       printf("\n        * *");
       printf("\n        * * *");
       printf("\n       * * * *\n");

       printf("************************************\n");
       printf("FINISHED DM RAND PRAGMA OMP CRITICAL\n");
       printf("************************************\n");
       printf("FINISHED RANDOM NESTED LOOPS! \n");
       printf("Counts: ");

       //Divides the Random counts by two and compute floor for integer.
       #pragma simd
       for(i =0; i< num_bins_preproc; ++i){
          randdistance_counts[i] = (long long) (floor(randdistance_counts[i]/2.)) ;
          //printf("%ld ", randdistance_counts[i]);

       }
    printf("\n") ;
     } //END IF NOT PERIODIC_BOUNDARY_CONDITIONS 


   // =========================================================
   /***
    *      ___           _     __  __      _   _
    *     |   \ __ _ _ _| |__ |  \/  |__ _| |_| |_ ___ _ _
    *     | |) / _` | '_| / / | |\/| / _` |  _|  _/ -_) '_|
    *     |___/\__,_|_| |_\_\ |_|  |_\__,_|\__|\__\___|_|
    */

   //static const char DM_datafile[] = "/home/chasonn/LargeScaleHW/HW3Codes/DMCorrelationFunctionParallel/DM.dat";
   //static const char DM_datafile[] = "/scratch/chasonnscr/positdirHPC/PositionXYZdatafile0.dat";
   char DM_datafile[100];
   int filenum = 0;
   int max_filenum = 1;

   /*
    MAIN LOOP TO ITERATE OVER MULTIPLE REDSHIFTS OR FILES AND COMPUTE Xi WITH SAME DM RANDOM.
   */
   for(filenum=0; filenum<max_filenum; filenum++){
     //sprintf(DM_datafile, "/scratch/chasonnscr/positdirHPC/PositionXYZdatafile%d.dat", filenum );
     //sprintf(DM_datafile, "/scratch/chasonnscr/HOSTSpositdirHPC/PositionXYZdatafileHOSTS%d.dat", filenum );

     sprintf(DM_datafile, "/scratch/chasonnscr/positdirHPC/PositionXYZdatafile%d.dat", filenum );
     //sprintf(DM_datafile, "/scratch/chasonnscr/known/DM_knownResearch.dat");


     FILE *myfile = fopen ( DM_datafile, "r" );
     if (myfile == NULL) {
        printf("DM.dat not opened, trying local filepath...DM.dat\n");
        static const char DM_datafile[] = "./DM.dat";
        myfile = fopen ( DM_datafile, "r" );
        if (myfile == NULL) {
          printf("DM.dat not opened, exiting...\n");
          exit(0);
        }
     }

     printf("Opened file - Begining assignment.\n");

     //long int N_data = 1000;
     long int N_data = FILELENGTHDM;
     double X_LIST[N_data];
     double Y_LIST[N_data];
     double Z_LIST[N_data];

     //==================
     // READ IN DM DATA.
     //==================
     //i = 0;
     //while(fscanf(myfile, "%f %f %f", &RA,&DEC,&Z) != EOF){
     //while(fscanf(myfile, "%f", &RA_LIST[i]) != EOF && i < N_data){
     //for(i = 0; i < (N_data); i++){
     i = 0;
     while( i<(N_data) && !feof(myfile)){
        fscanf(myfile, "%lf", &X_LIST[i]);
        fscanf(myfile, "%lf", &Y_LIST[i]);
        fscanf(myfile, "%lf", &Z_LIST[i]);
        //printf("i: %ld Data Z: %f \n", i, Z_LIST[i] );

        if (X_LIST[i] > (float) (box_side_length)) {
              box_side_length = (int) ceil(X_LIST[i]);
         }

        if (i >= (N_data-1)){
           printf("Close or exceeded N_data limit. X: %lf \n", X_LIST[i]);
           //break;

        }
        ++i;
     }


     fclose(myfile);
     printf("Closing File.\n");


     //CHECK STANDARD GRID_NUMBINS VALUE
     grid_binmax = box_side_length;
     grid_numbins = floor( grid_binmax / r_max );

     printf("Final Grid Bin Resolution: %d Bins......\n", grid_numbins);

     float grid_normalization = grid_numbins/(grid_binmax);
     float grid_binwidth = (grid_binmax/grid_numbins);
     printf("Final Grid Binwidth: %f width [mpc]\n", grid_binwidth);

     printf("correcting N_data size, due to different datafile. N_data: %ld i: %ld\n", N_data, i);
     N_data = i;

     #pragma simd
     for(i = 0; i < (N_rand); ++i){
         //MODEL EQ.
         //dist_index=(int)(floor((logD-logr_min)*((num_bins_preproc)/(logr_max-logr_min))));
         gridx_LIST[i]= (int) (floor(X_LIST[i]*(grid_normalization)));
         gridy_LIST[i]= (int) (floor(Y_LIST[i]*(grid_normalization)));
         gridz_LIST[i]= (int) (floor(Z_LIST[i]*(grid_normalization)));
      }

     //Begin Nested For Loop to Count Pairs. 

     //==========================
     // COMPUTE DM DATA COUNTS
     //==========================
     printf("Beginning Nested Loops...\n");
     printf("Note: using Distance Squared to avoid taking SQRT for index rank.\n");

     #pragma omp parallel shared( Z_LIST, Y_LIST, X_LIST, gridx_LIST, gridy_LIST, gridz_LIST, grid_numbins, N_data,\
     chunk, num_threads_set) private (D, logD, x1, y1, z1, x2, y2, z2, dist_index, i, j, tid, \
     grid_x1, grid_y1, grid_z1, grid_x2, grid_y2, grid_z2, grid_x_check, grid_y_check, grid_z_check)
     {
        //OMP_NUM_THREADS = 16; 
        //omp_set_num_threads(num_threads);
        omp_set_num_threads(num_threads_set);
        long int sum_local_counts[num_bins_preproc];
        memset(sum_local_counts, 0, num_bins_preproc * sizeof(sum_local_counts[0]) );
        tid = omp_get_thread_num();
        if (tid==0){
          nthreads = omp_get_num_threads();
          printf("thread id: %d, num threads: %d \n",tid, nthreads);
        }

        #pragma omp for schedule(guided, chunk)
        for(i=0; i < (N_data-1); ++i){
           x1 = X_LIST[i];
           y1 = Y_LIST[i];
           z1 = Z_LIST[i];

           grid_x1 = gridx_LIST[i];
           grid_y1 = gridy_LIST[i];
           grid_z1 = gridz_LIST[i];

           for(j=0; j < (N_data-1); ++j){
              if ( j!=i ){
                 x2 = X_LIST[j];
                 y2 = Y_LIST[j];
                 z2 = Z_LIST[j];

                 grid_x2 = gridx_LIST[j];
                 grid_y2 = gridy_LIST[j];
                 grid_z2 = gridz_LIST[j];

                /*if (grid_x1==15 && grid_x2==0){
                  printf("HALLEUIUEUA");
                }*/

                 grid_x_check = ((grid_x1==grid_x2 || grid_x1==(grid_x2+1) || grid_x1==(grid_x2-1)) ||\
                                     (grid_x1==(grid_numbins-1) && grid_x2==0)||(grid_x2==(grid_numbins-1) && grid_x1==0) );
                 grid_y_check = ((grid_y1==grid_y2 || grid_y1==(grid_y2+1) || grid_y1==(grid_y2-1)) ||\
                                     (grid_y1==(grid_numbins-1) && grid_y2==0)||(grid_y2==(grid_numbins-1) && grid_y1==0) );
                 grid_z_check = ((grid_z1==grid_z2 || grid_z1==(grid_z2+1)|| grid_z1==(grid_z2-1)) ||\
                                     (grid_z1==(grid_numbins-1) && grid_z2==0)||(grid_z2==(grid_numbins-1) && grid_z1==0) );
                 /*
                 grid_x_pbc = ((grid_x1==(grid_numbins-1) && grid_x2==0) || (grid_x2==(grid_numbins-1) && grid_x1==0));

                 grid_y_pbc = ((grid_y1==(grid_numbins-1) && grid_y2==0) || (grid_y2==(grid_numbins-1) && grid_y1==0));

                 grid_z_pbc = ((grid_z1==(grid_numbins-1) && grid_z2==0) || (grid_z2==(grid_numbins-1) && grid_z1==0));

                 if (grid_x_pbc||grid_y_pbc||grid_z_pbc) {
                    printf("Entered into some pbc parts with logic");
                 }*/

                 if (grid_x_check && grid_y_check && grid_z_check){
                    //D = distance_given_2points(*x1, *y1, *z1, *x2, *y2, *z2);
                    //Comment: actual D commented out, so not to comput sqrt just for index rank.
                     //D = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );


                   delx = fabs(x1 - x2);
                   dely = fabs(y1 - y2);
                   delz = fabs(z1 - z2);

                   //CORRECT FOR PERIOD BOUNDARY CONDITIONS
                  

                   if ((delx) > (2.*grid_binwidth)){
                     //printf("ENTERED THE PERIODIC BOUNDARY! X\n");
                     if (x1 > x2){
                       delx = ((x2+grid_binmax)-x1);
                     }
                     else {
                       delx = ((x1+grid_binmax)-x2);
                     }
                   }
                   if ((dely) > (2.*grid_binwidth)){
                     //printf("ENTERED THE PERIODIC BOUNDARY! Y\n");
                     if (y1 > y2){
                       dely = ((y2+grid_binmax)-y1);
                     }
                     else {
                       dely = ((y1+grid_binmax)-y2);
                     }

                   }
                  if ((delz) > (2.*grid_binwidth)){
                     //printf("ENTERED THE PERIODIC BOUNDARY! Z\n");
                     if (z1 > z2){
                       delz = ((z2+grid_binmax)-z1);
                     }
                     else {
                       delz = ((z1+grid_binmax)-z2);
                     }

                   }

                   D = delx*delx + dely*dely + delz*delz;
                   logD = log10(D);
                   //COMPUTE PSUEDO-DISTANCE INDEX.
                   dist_index = (int) (floor((logD - logr_min)*((num_bins_preproc)/(logr_max - logr_min))));
                   if (dist_index >= 0 && dist_index < num_bins_preproc){
                       //Increment the appropiate bin.
                       if (dist_index >= num_bins_preproc)
                          printf("YELLING!");
                       /*
                       //OLD NON-MULTITHREADED:
                       distance_counts[dist_index] += 1;
                       */
                       sum_local_counts[dist_index] += 1;
                   }
                 } //end grid check statement
              } //end i!=j if statement
           //end inner for loop
           }
        //end outer for loop
        }

        #pragma omp critical
        {
           //Sum up over the local counts on each thread to get total distance count.
           for(i=0 ; i < num_bins_preproc; ++i)
           {
              distance_counts[i] += sum_local_counts[i];
           }

        }

     //END PRAGMA PARALLEL BLOCK
     }
     
     printf("\n*");
     printf("\n   *");
     printf("\n     *");
     printf("\n       *");
     printf("\n      *");
     printf("\n     *");
     printf("\n    *");
     printf("\n   *        *");
     printf("\n   *");
     printf("\n     *");
     printf("\n      *");
     printf("\n       **");
     printf("\n        * *");
     printf("\n        * * *");
     printf("\n       * * * *\n");

     printf("****************************\n");
     printf("FINISHED DM PRAGMA OMP CRITICAL\n");
     printf("****************************\n");


     printf("FINISHED DM NESTED LOOP. \n");
     printf("Dividing Counts by two to correct double counting...");
     printf("Counts: ");
     #pragma simd
     for(i=0 ; i < num_bins_preproc; ++i)
     {
        distance_counts[i] = (long long) (floor(distance_counts[i]/2.)) ;
        printf("%ld ", distance_counts[i]);
     }
   
     printf("\n") ;

      if (periodic_boundary_conditions){
      //ANALYTIC EXPRESSION FOR PERIODIC BOX.   
      float pbc_rand_normalization = 0.5 * N_rand*N_rand/(box_side_length * box_side_length * box_side_length); //[n^2/V]
      float r1, r2; 
      printf("Random Counts Analytic:");
      for(i =0; i< num_bins_preproc; ++i){
          r1= sqrt(pow(10, ((logr_max - logr_min)/num_bins_preproc ) * i + logr_min));
          r2= sqrt(pow(10, ((logr_max - logr_min)/num_bins_preproc ) * (i+1) + logr_min));
          randdistance_counts[i] = (long long) (floor(pbc_rand_normalization*(4./3.)*PI*(pow(r2, 3) - pow(r1, 3)) ));
          printf("%ld ", randdistance_counts[i]);

       } 
    printf("\n") ;
    }//end periodic boundary conditions.
 




     //COUNTS COMPLETED.  {random, 21r, 20r, 20z}
     //==============================================================================
     /***
      *       ___ ___  _   _ _  _ _____ ___    ___ ___  __  __ ___ _    ___ _____ __ 
      *      / __/ _ \| | | | \| |_   _/ __|  / __/ _ \|  \/  | _ \ |  | __|_   _| __|
      *     | (_| (_) | |_| | .` | | | \__ \ | (_| (_) | |\/| |  _/ |__| _|  | | | _|
      *      \___\___/ \___/|_|\_| |_| |___/  \___\___/|_|  |_|_| |____|___| |_| |___|
      *      -------------------------------   ---------------------------------------
      */
     //==============================================================================

     //==============================================================================
     /***
      *       ___ ___  __  __ ___ _   _ _____ ___                                 
      *      / __/ _ \|  \/  | _ \ | | |_   _| __|                                
      *     | (_| (_) | |\/| |  _/ |_| | | | | _|                                 
      *      \___\___/|_|_ |_|_|__\___/  |_|_|___|___ ___  _  _                   
      *      / __/ _ \| _ \ _ \ __| |    /_\_   _|_ _/ _ \| \| |                  
      *     | (_| (_) |   /   / _|| |__ / _ \| |  | | (_) | .` |                  
      *      \___\___/|_|_\_|_\___|____/_/ \_\_|_|___\___/|_|\_|_   ___  __  __   
      *     | __| | | | \| |/ __|_   _|_ _/ _ \| \| |   | _ )/ _ \ / _ \|  \/  |  
      *     | _|| |_| | .` | (__  | |  | | (_) | .` |_  | _ \ (_) | (_) | |\/| |_ 
      *     |_|  \___/|_|\_|\___| |_| |___\___/|_|\_(_) |___/\___/ \___/|_|  |_(_)
      */
     //==============================================================================

     // =================================
     // COMPUTE THE CORRELATION FUNCTION
     // =================================

     //DM.
     printf("Calculating DM Correlation Function...\n");
     double ratio = (double) N_rand / (double) N_data;
     #pragma simd
     for(i=0; i<num_bins_preproc; i++){
        //Compute the Correlation Function: Xi = (Ngal/Nrand)^2 * (DD/RR  - 1)
        Xi_func[i] = ratio * ratio *  ( (double) distance_counts[i] / (double) randdistance_counts[i]) - 1.0 ;
        //CHECK VALUES:
        printf("%.2lf ", Xi_func[i]);
     }
     printf("\n");

     //----------------------------------------------------
     /***
      *      ___                 _         ___ ___ _    ___
      *     / __| __ ___ _____  | |_ ___  | __|_ _| |  | __|
      *     \__ \/ _` \ V / -_) |  _/ _ \ | _| | || |__| _|
      *     |___/\__,_|\_/\___|  \__\___/ |_| |___|____|___|
      */
     //----------------------------------------------------


     // ==========================
     // WRITING BINS RESULTS TO FILE
     // ==========================
     printf("Saving DM bins to file.\n");
     FILE *fp_out;
     fp_out = fopen("output_DM_allbins_Research.txt","w");
     if ( fp_out == NULL ) {
        printf("output_file.txt not opened, exiting...\n");
        exit(0);
     }
     for ( i=0 ; i <= num_bins_preproc ; i++ ) {
        bin_counts[i] = sqrt(pow(10, ((logr_max - logr_min)/num_bins_preproc ) * i + logr_min));
        fprintf(fp_out, "%lf \n", bin_counts[i]);
     }
     fclose(fp_out);


     /*
     // ==========================
     // WRITING X-GRID BINWIDTH RESULTS TO FILE
     // ==========================
     printf("Saving GRID BINWIDTH to file.\n");
     //FILE *fp_out;
     fp_out = fopen("output_DMGRIDBins.txt","w");
     if ( fp_out == NULL ) {
        printf("output_file.txt not opened, exiting...\n");
        exit(0);
     }
     for ( i=0 ; i < (N_data/2.) ; i++ ) {
        fprintf(fp_out,"%d %d %d\n", gridx_LIST[i], gridy_LIST[i], gridz_LIST[i]);
     }
     fclose(fp_out);
     */
  
     // ==========================
     // WRITING DM RESULTS TO FILE
     // ==========================
     printf("Saving DM counts to file.\n");
     //FILE *fp_out;
     fp_out = fopen("output_counts_DMdataResearchGRID.txt","w");
     if ( fp_out == NULL ) {
        printf("output_file.txt not opened, exiting...\n");
        exit(0);
     }
     for ( i=0 ; i < num_bins_preproc ; i++ ) {
        fprintf(fp_out,"%ld \n", distance_counts[i]);
     }
     fclose(fp_out);

     // ==========================
     // WRITING logr RESULTS TO FILE
     // ==========================
     printf("Saving logr counts to file.\n");
     fp_out = fopen("output_logrResearch.txt","w");
     if ( fp_out == NULL ) {
        printf("output_logr.txt not opened, exiting...\n");
        exit(0);
     }
     fprintf(fp_out,"%lf %lf %d \n", logrNOTSQUARED_min, logrNOTSQUARED_max, num_bins_preproc);
     fclose(fp_out);


     // ==========================
     // WRITING random RESULTS TO FILE
     // ==========================
     printf("Saving Random counts to file.\n");
     //FILE *fp_out;
     fp_out = fopen("output_counts_DMrandomResearchGRID.txt","w");
     if ( fp_out == NULL ) {
        printf("output_file.txt not opened, exiting...\n");
        exit(0);
     }
     for ( i=0 ; i < num_bins_preproc ; i++ ) {
        fprintf(fp_out,"%ld \n", randdistance_counts[i]);
	//fprintf(stderr, "%d", i);
     }
     fclose(fp_out);



     // ===============================
     // WRITING DM Xi RESULTS TO FILE
     // ===============================
     printf("Saving Xi-DM to file.\n");
     //FILE *fp_out;
     char Xifilename[50];
     sprintf(Xifilename, "output_Xi_DMResearchGRID_%d.txt", filenum);
     fp_out = fopen(Xifilename,"w");
     if ( fp_out == NULL ) {
        printf("output_file.txt not opened, exiting...\n");
        exit(0);
     }
     for ( i=0 ; i < num_bins_preproc ; i++ ) {
        fprintf(fp_out,"%f \n", Xi_func[i]);
     }
     fclose(fp_out);

   }
     printf("Done.\n");
     return 0;

  }



//END
