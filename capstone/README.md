#Capstone Project. Parallel Two-Point

This Directory has lots of files. This README will make it easy to find what you want. 

If you JUST got here, ignore many of the files. 

COMPILE:
First, you probably need to build the code: FIRST: make clean SECOND: make 
	(This will make three version of the code, a serial, a parallel, and a variable parallel)

RUNNING C CODE:
If you want to run a single iteration of the Parallel Two Point CF type: ./DM_CF_Parallel 32

If you want to run a single iteration with different threads use the .sh script included! type: ./run_TPCF.sh #

If you want to run multiple iteratins to go through the Parallel Benchmarks, type ./RUN_ALL_STATS.sh


PLOTTING:
Once the .c code has finished you can make plots using python plotplotHPCcapstonetimes.py. 


Additional Info:

The files to run the variable version are located in /scratch/chasonn/ and may not be accessible. 
Only DM.day and DM_random.dat are needed for the code to work. 

The output files generated include the raw counts for many files which may get overwritten. The other
output files of Xi are listed by either just pure DM or with an index showing the order generated which
corresponds to the redshift index value. 



TIMES:

-rw-r--r-- 1 chasonn lss 18 Apr 21 22:30 time_stats_TPCF_omp100.txt

-rw-r--r-- 1 chasonn lss 18 Apr 21 02:05 time_stats_TPCF_omp16.txt

-rw-r--r-- 1 chasonn lss 19 Apr 21 02:11 time_stats_TPCF_omp1.txt

-rw-r--r-- 1 chasonn lss 19 Apr 21 22:31 time_stats_TPCF_omp200.txt

-rw-r--r-- 1 chasonn lss 18 Apr 21 23:03 time_stats_TPCF_omp32.txt

-rw-r--r-- 1 chasonn lss 18 Apr 21 02:05 time_stats_TPCF_omp4.txt

-rw-r--r-- 1 chasonn lss 18 Apr 21 23:04 time_stats_TPCF_omp64.txt

-rw-r--r-- 1 chasonn lss 18 Apr 21 22:30 time_stats_TPCF_omp75.txt

-rw-r--r-- 1 chasonn lss 18 Apr 21 02:05 time_stats_TPCF_omp8.txt

[chasonn@bender capstone]$ cat time_stats_TPCF_omp*

14.04 424.16 4.07

18.90 297.90 0.14

253.73 253.75 0.01

14.44 426.68 10.70

14.91 453.43 0.63

64.12 256.04 0.02

14.04 421.59 2.45

14.24 423.79 3.57

33.11 263.68 0.07


