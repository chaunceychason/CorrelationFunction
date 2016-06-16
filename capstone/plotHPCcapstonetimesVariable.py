import sys
import numpy as np
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter

#THIS PROGRAM PLOTS DATA IN FORM 'x, y'

#filename = 'time_stats_run1.txt'

#From Serial with ICC
#TOTAL FLOPS: 1 processor W/ error. = 191243625051
#TOTAL FLOPS: 1 processor double counting = 106246175034*2
#TOTAL FLOPS: 1 processor 106246175034
#real	2m7.472s
#user	2m7.480s
#sys	0m0.011s

#From Serial with GCC:
#TOTAL FLOPS: 106246175034
#real	8m7.327s
#user	8m7.391s
#sys	0m0.012s



ParaFLOPScount = 106246175034*2
FLOPScount = 106246175034
#Efficent Serial Time
serial_time = 2. * 60. + 7.472  #seconds
slow_serial_time = 253.73

filename1 = 'time_stats_TPCF_omp1.txt'
filename2 = 'time_stats_TPCF_omp4.txt'
filename3 = 'time_stats_TPCF_omp8.txt'
filename4 = 'time_stats_TPCF_omp16.txt'
filename5 = 'time_stats_TPCF_omp32.txt'
filename6 = 'time_stats_TPCF_omp64.txt'

filenameXi   = 'output_Xi_DM.txt'
filenameXiserial   = 'output_Xi_DM_serial.txt'
filenamelogr = 'output_logr.txt'
#Xi = np.loadtxt(filenameXi, usecols=(0,), unpack=True)
Xiserial = np.loadtxt(filenameXiserial, usecols=(0,), unpack=True)


logr_min, logr_max, num_bins = np.loadtxt(filenamelogr, usecols=(0, 1, 2), unpack=True)
print logr_min
print logr_max
#r_array = np.logspace(logr_min, logr_max, num=num_bins)
r_array = np.linspace(logr_min, logr_max, num=num_bins)


#Threads:
process_list = np.array([ 1, 4, 8, 16, 32, 64])  #Note idx[0] is 1.


time_offset_secs = 0.1 #This is a standard time difference between recorded time
						#and the actual time reported by process 0. It varied only
						#slightly amoung different runs, and thus should not skew.

npdata1 = np.loadtxt(filename1, usecols=(0,), unpack=True) - time_offset_secs
npdata2 = np.loadtxt(filename2, usecols=(0,), unpack=True) - time_offset_secs
npdata3 = np.loadtxt(filename3, usecols=(0,), unpack=True) - time_offset_secs
npdata4 = np.loadtxt(filename4, usecols=(0,), unpack=True) - time_offset_secs
npdata5 = np.loadtxt(filename5, usecols=(0,), unpack=True) - time_offset_secs
npdata6 = np.loadtxt(filename6, usecols=(0,), unpack=True) - time_offset_secs


y_values = np.array([npdata1,npdata2,npdata3,npdata4,npdata5, npdata6])
print y_values

print("Data should be structured (Real, User, Sys)... [-1] sets Real.")


"""
=====================================
        Xi  PLOT 
=====================================
"""
plot_title="Correlation Func. (Xi) vs Radius"
x_axis="logR"
y_axis="Xi(r)"

plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

max_num_redshifts = 8
colors = cm.rainbow(np.linspace(0,1, max_num_redshifts))
#for idx in range(0, max_num_redshifts):
for idx, c in zip(range(0, max_num_redshifts), colors):
	filenameXi = 'output_Xi_DM_%d.txt' % idx
	Xi = np.loadtxt(filenameXi, usecols=(0,), unpack=True)
	plt.scatter(r_array, Xi  , marker='o', color=c, label=idx)
	plt.plot(r_array, Xi, color=c )


plt.scatter(r_array, Xiserial  , marker='x', s=50, color='b', label='DM_Full')
plt.yscale("log")
#plt.xscale("log")
#plt.xlim([1, 64])
#---------------	
#plt.ylim(0, 1400)	
plt.legend(loc='best')

#Saves Plot
tmp_filename = "HPC_capstone_plotXi.png"
plt.savefig(tmp_filename)
plt.clf()


"""
=====================================
         REAL TIME  PLOT 
=====================================
"""
plot_title="Real Time vs. OpenMP Threads"
x_axis="Threads"
y_axis="Real Time [sec]"

plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.scatter(process_list, y_values  , marker='o',  color='r' )

plt.plot(process_list, y_values  , c='r')
plt.xlim([1, 64])
#---------------	
#plt.ylim(0, 1400)	
#plt.legend(loc='best')

#Saves Plot
tmp_filename = "HPC_capstone_plotREAL.png"
plt.savefig(tmp_filename)
plt.clf()

"""
=====================================
        LOG REAL TIME  PLOT 
=====================================
"""
plot_title="LOG Real Time vs. OpenMP Threads"
x_axis="Threads"
y_axis="LOG Real Time [sec]"

plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.scatter(process_list, y_values  , marker='o',  color='r' )

plt.plot(process_list, y_values  , c='r')
plt.xlim([1, 64])
#---------------	
#plt.ylim(0, 1400)	
#plt.legend(loc='best')
plt.yscale('log')
#Saves Plot
tmp_filename = "HPC_capstone_plotREALlog.png"
plt.savefig(tmp_filename)
plt.clf()



"""
=====================================
         SPEEDUP  PLOT 
=====================================
"""
plot_title="SPEEDUP vs. OpenMP"
x_axis="Threads"
y_axis="Speedup"

#Times for SERIAL particles
#serial_time = 1

plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.scatter(process_list, (serial_time/y_values ) , marker='o',  color='g',  label='1 Node' )

plt.plot(process_list, (serial_time/y_values ) , c='g')

plt.xlim([1, 64])
#---------------	
#plt.ylim(0, 1400)	
#plt.legend(loc='best')

#Saves Plot
tmp_filename = "HPC_capstone_plotREALspeedup.png"
plt.savefig(tmp_filename)
plt.clf()


"""
=====================================
         Efficiency  PLOT 
=====================================
"""
plot_title="Efficency vs. OpenMP Threads"
x_axis="Threads"
y_axis="Efficiency"

plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.scatter(process_list, (serial_time/y_values/ process_list ) , marker='o',  color='r',  label='OpenMP' )

plt.plot(process_list, (serial_time/y_values / process_list) , c='b')
plt.xlim([1, 64])
#---------------	
#plt.ylim(0, 1400)	
plt.legend(loc='best')

#Saves Plot
tmp_filename = "HPC_capstone_plotREALefficiency.png"
plt.savefig(tmp_filename)
plt.clf()


"""
=====================================
         FLOPS  PLOT 
=====================================
"""
plot_title="FLOPS vs. OpenMP Threads"
x_axis="Threads"
y_axis="FLOPS"


plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

flopsarray = ParaFLOPScount / y_values / process_list 
effective_flopsarray = FLOPScount / y_values / process_list 
serial_flops = FLOPScount / serial_time 

plt.scatter(process_list, flopsarray , marker='o',  color='r',  label='openMP' )
plt.scatter(process_list, effective_flopsarray , marker='o',  color='b',  label='eff. openMP' )
plt.scatter(1.05, serial_flops , marker='D',  color='g',  label='serial' )

plt.plot(process_list, flopsarray , c='r')
plt.plot(process_list, effective_flopsarray , c='b')
plt.xlim([1, 64])
#---------------	
#plt.ylim(0, 1400)	
plt.legend(loc='best')

#Saves Plot
tmp_filename = "HPC_capstone_plotREAL_FLOPS.png"
plt.savefig(tmp_filename)
plt.clf()





print "The program has finished running. All files closed. \nThe results should be in your directory"



#End
