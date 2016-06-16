#!/bin/bash

#Script to read in snapshots and plot
#Initially plot single snapshot. Should be extending to plot multiple snaps.

echo "Hello, this program is working fine so far."
echo " Make sure the paths are correct for loading in files."
echo "This program uses python and unix commands."

#FILE NAME TYPE: tree_0_0_0.dat   0-10,0-10,0-10


#SET THESE VARIABLES TO the start and end snap.
#They will read sequentally unless programmed otherwise.

startsnap1=1
startsnap2=0
startsnap3=0

endsnap1=7
endsnap2=10
endsnap3=10

increment1=1 #not currently set
increment2=1 #not currently set
increment3=1 #not currently set


##scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) 
# upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) 
# rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) 
# x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) 
# Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29)
# Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID (32) 
# Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) 
# Tidal_Force(35) Tidal_ID(36) 
# Rs_Klypin M200b_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a 
# A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c)
# A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius

#Omega_M = 0.250000; Omega_L = 0.750000; h0 = 0.700000
#Full box size = 130.000000 Mpc/h

#PARAMETERS
#Scale: Scale factor of halo.
#ID: ID of halo (unique across entire simulation).
#Desc_Scale: Scale of descendant halo, if applicable.
#Descid: ID of descendant halo, if applicable.
#Num_prog: Number of progenitors.
#Pid: ID of least massive host halo (-1 if distinct halo).
#Upid: ID of most massive host halo (different from Pid when the halo is within two or more larger halos).
#Desc_pid: Pid of descendant halo (if applicable).
#Phantom: Nonzero for halos interpolated across timesteps.
#SAM_Mvir: Halo mass, smoothed across accretion history; always greater than sum of halo masses of contributing progenito
#rs (Msun/h).  Only for use with select semi-analytic models.
#Mvir: Halo mass (Msun/h).
#Rvir: Halo radius (kpc/h comoving).
#Rs: Scale radius (kpc/h comoving).
#Vrms: Velocity dispersion (km/s physical).
#mmp?: whether the halo is the most massive progenitor or not.
#scale_of_last_MM: scale factor of the last major merger (Mass ratio > 0.3).
#Vmax: Maxmimum circular velocity (km/s physical).
#X/Y/Z: Halo position (Mpc/h comoving).
#VX/VY/VZ: Halo velocity (km/s physical).
#JX/JY/JZ: Halo angular momenta ((Msun/h) * (Mpc/h) * km/s (physical)).
#Spin: Halo spin parameter.
#Breadth_first_ID: breadth-first ordering of halos within a tree.
#Depth_first_ID: depth-first ordering of halos within a tree.
#Tree_root_ID: ID of the halo at the last timestep in the tree.
#Orig_halo_ID: Original halo ID from halo finder.
#Snap_num: Snapshot number from which halo originated.
#Next_coprogenitor_depthfirst_ID: Depthfirst ID of next coprogenitor.
#Last_progenitor_depthfirst_ID: Depthfirst ID of last progenitor.
#Last_mainleaf_depthfirst_ID: Depthfirst ID of last progenitor on main progenitor branch.
#Tidal_Force: Strongest tidal force from any nearby halo, in dimensionless units (Rhalo / Rhill).
#Tidal_ID: ID of halo exerting strongest tidal force.
#Rs_Klypin: Scale radius determined using Vmax and Mvir (see Rockstar paper)
#M200b_all: Mass enclosed within the specified overdensity, including unbound particles (Msun/h)
#M200b--M2500c: Mass enclosed within specified overdensities (Msun/h)
#Xoff: Offset of density peak from average particle position (kpc/h comoving)
#Voff: Offset of density peak from average particle velocity (km/s physical)
#Spin_Bullock: Bullock spin parameter (J/(sqrt(2)*GMVR))
#b_to_a, c_to_a: Ratio of second and third largest shape ellipsoid axes (B and C) to largest shape ellipsoid axis (A) (di
#mensionless).
#  Shapes are determined by the method in Allgood et al. (2006).
#  (500c) indicates that only particles within R500c are considered.
#A[x],A[y],A[z]: Largest shape ellipsoid axis (kpc/h comoving)
#T/|U|: ratio of kinetic to potential energies
#M_pe_*: Pseudo-evolution corrected masses (very experimental)




##scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) 
# upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) 
# rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) 
# x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) 
# Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29)
# Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID (32) 
# Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) 
# Tidal_Force(35) Tidal_ID(36) 
# Rs_Klypin M200b_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a 
# A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c)
# A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius


#SET THIS VARIABLE TO THE NUMBER OF PARTICLES REQUIRED TO CONSIDER HALO.
#Note: Anything lower than 20 is redundant since numparticles > 20 for all halos.
#Need to import to awk function.
halo_particle_min=200

#fPATH="/zpool0/fs2/shared/new_130_mpc_box/hires/5045/trees/so_m200b/"
#fFILE="tree_0_0_0.dat"

var1=0
var2=0
var3=0

#Loop is not well definied yet for proper looping over .list files
for ((var1=$startsnap1; var1 <= $endsnap1; var1=var1+$increment))
  do

	for ((var2=$startsnap2; var2 <= $endsnap2; var2=var2+$increment2))
	  do

		for ((var3=$startsnap3; var3 <= $endsnap3; var3=var3+$increment3))
  		  do
    	    echo "The loop is working. Current snap: $var1 $var2 $var3"
            echo "Trying to open file: "
		 	echo "Creating Data file ${var}"
			
			#1.
			#Open the tree files and remove all comment lines denoted with '#'
			#cat /zpool0/fs2/shared/new_130_mpc_box/hires/5045/rockstar_halos/so_m200b/trees/tree_${var1}_${var2}_${var3}.dat \
			#| grep -v "#" | awk '($2>=0) && ($18>=0) && ($19>=0) && ($20>=0) {print $0 " " $2 " " $10 " " $11}' > datafile${var}.txt 
			#awk '(NR>= 50) {print $1 " " $2 " " $10}' /zpool0/fs2/shared/new_130_mpc_box/hires/5045/rockstar_halos/so_m200b/trees/tree_${var1}_${var2}_${var3}.dat  > datafile${var}.txt 
			gunzip -c /zpool0/fs2/shared/new_130_mpc_box/hires/5045/rockstar_halos/so_m200b/full_res/hlists/hlist_${var1}.dat.gz | \
			   awk '(NR>=61 && $11 >= 2750000000) {print $1 " " $2 " " $11 " " $12 " " $18 " " $19 " " $20 " " $17 " " $6 " " $27 }' \
			    > /scratch/chasonnscr/positdir/MAdatafile${var1}.txt 

			#2=ID; 10=Mvir; 11=Rvir; 18,19,20=X,Y,Z 
			#cat initdata${var}.txt | awk '($2>=0) && ($18>=0) && ($19>=0) && ($20>=0) {print $2 " " $10 " " $11 " " $18 " " $19 " " $20}' > datafile${var}.txt 
		        #cat datafile9$var.txt | awk '{print $2 " " $3}' > pyplotdata9$var.txt
			#let "var1++"
			#let "var2++"
			#let "var3++"
		  done
	  done
  done 	  

echo "The script has finished running."
echo
echo "In order to plot the data points, run haloplotter.py"
echo "Make certain to check the PATH"
echo "The results should print out to .png files"


