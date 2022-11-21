# Comparison-of-tricubic-treecodes-regularized-stokeslet-kernel
Code for computing the regularized Stokeslet using three different tricubic treecodes: C-1 continuous, C-0 continuous and discontinuous tribubic treecodes

!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Authors:

  	Henry A. Boateng  (boateng@sfsu.edu) 
  	Department of Mathematics
  	San Francisco State University
  	San Francisco, CA
     
  	Svetlana Tlupova (tlupovs@farmingdale.edu)
  	Department of Mathematics
  	Farmingdale State College, SUNY
  	Farmingdale, NY
  
  This material is based upon work partially supported by the 
  National Science Foundation under Grant Nos. CHE-1800181 and DMS-2012371
  and U.S. Department of Energy, Office of Science, Office of Work- force 
  Development for Teachers and Scientists (WDTS) under the Visiting Faculty Program.
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   NOTE: Please include the following references in any work that
         utilizes this code:
         
         (1) Boateng. H. A., Tlupova, S.: The effect of smoothness on
         the accuracy of treecode approximations
            Communications in Computer Physics, (2022)  
		 
        (2) Boateng. H. A., Tlupova, S.: A treecode algorithm based on 
            tricubic interpolation
            Journal of Computational Mathematics and Data Science, (2022)  
	    
Summary of files :
------------------

      rand*.txt       : Files containing the coordinates of the particles
                        distributed in [-5.0,5.0]x[-5.0,5.0]x[0.0,10.0].
                        The files rand_10000.txt, rand_80000.txt and rand_640000.txt
                        correspond to 10000, 80000, and 640000 particles respectively.
                        
      lambda*.txt      : Files containing the vector weights of the particles in the
                         corresponding rand*.txt
                        
                      

!!!!!!!
! .cpp and .h FILES 
      direct_sum.cpp : C++ file for exact computation of the potential and field. 
      
      Tricubic_RS_C1.cpp       : Main C++ program for the global-C1 tricubic treecode.
      			                     This version uses analytical derivatives of the kernel.
                   
      Tricubic_RS_C0.cpp       : Main program for the global-C0 tricubic treecode
                                  This version uses analytical derivatives of the kernel.
			
      Tricubic_RS_DC.cpp       : Main program for the discontinuous tricubic treecode. 
                                 This version uses analytical derivatives of the kernel.
                             
      Tricubic_RS_DC_fdiff.cpp : Main program for the discontinuous tricubic treecode which uses
                                 finite differences for the derivatives of the kernel
 
			  
                       The main programs direct_sum.cpp, Tricubic_RS_C1.cpp, Tricubic_RS_C0.cpp,
		                   Tricubic_RS_DC.cpp and Tricubic_RS_DC_fdiff.cpp depend on the following helper files 
		                   (utilities.h, utilities.cpp, kernel_utils.h,kernel_utils.cpp, tree_utils.h, 
                       tree_utils.cpp, tricubic_utils.h, tricubic_utils.cpp)
      
      makefile       : A makefile for compiling the direct sum and tricubic treecodes. 
                       It produces the executables direct_sum and Tricubic_RS_C1,
		                    Tricubic_RS_C0, Tricubic_RS_DC, and Tricubic_RS_DC_fdiff
                        
      input_direct.txt   : Input file for the direct_sum
      input_params.txt   : Input file for the tricubic treecodes

Input for direct_sum  :
-----------------------

     input_direct.txt specifies the following two required options on the same line
     
     number of particles, regularization parameter for regularized stokeslet

Input for tricubics :
-----------------------------------

      input_params.txt specifies the following required options in
      the given order:
      
#// Treecode method to be used: Particle-Cluster

	Particle-Cluster    #// (This option picks particle-cluster)
 
// number of particles 

	10000   // (If the system has 10000 particles)
 
// N0, leaf size 

	1000     // Pick N0 = 1000
 
// parameter to choose the type of MAC (0 - Regular MAC, 1 - Spherical MAC)

	0   // Use Regular MAC
 
// theta, MAC parameter (Should be less than 1 for Regular MAC)

	0.5  //theta = 0.5
  
// epsilon, regularization parameter for regularized stokeslet

0.02 
 
// 1 if finite differences should be used, 0 otherwise (Should be 1 for Tricubic_RS_DC_fdiff and 0 for the other three)

	0  // Use analytical derivatives
 
// mgrid, any value when finite differences are not used. Do not delete.

	20  # Grid size for 3-point interpolation of the kernel. Can be varied for varying accuracy
 
// hgrid, any value when finite differences are not used. Do not delete.

	0.01 # step-size for finite differences (can be varied as well)

Output for the executable direct_sum :
-------------------------------------

The output is the file exact_sum_RegStokeslet_Nnumberofparticles. 
Example: The output for  10000 particles is exact_sum_RegStokeslet_N10000

For a system of size N, the file has N+1 lines. The first N lines correspond to data for
the N particles. Each of the first N lines has 5 entries:

	Particle number,  x-component of velocity, y-component of velocity, z-component of velocity
   
The N+1 line is the time the computation took in 10^(-2) seconds (i.e. centiseconds)


Output for the executables Tricubic_RS_C1, Tricubic_RS_C0, Tricubic_RS_DC  :
-------------------------------------------------------------------

Running Tricubic_RS_C1 generates two files : output_RS_C1.txt and tricubic_sum_RegStokeslet_RS_C1_Nnumberofparticles
Running Tricubic_RS_C0 generates two files : output_RS_C0.txt and tricubic_sum_RegStokeslet_RS_C0_Nnumberofparticles
Running Tricubic_RS_DC generates two files : output_RS_DC.txt and tricubic_sum_RegStokeslet_RS_DC_Nnumberofparticles
Running Tricubic_RS_DC_fdiff generates two files : output_RS_DC_fdiff.txt and tricubic_sum_RegStokeslet_RS_DC_fdiff_Nnumberofparticles

The file output_RS_*.txt has the following entries:

	Number of particles (N_cube), Maximum number of particles in a leaf (N0), MAC (theta), relative 2-norm error in velocity, divergence, per-particle divergence,  time for direct sum in centiseconds, treecode time in centiseconds

For the files tricubic_sum_RegStokeslet_RS_C1_Nnumberofparticles, tricubic_sum_RegStokeslet_RS_C0_Nnumberofparticles,  
tricubic_sum_RegStokeslet_RS_DC_Nnumberofparticles, or tricubic_sum_RegStokeslet_RS_DC_fdiff_Nnumberofparticles:

    The first line has the data:
    cpu time for direct sum (in centiseconds), cpu time for treecode (in centiseconds)
    
    Following the first line are direct sum and the treecode approximation for the velocity and per-particle divergence in the sample format:
    
    particle number,  exact x-component of velocity, exact y-component of velocity, exact z-component of velocity, 0 
    particle number,  tree approximation of x-component of velocity, tree approximation of y-component of velocity, tree approximation of z-component of velocity, tree approximation of per-particle divergence
 


		 
