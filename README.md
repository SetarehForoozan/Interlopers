# Interlopers-in-BAO

If you store this is some folder called /STORING_ADRESS, then 

This project consist of these files:

1- default_config:
  Includes fitting parameters and cosmology parameters used when generating the correlation functions. You can change these parameters according to your correlation functions and the fitting parameters you want to use.
  
 
2- fitting_params.py:

  Takes some of the simulation and fitting parameters from you to use in the main code
  
3- mono_functions.py:
  
  Contains all functions needed for monopole only fit
  
 4- quad_functions.py:

  contains all functions needed for mono+quad fit
  
 5- fitting_mono.py:
 
  This is the main file to run. 
  Enter:
  # STORING_ADRESS/codes/new_fitting_codes/fit_mono.py $nsim $case $A $frac $dof_a $dof_q 0 $rsd
  
  nsim:
    The number of simulations
  case:
    0: no correction
    1: cross term measured directly from the simulations
    8: cross term is calculated from the convolution of the contaminated auto correlation. (Needs to be divided by a prefactor which is a function of fi,         that will be fitted in MCMC)
  A:
    Is the interloper displacement in Mpc/h
    
  frac:
  
    This is the true fraction of interlopers, used just for reading the catalogues that are contaminated by fraction "frac".
    
  dof_a, dof_q:
  
    Polynomials that need to be used in monopole, and quadrupole respectively. Usually chosen from 1,2,3.
    
  rsd:
    Whether or not measured in redshift space ot real space. for redshift space use "rsd", for real space use "norsd".
    
  The output of this code for example for rsd = "rsd" and A = "97", case = 8, 1000 sims, frac = 2, dega = 3, would be stored in this folder
  STORING_ADRESS/codes/rsd_97/mono_case8_1000sims_Percival
  
  There are three types of files that this code stores:
  1- A pdf file called frac2_dega3_sims0to999.pdf: Contains the visualizations. 
 
  2- AVG-STD/AVG_frac2_dega3_sims0to999.txt and STD_frac2_dega3_sims0to999.txt: These contain the average and standard deviations calculated from the best      fit MCMC code. The numbers stored are a1, a2, a3, b, alpha, fi, respectively.
 
  3- paper/FINAL_frac2_dega3_dega3_sims0to999.txt: The above numbers, but in scientific format to be copied and pasted in your Latex code. The last number     is reduced xi^2
  
  
  
  
  
  
  
    
  
  
