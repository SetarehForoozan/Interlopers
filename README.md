# Interlopers-in-BAO

Cloning this project would creat a folder called ./Interlopers, consisting of two examplary folders called "rsd_97" and "norsd_85" and a folder called "codes". If you measure your own correlation functions, you need to create another folder here and name it in the format of X_Y, where X is either "rsd" or "norsd", and Y is the displacement in Mpc/h. In which you must create two other folders: (In the following we assumed frac_2 and d = 97 Mpc/h)

### "rsd_97"
1- CovMats: 
  covmat_multi_t_frac_2_nsim1000_0to999.txt: the measured covariance matrix of monopole+quadrupole for fraction 2 and 1000 sims
  covmat_t_frac_2_nsim1000.txt: the measured covariance of monopole for fraction 2 and 1000 simulations

2- mean0to999 (sim 0 to sim 999):
  + poles_t_frac_2.txt: First column must be r, second, third, and forth columns are monopole, quadrupole and hexadecapoles of the measured correlation function.
  + xi_cross_L0_d97_frac2.txt: First column is r, second column is calculated monopole of the cross correlation from the convolution of the auto correlation, and third column is its error bar.
  + xi_cross_L2_d97_frac2.txt: Same as above, but for quadupole of the cross correlation.
3- Other folders will be created here as soon as you run the code. (Will be described later on)

  
### "codes"
This folder consist of the codes to fit BAO to the above-mentioned contaminated correlation functions:

1- default_config:

  Includes fitting parameters and cosmology parameters used when generating the correlation functions. You can change these parameters according to your correlation functions and the fitting parameters you want to use.
  
2- fitting_params.py:

  Takes some of the simulation and fitting parameters from you to use in the main code
  
3- mono_functions.py:
  
  Contains all functions needed for monopole only fit
  
4- quad_functions.py:

  contains all functions needed for mono+quad fit
  
5- make_covmat.py:

  There are two functions here just to adjust the fitting range of covariance matrices. The size of your saved multipole covariance matrix is [rmin to rmax + rmin to rmax]. But you need to change this to [rfitmin to rfitmax + rfitmin to rfitmax]. 
  

6- fit_mono.py: (main code for monolope only)
 
  This is the main file to run. 
  running this in your terminal:
  
  ```
  python3 ./Interlopers/codes/fit_mono.py $nsim $case $A $frac $dof_a $dof_q 0 $rsd
  ```
  
  nsim:
    The number of simulations (exapmle files are 1000)
  case:
    0: no correction
    1: cross term measured directly from the simulations (not active now)
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
  ./Interlopers/rsd_97/mono_case8_1000sims
  
  There are three types of files that this code stores in this location:
  + ./frac2_dega3_sims0to999.pdf: Contains the visualizations in a pdf file. 
 
  + ./AVG-STD/AVG_frac2_dega3_sims0to999.txt and STD_frac2_dega3_sims0to999.txt: These contain the average and standard deviations calculated from the best fit MCMC code. The numbers stored are a1, a2, a3, b, alpha, fi, respectively.
 
  + ./paper/FINAL_frac2_dega3_dega3_sims0to999.txt: The above numbers, but in scientific format to be copied and pasted in your Latex code. The last number is reduced xi^2
  
7- fit_multi.py: (main code for mono+quad)

  This is very similar to above, with the following changes:
  
  + Everything is stored in "./Interlopers/rsd_97/**multi**_case8_1000sims"
  + ./AVG-STD/AVG_frac2_dega3_sims0to999.txt and STD_frac2_dega3_sims0to999.txt: The numbers stored are a1, a2, a3, b, alpha, fi, **q1, q2, q3, epsilon** respectively.

  ```
  python3 ./Interlopers/codes/fit_multi.py $nsim $case $A $frac $dof_a $dof_q 0 $rsd
  ```

For questions you can email me: s2forooz@uwaterloo.ca
  
  
