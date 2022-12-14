# Correcting for small-displacement interlopers in BAO analyses
You can use this code for modelling interlopers in BAO analysis.
This code was used for generating the results of our paper https://arxiv.org/abs/2208.05001.

# Installation
To clone this project type the following command in your command-line:

```
git clone https://github.com/SetarehForoozan/Interlopers.git
```

This command would create a folder named _./Interlopers_, consisting of two examplary folders: _rsd_97_ and _norsd_85_. And another folder containing the codes named _codes_. 
If you measure your own correlation functions, you need to create another folder in _./Interlopers_ and name it in the format of X_Y, where X is either "rsd" or "norsd", and Y is the displacement in Mpc/h. In this folder, you must create two other folders such as ones in _rsd_97_ below: (In the following we assumed frac = 2 and d = 97 Mpc/h)

## _rsd_97_
1- CovMats: 
  
  + _covmat_multi_t_frac_2_nsim1000_0to999.txt_: the measured covariance matrix of monopole + quadrupole for fraction 2 and 1000 simulations
  
  + _covmat_t_frac_2_nsim1000.txt_: the measured covariance of monopole for fraction 2 and 1000 simulations

2- _mean0to999_ (sim 0 to sim 999):
  + _poles_t_frac_2.txt_: First column must be r, second, third, and forth columns are monopole, quadrupole and hexadecapoles of the measured correlation function.
  + _xi_cross_L0_d97_frac2.txt_: First column is r, second column is calculated monopole of the cross correlation from the convolution of the auto correlation, and third column is its error bar.
  + _xi_cross_L2_d97_frac2.txt_: Same as above, but for quadupole of the cross correlation.
  + 
3- Other folders will be created here as soon as you run the code. (Will be described later on)

  
## _codes_
This folder consist of the codes to fit BAO to the above-mentioned contaminated correlation functions:

### 1- _default_config.py_:

  Includes fitting parameters and cosmology parameters used when generating the correlation functions. You can change these parameters according to your correlation functions and the fitting parameters you want to use.
  
### 2- _fitting_params.py_:

  Takes some of the simulation and fitting parameters from you to use in the main code
  
### 3- _mono_functions.py_:
  
  Contains all functions needed for monopole only fit
  
### 4- _quad_functions.py_:

  contains all functions needed for mono+quad fit
  
### 5- _make_covmat.py_:

  There are two functions here just to adjust the fitting range of covariance matrices. The size of your saved multipole covariance matrix is [rmin to rmax + rmin to rmax]. But you need to change this to [rfitmin to rfitmax + rfitmin to rfitmax]. 
  

### 6- _fit_mono.py_: (main code for monolope only)
  #### Description 
 
  This is the main file to run. 
  You can give this command in your terminal to run the code:
  
  ```
  python3 ./Interlopers/codes/fit_mono.py $nsim $case $A $frac $dof_a $dof_q 0 $rsd
  ```
  
  nsim:
  + The number of simulations (in the exapmle files is 1000)
  
  case:
  + 0: no correction
  + 1: cross term measured directly from the simulations (not active now)
  + 8: cross term is calculated from the convolution of the contaminated auto correlation. (Needs to be divided by a prefactor which is a function of fi, that will be fitted in MCMC)
  
  A:
  + Is the interloper displacement in Mpc/h
    
  frac:
  + This is the true fraction of interlopers, used just for reading the catalogues that are contaminated by fraction "frac".
    
  dof_a, dof_q:
  + Polynomials that need to be used in monopole, and quadrupole respectively. Usually chosen from 1,2,3.
    
  rsd:
  + Whether or not measured in redshift space or real space. For redshift space use "rsd", for real space use "norsd".
    
    
### 7- _fit_multi.py_: (main code for mono+quad)

  This is very similar to number 6, with the following changes:
  + Everything is stored in _./Interlopers/rsd_97/**multi**_case8_1000sims_
  + _./AVG-STD/AVG_frac2_dega3_sims0to999.txt_ and _STD_frac2_dega3_sims0to999.txt_: The numbers stored are a1, a2, a3, b, alpha, fi, **q1, q2, q3, epsilon** respectively.

  ```
  python3 ./Interlopers/codes/fit_multi.py $nsim $case $A $frac $dof_a $dof_q 0 $rsd
  ```
  
  # Test
  
  If you run this code for nsim = 1000, case = 8, A = "97", frac = 2, dof_a = 3, dof_q = 3, and rsd = "rsd",
  
  ```
  python3 ./Interlopers/codes/fit_mono.py 1000 8 97 2 3 3 0 "rsd"
  ```
  
  you should see this folder appear: ./Interlopers/rsd_97/mono_case8_1000sims, which contains three types of files:
  
  + _./frac2_dega3_sims0to999.pdf_: Contains the visualizations in a pdf file. 
  + _./AVG-STD/AVG_frac2_dega3_sims0to999.txt_ and STD_frac2_dega3_sims0to999.txt: These contain the average and standard deviations calculated from the best fit MCMC code. The numbers stored are a1, a2, a3, b, alpha, fi, respectively.
  + _./paper/FINAL_frac2_dega3_dega3_sims0to999.txt_: The above numbers, but in scientific format to be copied and pasted in your Latex code. The last number is the reduced xi^2.
  


For questions you can email me at s2forooz@uwaterloo.ca
  
  
