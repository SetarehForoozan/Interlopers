
r"""
  This file contains some default parameters:
   - Fitting parameters
   - cosmology parameters used when generating the correlation functions. 
"""

import numpy as np

z_med = 1

######## Fitting range ########
rfitmin = 50
rfitmax = 150
sig_s_fog = 3.0

####### Measured correlation function seperations ########
rstep = (150-0.5)/101
rmin = 0.5
rmax = 150
r_old = np.arange(rmin,rmax,rstep)
r_list = (r_old[1:]+r_old[:-1])/2.0 # central
dr = r_old[1]-r_old[0]

###### MCMC walker information #######
nwalkers = 30
nsteps_mono = 4000
ndiscard_mono = 2000
nsteps_multi = 5000
ndiscard_multi = 2500

##### Cosmology ########
ns = 0.9624
omega_m = 0.3175
omega_b = 0.049
omega_c = omega_m - omega_b

z_med_p = '1'
#Add the redshift information here if it's not 1
if z_med == 1:
  bias = 2.7
  sigma8_lin =  0.834 * 0.6059
  signl_ap = 4.771
  z_med_p = '1'
  
snl = signl_ap

