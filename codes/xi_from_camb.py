import numpy as np
from spl import *
from readcamb import *
from apodize import *
from p2xi import *
from default_config import *
from stream import *
from read_data import *
from fitting_params import *

def xi_camb_func(ps_model):
    print("mono rsd:", rsd)
    #setare: Transfer function with wigglez from Graham.
    knw, pnw = np.loadtxt("./Interlopers/codes/fiducial/CAMB_nw_matterpow_%s.dat"%(z_med), unpack=True)
    fPlin   = "./Interlopers/codes/fiducial/%s_matterpow_%s.dat"%(ps_model,z_med)
    kr, J0 = np.loadtxt('./Interlopers/codes/fiducial/J0filt.dat', unpack=True)
    ksp, psp = np.loadtxt(fPlin, unpack=True)
    tspnw = spllog(ksp, pnw, knw) # where knw and pnw are the k and P from the no-wiggle file
    pm = apodize_Pk(tspnw, psp, ksp, snl) #power spectrum with damped bao #setare: changed pnw to tspnw
    pm *= bias**2
    #UNCOMMENT if want model with RSD
    if rsd == 'rsd':
        omega_mz = omega_m*(1.+z_med)**3/(omega_m*(1.+z_med)**3+1.-omega_m)
        nk = np.size(ksp)
        beta = omega_mz**0.55/bias
        sig_s = sig_s_fog #fit for this parameter?
        rsd_factor = calc_stream(ksp, nk, beta, sig_s, 'exp')
        pm = pm*rsd_factor

    xi_mod = p2xi_lin(ksp, pm, kr, J0, r_list)
    return xi_mod

def xi_camb_NL_func(ps_model):
    fPlin   = "./Interlopers/codes/fiducial/%s_matterpow_%s.dat"%(ps_model,z_med)
    kr, J0 = np.loadtxt('./Interlopers/codes/fiducial/J0filt.dat', unpack=True)
    k_NL, pm_NL = np.loadtxt(fPlin, unpack=True)
    pm_NL *= bias**2
    if rsd == 'rsd':
        omega_mz = omega_m*(1.+z_med)**3/(omega_m*(1.+z_med)**3+1.-omega_m)
        nk = np.size(k_NL)
        beta = omega_mz**0.55/bias
        sig_s = sig_s_fog
        rsd_factor = calc_stream(k_NL, nk, beta, sig_s, 'Gauss')
        pm_NL = pm_NL*rsd_factor
    xi_mod = p2xi_lin(k_NL, pm_NL, kr, J0, r_list) 
    return xi_mod



