import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import json
#from read_data import *
#import default_config as cfg
l = 2
def Legendre(mu, l):
    coeff = (2*l+1)/2.
    if l == 0:
        Legendre = np.ones(shape = mu.shape)
    elif l == 2:
        Legendre = 1/2*(3*mu**2-1)
    elif l == 4:
        Legendre = 1/8*(35*mu**4-30*mu**2+3)
    return Legendre * coeff

def double_mirror(x):
    xx = np.concatenate((x[::,::-1],x), axis = 1)
    return xx

def mup_rp(mup, rp, A):
    r = (rp**2+A**2-2*A*rp*mup)**(1/2)
    mu = (rp*mup-A)/r
    return mu, r

def mapping(sim, frac, dirname, A, l):
    rp_list = np.linspace(0.5,150,100)
    bins_r = len(rp_list) #number of bins in r                                                                                                                                                                                                                                                          
    bins_mu = 300

    mu_arr = np.linspace(-1,1,int(bins_mu*2))
    mup_array = np.array([mu_arr]*bins_r).flatten()
    
    f = open(dirname+'/%s/xi_vm_%s_frac_%s.txt'%(sim, 't', frac))
    data = json.load(f)
    xi_arr = np.zeros((bins_mu, bins_r))
    Rp_array = np.zeros((bins_mu, bins_r))

    for mu_ind in range(bins_mu):
        x = []
        y = []
        for arr,ind_r in zip(data['corr']['__data__'],range(bins_r)):
            x += [arr[mu_ind][1]]
            y += [arr[mu_ind][0]]
        xi_arr[mu_ind,:] = y
        Rp_array[mu_ind,:] = x

    Rp_array = double_mirror(Rp_array.T).flatten()
    xi_corr = double_mirror(xi_arr.T).flatten()
    
    mu_array, R_array = mup_rp(mup_array, Rp_array, A)
    
    sortedinds = R_array.argsort()
    R_array  = R_array[sortedinds[::-1]]
    mu_array = mu_array[sortedinds[::-1]]
    xi_corr  = xi_corr[sortedinds[::-1]]
    
    R_array  = R_array.reshape(bins_r, bins_mu * 2)
    mu_array = mu_array.reshape(bins_r, bins_mu * 2)
    xi_corr  = xi_corr.reshape(bins_r, bins_mu * 2)

    xi_r = []
    for i in range(0, bins_r):
        sortinds = mu_array[i, :].argsort()
        R_array[i,:]  = R_array[i,sortinds]
        mu_array[i,:] = mu_array[i,sortinds]
        xi_corr[i,:]  = xi_corr[i,sortinds]
        mu_mid = (mu_array[i,1:]+mu_array[i,:-1])/2
        xi_r += [np.sum((mu_array[i,1:]-mu_array[i,:-1])*(xi_corr[i,1:]+xi_corr[i,:-1]) * Legendre(mu_mid, l))]
    xi_r = np.array(xi_r)/2 #divide by 2 because of mean of Xi_corr
    r_avg = np.mean(R_array, axis = 1)
    return r_avg[::-1], xi_r[::-1] 


sims = np.arange(400,600)
meanfile_name = 'mean%dto%d'%(np.min(sims), np.max(sims))
for rsd in ['norsd','rsd']:
    for A in [85,90,97.4]:
        dirname_f = '/home/setareh/projects/rrg-wperciva/setareh/results/%s_%0.0f'%(rsd,A)
        dirname_0 = '/home/setareh/projects/rrg-wperciva/setareh/results/%s_%0.0f'%(rsd,0)
        fracs = [0,2,5,10,15]
        for frac,u in zip(fracs,range(len(fracs))):
            if frac == 0: dirname = dirname_0
            else: dirname = dirname_f
            xc_exp = []
            yc_exp = []
            for sim in sims:
                xcn, ycn =  mapping(sim, frac, dirname, A, 2)
                xc_exp += [xcn]
                yc_exp += [ycn]
            xc_exp = np.average(xc_exp, axis = 0)
            yc_exp = np.average(yc_exp, axis = 0)
            data = np.concatenate((xc_exp.reshape(len(xc_exp),-1), yc_exp.reshape(len(xc_exp),-1)),axis =  1)
            np.savetxt('/home/setareh/scratch/results/%s_%d/%s/xi_cross_L%d_d%d_frac%d.txt'%(rsd,A,meanfile_name,l,A,frac),data)
            