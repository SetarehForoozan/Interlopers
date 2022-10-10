import os
import numpy as np
from read_data import *
import default_config as cfg
import sys
import fitting_params as fprm

def read_covmat(ind_x, frac, nsim_covmat):
    #we have made the covariances in graham
    covpath = fprm.hhdir+'/CovMats/covmat_%s_frac_%s_nsim%s.txt'%(ind_x,frac,nsim_covmat)
    if os.path.isfile(covpath):
        print('\n Covariance matrix exists. \n')
        covmat = np.loadtxt(covpath)
    else:
        print('\n Covariance matrix Does not Exist. \n')
        sys.exit()
    return covmat/len(fprm.sims)

def read_covmat_multi(ind_x, frac, nsim_covmat, indmin, indmax, nr):
    #we have made the covariances in graham
    covpath = fprm.hhdir+'/CovMats/covmat_multi_%s_frac_%s_nsim%s_0to999.txt'%(ind_x,frac,nsim_covmat)
    if os.path.isfile(covpath):
        print('\n Covariance matrix exists. \n')
        covmat = np.loadtxt(covpath)
    else:
        print('\n Covariance matrix Does not Exist. \n')
        sys.exit()
    
    covmat1 = covmat[indmin:indmax,indmin:indmax] 
    covmat2 = covmat[indmin:indmax,indmin+nr:indmax+nr] 
    covmat3 = covmat[indmin+nr:indmax+nr,indmin:indmax] 
    covmat4 = covmat[indmin+nr:indmax+nr,indmin+nr:indmax+nr]
    
    covmat1 = np.concatenate((covmat1, covmat2), axis = 0)
    covmat3 = np.concatenate((covmat3, covmat4), axis = 0)
    covmatfin = np.concatenate((covmat1, covmat3), axis = 1)
    
    return covmatfin/len(fprm.sims)


def cov_to_err(cov):
    return np.sqrt(np.diagonal(cov))



r"""
import default_config as cfg
######This part of the code saves the covariance matrix and correlation matrix######
import matplotlib.pyplot as plt
font = {'size'   : 15 }
plt.rc('font', **font)

rstep = (150-0.5)/101
rmin = 0.5
rmax = 150
r_old = np.arange(rmin,rmax,rstep)
r_list = (r_old[1:]+r_old[:-1])/2.0 # central
rfitmin = 50
rfitmax = 150
r0_data = r_list
indmin, indmax = np.min(np.where(np.logical_and(r0_data >= rfitmin, r0_data <= rfitmax))[0]),np.max(np.where(np.logical_and(r0_data >= rfitmin, r0_data <= rfitmax))[0])
r = r_list[indmin:indmax]
xx, yy = np.meshgrid(r, r)

fig, axes = plt.subplots(2,2, figsize = (18,15), sharex = True, sharey = True)
fig.suptitle('Correlation Matrices of Auto-Correlations')
ind_x = 't'
for frac,i,j in zip([0,5,10,15],[0,0,1,1],[0,1,0,1]):
    covmat = read_covmat(ind_x, frac, 1000)
    covmat = covmat[indmin:indmax,indmin:indmax]
    cormat = np.copy(covmat)
    for k in range(len(r)):
        for p in range(len(r)):
            cormat[k,p] = covmat[k,p] / (covmat[k,k]*covmat[p,p])**(0.5)
    ax = axes[i,j]
    ax.set_title(r'$\xi_{%s}, f_i = %d$'%(ind_x, frac))
    pcm = ax.pcolormesh(xx,yy,cormat)
    fig.colorbar(pcm, ax=ax)
axes[0,0].set_ylabel(r'$r(Mpc/h)$')
axes[1,0].set_ylabel(r'$r(Mpc/h)$')
axes[1,1].set_xlabel(r'$r(Mpc/h)$')
axes[1,0].set_xlabel(r'$r(Mpc/h)$')
plt.savefig(cfg.hhdir+'/CorrelationMatrices_t.png')


fig, axes = plt.subplots(2,2, figsize = (18,15), sharex = True, sharey = True)
fig.suptitle('Correlation Matrices of Cross-Correlations')
ind_x = 'c'
for frac,i,j in zip([5,10,15],[0,1,1],[1,0,1]):
    covmat = read_covmat(ind_x, frac, 1000)
    covmat = covmat[indmin:indmax,indmin:indmax]
    cormat = np.copy(covmat)
    for k in range(len(r)):
        for p in range(len(r)):
            cormat[k,p] = covmat[k,p] / (covmat[k,k]*covmat[p,p])**(0.5)
    ax = axes[i,j]
    ax.set_title(r'$\xi_{%s}, f_i = %d$'%(ind_x, frac))
    pcm = ax.pcolormesh(xx,yy,cormat)
    fig.colorbar(pcm, ax=ax)
axes[0,0].set_ylabel(r'$r(Mpc/h)$')
axes[1,0].set_ylabel(r'$r(Mpc/h)$')

axes[1,1].set_xlabel(r'$r(Mpc/h)$')
axes[1,0].set_xlabel(r'$r(Mpc/h)$')

plt.savefig(cfg.hhdir+'/CorrelationMatrices.png')

"""