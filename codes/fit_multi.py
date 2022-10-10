from default_config import *
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import emcee
import corner
from IPython.display import display, Math
import mono_functions as monofunc
import sys
from make_covmat import *
from xi_from_camb import *
from quad_functions import *

################## General plotting parameters #######################
plt.rc('lines', linewidth=2, markersize = 7)
plt.rc('font', size = 13)

color_fit   = '#ff7f00'
color_data  = '#377eb8'
color_poly  = '#e41a1c'
color_fid_1 = '#f781bf'
color_fid_2 = '#4daf4a'
posterior   = 'Percival'

theta_mono_0, theta_quad_0, params = set_initialparams(dof_a, dof_q, case)
parmax = len(params)
ndim = len(theta_mono_0)+len(theta_quad_0)

################## Read the cross terms and fiducial terms ##################
if case != 0:
    r_cross_L0, xi_cross_L0 = monofunc.xi_cross_func(frac,A,hhdir)
    r_cross_L2, xi_cross_L2 = monofunc.xi_cross_func(frac,A,hhdir,l=2) 
    r_cross_sim, xi_crossL0_sim, xi_crossL2_sim = read_mean_data(sims, meanfile_name, 'c', frac, hhdir, multi = True) #this is the average file

xi0_fid = np.loadtxt('./Interlopers/codes/xi0_cosmo_%s.txt'%(rsd))*bias**2
xi2_fid = np.loadtxt('./Interlopers/codes/xi2_cosmo_%s.txt'%(rsd))*bias**2
xi4_fid = np.loadtxt('./Interlopers/codes/xi4_cosmo_%s.txt'%(rsd))*bias**2

r0_data, xi0_data, xi2_data = read_mean_data(sims, meanfile_name, 't', frac, hhdir, multi = True)
if case!= 0:    
    r0_data_avg, xi0_avg0, xi2_avg0 =  read_mean_data(sims, meanfile_name, 't', 0, hhdir, True) #this is the average file for f = 0
 
"""import xi_from_camb as xi
xi_mono = xi.xi_camb_func('CAMB')
xi0_fid_mine = np.loadtxt('/Users/setarehforoozan/Desktop/Interloper-project/codes/new_fitting_codes/xi%d_%s_fid.txt'%(0,rsd))[:,1]*bias**2
plt.plot(r_list, r_list**2*xi_mono, label = 'Monpole code')
plt.plot(r_list, r_list**2*xi0_fid, label = 'Cosmocode Multipole')
plt.plot(r_list, r_list**2*xi0_fid_mine, label = 'My multipole code')
plt.legend()
plt.show()
sys.exit()"""

################## Set the saving locations ##################
specific_name = 'frac%d_dega%d_dega%d_sims%dto%d_%dto%d'%(frac, dof_a, dof_q, np.min(sims), np.max(sims), rfitmin,rfitmax)
casename = 'multi_case%d_%dsims'%(case, len(sims))
path1 = hhdir+'/'+casename+'/'
path2 = hhdir+'/'+casename+'/AVG-STD'
path3 = hhdir+'/'+casename+'/paper'

make_path(path1)
make_path(path2)
make_path(path3)
    
figout = PdfPages(hhdir+"/%s/%s.pdf"%(casename,specific_name))
fout_avg_name = hhdir+"/%s/AVG-STD/AVG_%s.txt"%(casename,specific_name)
fout_std_name = hhdir+"/%s/AVG-STD/STD_%s.txt"%(casename,specific_name)
fout_final_name = hhdir+"/%s/paper/FINAL_%s.txt"%(casename,specific_name)
fout_avg = open(fout_avg_name, "w")
fout_std = open(fout_std_name, "w")

############################# READ THE DATA AND FIDUCAL IN R FIT RANGE ###################################

indmin, indmax = np.min(np.where(np.logical_and(r0_data >= rfitmin, r0_data <= rfitmax))[0]),np.max(np.where(np.logical_and(r0_data >= rfitmin, r0_data <= rfitmax))[0])
covmat = read_covmat_multi('t', frac, nsim_covmat, indmin, indmax, 100) #change this if you change the number of bins in r
r0_data = r0_data[indmin:indmax]
xi0_data = xi0_data[indmin:indmax]
xi2_data = xi2_data[indmin:indmax]
xi0_data_error = cov_to_err(covmat)[0:len(r0_data)]
xi2_data_error = cov_to_err(covmat)[len(r0_data):]
xi_data = np.concatenate((xi0_data.reshape(len(xi0_data), -1),xi2_data.reshape(len(xi2_data), -1)),axis = 1)

xi0_fid = np.interp(r0_data, r_list, xi0_fid)
xi2_fid = np.interp(r0_data, r_list, xi2_fid)
xi4_fid = np.interp(r0_data, r_list, xi4_fid)
xi_fid = np.concatenate((xi0_fid.reshape(len(xi0_fid), -1),xi2_fid.reshape(len(xi2_fid), -1)),axis = 1)
xi_fid = np.concatenate((xi_fid,xi4_fid.reshape(len(xi4_fid ), -1)),axis = 1)

if case != 0:
    xi_cross_L0 = np.interp(r0_data, r_cross_L0, xi_cross_L0)
    xi_cross_L2 = np.interp(r0_data, r_cross_L2, xi_cross_L2)
    xi_cross = np.concatenate((xi_cross_L0.reshape(len(xi_cross_L0), -1),xi_cross_L2.reshape(len(xi_cross_L2), -1)),axis = 1)
else:
    xi_cross = 0
mu_fid, weight = np.loadtxt("./Interlopers/codes/gauss.txt", unpack=True)

"""
#plot cormat and investigate whether binning is right
cormat = np.copy(covmat)
double_r = np.concatenate((r0_data,r0_data+r0_data[-1]-r0_data[0]), axis =0)
xx = np.zeros(shape = covmat.shape)
yy = np.zeros(shape = covmat.shape)
for k in range(len(double_r)):
    for p in range(len(double_r)):
        cormat[k,p] = covmat[k,p] / (covmat[k,k]*covmat[p,p])**(0.5)
        xx[k,p] = double_r[k]
        yy[k,p] = double_r[p]

pcm = plt.pcolormesh(xx,yy,cormat)
plt.show()
sys.exit
plt.errorbar(r0_data, r0_data**2 * xi0_data, yerr = r0_data**2* xi0_data_error, label = 'Data', marker = 'o', ms = 2)
plt.show()
plt.errorbar(r0_data, r0_data**2 * xi2_data, yerr = r0_data**2* xi2_data_error, label = 'Data', marker = 'o', ms = 2)
plt.show()
sys.exit()"""
    
r"""####plot the thing    
plt.scatter(r0_data, xi0_data * r0_data**2, label = 'Data')
theta_fiducial = [0,0,0,1,1,frac/100]
fiducial_model = xi0_model(theta_fiducial, r0_data, xi_mod, r_cross, xi_cross)[0]
plt.plot(r0_data, fiducial_model * r0_data**2, label = 'Fiducial (0,0,0,1,1,%d)'%(frac))
plt.legend()
plt.ylabel(r'$r^2 \xi_0(r)$')
plt.xlabel(r'r(Mpc/h)')
plt.show()    
sys.exit()"""

####################### MONTE CARLO PREPARATION #########################

theta_0_walkers = np.array(np.concatenate((theta_mono_0,theta_quad_0), axis = 0)) + 1e-2 * np.random.randn(nwalkers, ndim)
#mu_fid in the next line is not necessary since it's always read from the quad_functions file,x, weights.
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(r0_data, mu_fid, weight, xi_fid, xi_cross, xi_data, covmat, case, posterior))
sampler.run_mcmc(theta_0_walkers, nsteps_multi, progress=True)
samples = sampler.get_chain()
flat_samples = sampler.get_chain(discard=ndiscard_multi, thin=15, flat=True)

#####################  PLOT THE FIT AND DATA  ############################
fig, axes = plt.subplots(1, figsize = (10,7))
theta_avg = np.zeros(ndim)
theta_median = np.zeros(ndim)
theta_std = np.zeros(ndim)
plt.errorbar(r0_data, r0_data**2 * xi0_data, yerr = r0_data**2* xi0_data_error, label = 'Data', marker = 'o', ms = 2, color = color_data)
for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    theta_median[i] = mcmc[1]
    theta_avg[i] = np.average(flat_samples[:, i])
    theta_std[i] = np.std(flat_samples[:, i])
    txt_meidan = r'%s'%(params[i])+r'$= %0.3f_{-%0.3f}^{%0.3f}$'%(mcmc[1], q[0], q[1])
    txt_avg = r'%s'%(params[i])+r'$= %0.3f \pm %0.3f$'%(theta_avg[i], theta_std[i])
    fout_avg.write(str(r'%0.5f'%(theta_avg[i]))+' ')
    fout_std.write(str(r'%0.5f'%(theta_std[i]))+' ')
    plt.gcf().text(0.75,0.7-0.03*i, txt_avg, fontsize = 14)

if case != 0: uuu = 3 
else: uuu = 2
theta_mono_avg = theta_avg[:dof_a+uuu]
theta_quad_avg = theta_avg[dof_a+uuu:]
chi2 = chi2_func(theta_mono_avg, theta_quad_avg, r0_data, mu_fid, weight, xi_fid, xi_cross, xi_data, covmat, case)
chi2_reduced = chi2/(len(r_list)-ndim)
plt.gcf().text(0.75,0.7-0.03*(i+1), r'$\chi^2$: %0.2f'%(chi2), fontsize = 14)
plt.gcf().text(0.75,0.7-0.03*(i+2), r'Reduced $\chi^2$: %0.3f'%(chi2_reduced), fontsize = 14)    
plt.subplots_adjust(right = 0.7)

if case != 0:
    xi0_fit, xi2_fit, xicL0_fit, xicL2_fit, xi0_temp, xi2_temp = xi_poles_obs(theta_mono_avg, theta_quad_avg, r0_data, mu_fid, weight, xi_fid, xi_cross, case)
    plt.plot(r0_data, r0_data**2 * xi0_fit ,label = 'Fit', marker = 'o', ms = 4, color = color_fit)
else:
    xi0_fit, xi2_fit = xi_poles_obs(theta_mono_avg, theta_quad_avg, r0_data, mu_fid, weight, xi_fid, xi_cross, case)
    plt.plot(r0_data, r0_data**2 * xi0_fit ,label = 'Fit', marker = 'o', ms = 4, color = color_fit)

plt.plot(r0_data, r0_data**2 * xi0_fid ,label = 'Template', marker = 'o', ms = 4, color = 'black')
plt.plot(r0_data, r0_data**2* polynomial_a(theta_mono_avg,r0_data), label = 'Poly term', color = color_poly)

if case !=0: 
    theta_fid_m = np.array([1,1,frac/100]) #b,alpha,frac
else:
    theta_fid_m = np.array([1,1]) #b,alpha
    
theta_fid_m_1 = np.insert(theta_fid_m, 0, np.zeros(dof_a))  #a1,a2,a3 = 0
theta_fid_q = np.array([0]) #epsilon
theta_fid_q_1 = np.insert(theta_fid_q, 0, np.zeros(dof_q))  #q1,q2,q3 = 0
if case !=0: 
    theta_fid_m = np.array([1,frac/100]) #alpha,frac
else:
    theta_fid_m = np.array([1]) #alpha,frac
theta_fid_m_2 = np.insert(theta_fid_m, 0, theta_mono_avg[:dof_a+1]) #a1,a2,a3 = a1,a2,a3 best
theta_fid_q = np.array([0]) #epsilon 
theta_fid_q_2 = np.insert(theta_fid_q, 0, theta_quad_avg[:dof_q]) #a1,a2,a3 = a1,a2,a3 best

xi_fiducial_1 = xi_poles_obs(theta_fid_m_1, theta_fid_q_1, r0_data, mu_fid, weight, xi_fid,xi_cross, case)[0]
xi_fiducial_2 = xi_poles_obs(theta_fid_m_2, theta_fid_q_2, r0_data, mu_fid, weight, xi_fid,xi_cross, case)[0]

plt.plot(r0_data,r0_data**2 * xi_fiducial_1, label = r'Fid ($A=0,B=1$)', color = color_fid_1)
plt.plot(r0_data,r0_data**2 * xi_fiducial_2, label = r'Fid ($A=B=best$)', color = color_fid_2)
plt.title(r'Monopole of the Auto Correlation, $f_i^{%s} = %d$%s'%('true', frac, '%'))
plt.legend()
plt.ylabel(r'$r^2 \xi_0(r)$')
plt.xlabel(r'r(Mpc/h)')
plt.xlim(rfitmin,rfitmax)
figout.savefig(fig)
plt.clf()    
fout_avg.write('\n')
fout_std.write('\n')

"""####################### FIT FOR ONLY POLYNOMIALS ########################
theta_3 = theta_avg[:dof_a] #a1,a2,a3 = a
fid_theta_0_walkers = np.array(theta_3) + 1e-2 * np.random.randn(nwalkers, dof_a)
fid_sampler = emcee.EnsembleSampler(nwalkers, dof_a, log_probability, args=(r0_data, xi0_data, covmat, xi_mod, r_cross, xi_cross,True))
fid_sampler.run_mcmc(fid_theta_0_walkers, nsteps, progress=True)
fid_samples = fid_sampler.get_chain()
fid_flat_samples = fid_sampler.get_chain(discard=ndiscard, thin=15, flat=True)
fid_theta_avg = np.zeros(dof_a)
for i in range(dof_a):
    fid_theta_avg[i] = np. average(fid_flat_samples[:, i])
xi_fiducial = xi0_model_fix(fid_theta_avg,r0_data,xi_mod,r_cross, xi_cross)
plt.plot(r0_data,r0_data**2 * xi_fiducial, label = 'Fid+Fit Poly')
"""
    
################## Plot the quadrupole ##############
plt.errorbar(r0_data, xi2_data, yerr = xi2_data_error, label = 'Data', marker = 'o', ms = 2,color = color_data)
plt.plot(r0_data, xi2_fit, label = 'Fit', marker = 'o', ms = 4, color = color_fit)
plt.plot(r0_data, polynomial_q(theta_quad_avg,r0_data), label = 'Poly term', color = color_poly)

xi_fiducial_1 = xi_poles_obs(theta_fid_m_1, theta_fid_q_1, r0_data, mu_fid, weight, xi_fid,xi_cross, case)[1]
xi_fiducial_2 = xi_poles_obs(theta_fid_m_2, theta_fid_q_2, r0_data, mu_fid, weight, xi_fid,xi_cross, case)[1]

plt.plot(r0_data,xi_fiducial_1, label = r'Fid ($A=0,B=1$)', color = color_fid_1)
plt.plot(r0_data,xi_fiducial_2, label = r'Fid ($A=B=best$)', color = color_fid_2)
#plt.plot(r0_data,r0_data**2 * xi2_fid, label = r'Template')
 
plt.legend()
plt.title(r'Quadrupole of the Auto Correlation, $f_i^{%s} = %d$%s'%('true', frac, '%'))
plt.ylabel(r'$\xi_2(r)$')
plt.xlabel(r'r(Mpc/h)')
figout.savefig(fig)
plt.clf()

###################### PLOT L0 CROSS TERMS ############################
if case != 0:
    fig, axes = plt.subplots(1, figsize = (10,7))
    if r_cross_sim != 'False':
        plt.plot(r_cross_sim, xi_crossL0_sim*r_cross_sim**2, label ='Data', color = color_data)
    plt.plot(r0_data, xi_cross_L0*r0_data**2,label = r'$M(\xi)_0$', color = color_fid_1, ls = 'dashed')
    plt.plot(r0_data, r0_data**2 * xicL0_fit, ls = 'dashed', label =  r'$M(\xi)_{0}/(1+2f_i^2-2f_i)$', color = color_fit)
    plt.title(r'Monopole of the Cross Correlation, $f_i^{%s} = %d$%s'%('true', frac, '%'))
    plt.ylabel(r'$r^2 \xi_{0,hi}(r)$')
    plt.xlabel(r'r(Mpc/h)')
    plt.xlim((rfitmin,rfitmax))
    plt.legend()
    figout.savefig(fig)
    plt.clf()

###################### PLOT L2 CROSS TERMS ############################    
if case != 0:
    fig, axes = plt.subplots(1, figsize = (10,7))
    if r_cross_sim != 'False':
        plt.plot(r_cross_sim, xi_crossL2_sim, label ='Data', color = color_data)
    plt.plot(r0_data, xi_cross_L2,label = r'$M(\xi)_2$', color = color_fid_1, ls = 'dashed')
    plt.plot(r0_data, xicL2_fit, ls = 'dashed', label =  r'$M(\xi)_{2}/(1+2f_i^2-2f_i)$', color = color_fit)
    plt.title(r'Quadrupole of the Cross Correlation, $f_i^{%s} = %d$%s'%('true', frac,'%'))
    plt.ylabel(r'$\xi_{2,hi}(r)$')
    plt.xlabel(r'r(Mpc/h)')
    plt.xlim((rfitmin,rfitmax))
    plt.legend()
    figout.savefig(fig)
    plt.clf()

if case != 0: 
    ####################### AUTO CORRS FOR UNCONCTAMINATED ###################
    fig, axes = plt.subplots(1, figsize = (10,7))
    plt.plot(r0_data, xi0_fid*r0_data**2,  ls = 'dashed',label=r'CAMB, $f^{\rm true} = 0$', color = color_fid_1)
    plt.plot(r0_data_avg, xi0_avg0*r0_data_avg**2, label=r'Data, $f^{\rm true} = 0$', color = color_data)
    xi_est = xi0_temp
    plt.plot(r0_data,r0_data**2* xi_est, label = 'Estimated', color = color_fid_1 )
    plt.plot(r0_data,r0_data**2* polynomial_a(theta_mono_avg,r0_data), label = 'Poly term', color = color_poly)
    plt.xlim((rfitmin,rfitmax))
    plt.legend()
    plt.title(r'Monopole of the Auto Correlation, Uncontaminated')
    plt.ylabel(r'$r^2 \xi_0(r)$')
    plt.xlabel(r'r(Mpc/h)')
    figout.savefig(fig)
    plt.clf()

    ####################### QUADRUPOLE FOR UNCONTAMINATED ###################
    fig, axes = plt.subplots(1, figsize = (10,7))
    plt.plot(r0_data, xi2_fid,  ls = 'dashed',label=r'CAMB, $f^{\rm true} = 0$', color = color_fid_1)
    plt.plot(r0_data_avg, xi2_avg0, label=r'Data, $f^{\rm true} = 0$', color = color_data)

    xi_est = xi2_temp
    plt.plot(r0_data, xi_est, label = 'Estimated', color = color_fid_1)
    plt.plot(r0_data, polynomial_q(theta_quad_avg,r0_data), label = 'Poly term', color = color_poly)
    plt.xlim((rfitmin,rfitmax))
    plt.ylim((-0.04,0.04))
    plt.legend()
    plt.title(r'Quadrupole of the Auto Correlation, Uncontaminated')
    plt.ylabel(r'$\xi_2(r)$')
    plt.xlabel(r'r(Mpc/h)')
    figout.savefig(fig)
    plt.clf()

###################### PLOT THE WALKERS############################
fig, axes = plt.subplots(ndim, sharex=True, figsize = (10,7))
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(params[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("Step Number")
figout.savefig(fig)
plt.clf()

###################### PLOT THE CORNER PLOTS ############################
fig = corner.corner(flat_samples, labels=params, truths=theta_median, figsize = (3,3))
figout.savefig(fig)
plt.clf()    
    
fout_avg.close()
fout_std.close()
figout.close()
print('Done')

avg_read = np.loadtxt(fout_avg_name)
std_read = np.loadtxt(fout_std_name)
fout_final = open(fout_final_name, "w")

def sci_not(v,err,rnd=1):
    power = int(('%E' % v)[-3:])
    return '({0} \pm {1})e{2}'.format(
            round(v/10**power,rnd),round(err/10**power,rnd),power)

avg = np.copy(avg_read)
std = np.copy(std_read)

for i in range(parmax):
    if 0+dof_a <= i < 0+3: #or dof_c+1 <= i < 4:
        avg = np.insert(avg, i, np.inf)       
        std = np.insert(std, i, np.inf)       
    if 3+uuu + dof_q <= i < 6+uuu: #or dof_c+1 <= i < 4:
        avg = np.insert(avg, i, np.inf)       
        std = np.insert(std, i, np.inf)  

if case != 0:
    a1,a2,a3,b,al,fi,w1,w2,w3,ep = avg
    a1_err,a2_err,a3_err,b_err,al_err,fi_err,w1_err,w2_err,w3_err,ep_err = std
    sort_par = np.array([al, ep, fi, b, a1, a2, a3, w1, w2, w3])
    sort_par_err = np.array([al_err, ep_err, fi_err, b_err, a1_err, a2_err, a3_err, w1_err, w2_err, w3_err])
else:
    a1,a2,a3,b,al,w1,w2,w3,ep = avg
    a1_err,a2_err,a3_err,b_err,al_err,w1_err,w2_err,w3_err,ep_err = std
    sort_par = np.array([al, ep, b, a1, a2, a3, w1, w2, w3])
    sort_par_err = np.array([al_err, ep_err, b_err, a1_err, a2_err, a3_err, w1_err, w2_err, w3_err])

for i, v, err in zip(range(parmax),sort_par, sort_par_err):
    fout_final.write('$')
    if np.isfinite(v):
        if i != 2:
            fout_final.write(sci_not(v,err,rnd=2))
        if i == 2 and case != 0:
            fout_final.write('%0.1f'%(v*100) + ' \pm ' + '%0.1f'%(err*100))
    else:
        fout_final.write('-')
    fout_final.write('$')
    fout_final.write('' + '\n')
#chi2
fout_final.write('$%0.2f$'%(chi2))        
fout_final.write('\\'+ '\\' '\n')
fout_final.write('$%0.3f$'%(chi2_reduced))
fout_final.write('\\'+ '\\' '\n')
fout_final.close()
