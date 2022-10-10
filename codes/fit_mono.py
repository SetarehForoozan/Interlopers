from default_config import *
import os
import matplotlib.pyplot as plt
import math
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import emcee
import corner
from IPython.display import display, Math
from mono_functions import *
import sys
from make_covmat import *
from xi_from_camb import *

if case != 0:
    params = ["$b$", r"$\alpha$", r"$f_{i}$"]
    theta_0 = [0.8, 1.0, 0.1]
    parmax = 6
else: 
    params = ["$b$", r"$\alpha$"]
    theta_0 = [0.8, 1.0]
    parmax = 5
    
for i in range(1,dof_a+1):
    title = str(r"$a_%d$"%(i))
    params = np.insert(params, i-1, title)
    if i == 1:
        theta_0 = np.insert(theta_0, i-1, 0)
    if i == 2:
        theta_0 = np.insert(theta_0, i-1, 0)
    if i == 3:
        theta_0 = np.insert(theta_0, i-1, 1e-3)

xi_camb = xi_camb_func('CAMB')
xi_camb_NL = xi_camb_NL_func('CAMB_NL')
#if case != 0 : r_cross, xi_cross = xi_cross_func(frac,A,hhdir) 
r_cross, xi_cross = xi_cross_func(frac,A,hhdir) 
ndim = len(params) #for a1, a2, a3, b, alpha and fi

################## Set the saving locations ##################
specific_name = 'frac%d_dega%d_sims%dto%d'%(frac, dof_a, np.min(sims), np.max(sims))
casename = 'mono_case%d_%dsims'%(case, len(sims))
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
################################################################



plt.rc('lines', linewidth=2, markersize = 7)
plt.rc('font', size = 13)

############################# MAIN PART OF THE CODE###################################

int_x = 't'
r0_data, xi0_data = read_mean_data(sims, meanfile_name, int_x, frac, hhdir)
covmat = read_covmat('t', frac, nsim_covmat)

indmin, indmax = np.min(np.where(np.logical_and(r0_data >= rfitmin, r0_data <= rfitmax))[0]),np.max(np.where(np.logical_and(r0_data >= rfitmin, r0_data <= rfitmax))[0])
#plt.plot(r0_data, (xi0_data-np.interp(r0_data, r_list, xi_mod_NL))*r0_data**2)
r0_data = r0_data[indmin:indmax]
xi0_data = xi0_data[indmin:indmax]
covmat = covmat[indmin:indmax,indmin:indmax]

"""
xi0_data_error = cov_to_err(covmat)
plt.errorbar(r0_data, r0_data**2 * xi0_data, yerr = r0_data**2* xi0_data_error, label = 'Data', marker = 'o', ms = 2)
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
sys.exit()    """

theta_0_walkers = np.array(theta_0) + 1e-2 * np.random.randn(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(r0_data, xi0_data, covmat, xi_camb, r_cross, xi_cross))
sampler.run_mcmc(theta_0_walkers, nsteps_mono, progress=True)
samples = sampler.get_chain()
flat_samples = sampler.get_chain(discard=ndiscard_mono, thin=15, flat=True)

#####################  PLOT THE FIT AND DATA  ############################
fig, axes = plt.subplots(1, figsize = (10,7))
theta_avg = np.zeros(ndim)
theta_median = np.zeros(ndim)
theta_std = np.zeros(ndim)
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
chi2 = -2 * log_likelihood(theta_avg, r0_data, xi0_data, covmat, xi_camb, r_cross, xi_cross)
chi2_reduced = chi2/(len(r_list)-ndim)
plt.gcf().text(0.75,0.7-0.03*(i+1), r'$\chi^2$: %0.2f'%(chi2), fontsize = 14)
plt.gcf().text(0.75,0.7-0.03*(i+2), r'Reduced $\chi^2$: %0.3f'%(chi2_reduced), fontsize = 14)    
plt.subplots_adjust(right = 0.7)
xi0_data_error = cov_to_err(covmat)
plt.errorbar(r0_data, r0_data**2 * xi0_data, yerr = r0_data**2* xi0_data_error, label = 'Data', marker = 'o', ms = 2)

if case == 0: 
    plt.plot(r0_data, r0_data**2 * xi0_model_nonfrac(theta_avg, r0_data, xi_camb) ,label = 'fit')
else: 
    plt.plot(r0_data, r0_data**2 * xi0_model(theta_avg, r0_data, xi_camb, r_cross, xi_cross)[0] ,label = 'fit')

plt.plot(r0_data,r0_data**2* polynomial_a(theta_avg,r0_data), label = 'Poly term')

if case !=0: 
    theta_fid = np.array([1,1,frac/100]) #b,alpha,frac
    theta_fid_1 = np.insert(theta_fid, 0, np.zeros(dof_a))  #a1,a2,a3 = 0
    theta_fid = np.array([1,frac/100]) #alpha,frac
    theta_fid_2 = np.insert(theta_fid, 0, theta_avg[:dof_a+1]) #a1,a2,a3 = a1,a2,a3 best
    xi_fiducial_1 = xi0_model(theta_fid_1, r0_data, xi_camb, r_cross, xi_cross)[0]
    xi_fiducial_2 = xi0_model(theta_fid_2, r0_data, xi_camb, r_cross, xi_cross)[0]

elif case == 0:
    theta_fid = np.array([1.00,1.00]) #b,alpha
    theta_fid_1 = np.insert(theta_fid, 0, np.zeros(dof_a))  #a1,a2,a3 = 0, b = alpha = 1
    print(theta_fid_1)
    theta_fid = np.array([1.00]) #alpha
    theta_fid_2 = np.insert(theta_fid, 0, theta_avg[:dof_a+1]) #a1,a2,a3,b = a1,a2,a3,b best and alpha = 1
    print(theta_fid_2)
    xi_fiducial_1 = xi0_model_nonfrac(theta_fid_1, r0_data, xi_camb)
    xi_fiducial_2 = xi0_model_nonfrac(theta_fid_2, r0_data, xi_camb)

plt.plot(r0_data,r0_data**2 * xi_fiducial_1, label = 'Fiducial')
plt.plot(r0_data,r0_data**2 * xi_fiducial_2, label = 'Fid+A+B')

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
    
plt.legend()
plt.ylabel(r'$r^2 \xi_0(r)$')
plt.xlabel(r'r(Mpc/h)')
figout.savefig(fig)
plt.clf()    
fout_avg.write('\n')
fout_std.write('\n')


###################### PLOT CROSS TERMS ############################
if case != 0:
    #### CROSS CORRS ####
    fig, axes = plt.subplots(1, figsize = (7,7))
    plt.plot(r_cross, xi_cross*r_cross**2,label = r'Cross: M($\xi$)')
    xi_est = xi0_model(theta_avg, r0_data, xi_camb, r_cross, xi_cross)[1]
    plt.plot(r0_data, r0_data**2 * xi_est, ls = 'dashed', label =  r'Cross: M($\xi$) / p')
    
    if read_mean_data(sims, meanfile_name, 'c', 10, hhdir, multi = False) != 'False':
        r_cross_sim, xi_cross_sim =  read_mean_data(sims, meanfile_name, 'c', 10, hhdir, multi = False) #this is the average file
        plt.plot(r_cross_sim, xi_cross_sim*r_cross_sim**2, label ='Cross: CAT')
    plt.ylabel(r'$r^2 \xi_0(r)$')
    plt.xlabel(r'r(Mpc/h)')
    plt.xlim((0,150))
    plt.legend()
    figout.savefig(fig)
    plt.clf()

####################### AUTO CORRS ###################
fig, axes = plt.subplots(1, figsize = (7,7))
plt.plot(r_list, xi_camb*r_list**2,  ls = 'dashed',label='Auto: CAMB, f = 0')
#plt.plot(r_list, xi_mod_NL*r_list**2, ls = 'dotted' , label='Auto: CAMB_NL, f = 0')
r0_data, xi_avg =  read_mean_data(sims, meanfile_name, 't', 0, hhdir) #this is the average file
plt.plot(r0_data, xi_avg*r0_data**2, label='Auto: CAT, f = 0')
r0_data, xi_avg =  read_mean_data(sims, meanfile_name, 't', frac ,hhdir)  #this is the average file
#plt.plot(r0_data, xi_avg*r0_data**2, label='Auto: CAT, f = frac')
if case == 0: 
    xi_est = xi0_model_nonfrac(theta_avg, r0_data, xi_camb)
else: 
    xi_est = xi0_model(theta_avg, r0_data, xi_camb, r_cross, xi_cross)[2]
plt.plot(r0_data,r0_data**2* xi_est, label = 'Auto: Estimated')
plt.plot(r0_data,r0_data**2* polynomial_a(theta_avg,r0_data), label = 'Poly term')
plt.xlim((0,150))
plt.legend()
plt.ylabel(r'$r^2 \xi_0(r)$')
plt.xlabel(r'r(Mpc/h)')
figout.savefig(fig)
plt.clf()

###################### PLOT THE WALKERS############################
fig, axes = plt.subplots(ndim, sharex=True, figsize = (7,7))
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
    return '[{0} \pm {1}]e{2}'.format(
            round(v/10**power,rnd),round(err/10**power,rnd),power)

avg = np.copy(avg_read)
std = np.copy(std_read)

for i in range(parmax):
    if 1 + dof_a <= i < 4: #or dof_c+1 <= i < 4:
        avg = np.insert(avg, i, np.inf)       
        std = np.insert(std, i, np.inf)       

for i in range(parmax):
    if np.isfinite(avg[i]):
        if i == 0: #bc
            fout_final.write('%0.2f'%(avg[i]) + ' \pm ' + '%0.2f'%(std[i]))
            fout_final.write('' + '\n')

        if i == 1: #a1
            fout_final.write('%0.1f'%(avg[i]) + ' \pm ' + '%0.1f'%(std[i]))
            fout_final.write('' + '\n')

        elif i == 2: #a2
                fout_final.write('%0.2f'%(avg[i]) + ' \pm ' + '%0.2f'%(std[i]))
                fout_final.write('' + '\n')

        elif i == 3:#a3
                fout_final.write('%0.5f'%(avg[i]) + ' \pm ' + '%0.5f'%(std[i]))
                fout_final.write('' + '\n')

        elif i == 4: #b
            fout_final.write('%0.2f'%(avg[i]) + ' \pm ' + '%0.2f'%(std[i]))            
            fout_final.write('' + '\n')

        elif i == 5: #alpha
            fout_final.write('%0.4f'%(avg[i]) + ' \pm ' + '%0.4f'%(std[i]))
            fout_final.write('' + '\n')

        elif i == 6: #f_i
            fout_final.write('%0.1f'%(avg[i]*100) + ' \pm ' + '%0.1f'%(std[i]*100))
            fout_final.write('' + '\n')

    else:
        fout_final.write('-')
        fout_final.write('' + '\n')

#chi2        
fout_final.write('$%0.3f$'%(chi2_reduced))
fout_final.write('\\'+ '\\' '\n')
fout_final.close()
