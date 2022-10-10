import numpy as np
import os
import sys
def read_mean_data(sims, meanfile_name, int_x, frac, dirname, multi = False):
    path = dirname+'/%s/poles_%s_frac_%s.txt'%(meanfile_name,int_x, frac)
    if os.path.isfile(path):
        print('\n The file exists. sim = %s \n'%(meanfile_name))
        data = np.loadtxt(path)
        if multi == False:
            return data[:,0],data[:,1]
        else:
            return data[:,0], data[:,1], data[:,2]
    else:
        print('Creating the average file...')
        numsim = len(sims)
        path_solo = dirname+'/%s/poles_%s_frac_%s.txt'%(sims[0],int_x,frac)
        if os.path.isfile(path_solo) and len(sims) != 1000:
            poles0 = np.loadtxt(path_solo)
            poles_avg = np.zeros(poles0.shape)
            for sim in sims:
                poles = np.loadtxt(dirname+'/%s/poles_%s_frac_%s.txt'%(sim,int_x,frac))
                poles_avg += poles
            poles_avg /= numsim
            np.savetxt(path, poles_avg)
            if multi == False:
                return poles_avg[:,0],poles_avg[:,1]
            else:
                return poles_avg[:,0],poles_avg[:,1], poles_avg[:,2]
        else:
            if int_x == 'c':
                print(path)
                print('Warning: measured cross term file does not exist ... Won\'t plot that on the second page')
                if multi == False:
                    return 'False'
                else:
                    return 'False', 'False', 'False'
   
   

    
r"""
dirname = '/Users/setarehforoozan/Desktop/Interloper-project/codes/%s_%0.0f'%('norsd',97)
r_t, xi_t = read_mean_data(np.arange(0,10), 'means0to9', 't', 5, dirname)
r_c, xi_c = read_mean_data(np.arange(0,10), 'means0to9', 'c', 5, dirname)

alpha = 1.1
xi_t_alpha = np.interp(r_t*alpha, r_t, xi_t)
xi_c_alpha = np.interp(r_c*alpha, r_c, xi_c)
c1 = '#377eb8'
c2 = '#ff7f00'
import matplotlib.pyplot as plt
plt.plot(r_t,xi_t*r_t**2, label = r'auto, $\alpha = 1$', color = c1, ls = 'dashed')
plt.plot(r_t,xi_t_alpha*r_t**2, label = r'auto, $\alpha = %0.2f$'%alpha, color = c2, ls = 'dashed')

plt.plot(r_c,xi_c*r_c**2, label = r'cross, $\alpha = 1$', color = c1)
plt.plot(r_c,xi_c_alpha*r_c**2, label = r'cross, $\alpha = %0.2f$'%alpha, color = c2)
plt.ylabel(r'$r^2 \xi_0(r)$')
plt.xlabel(r'r(Mpc/h)')
plt.legend()

plt.show()"""