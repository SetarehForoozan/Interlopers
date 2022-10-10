import numpy as np
from spl import *
from default_config import *
from read_data import *
from fitting_params import *

def ratio_func(f):
    """ _summary_
    This is the prefactor in front of the cross-correlation that comes from auto-correlation
    Returns: prefactor
        _type_: float
    """
    f = f
    return (1+2*f**2-2*f)

def xi_cross_func(frac,A,hhdir, l = 0):
    rp_list = np.linspace(0.5,150,101)
    if case == 1: #### cross term read form directly mesuring from the simulations
        r_cross, xi_cross = read_mean_data(sims, meanfile_name, 'c', frac, hhdir)
        return r_cross, xi_cross

    if case == 6 or case == 8: #### cross term from convolution of the auto correlation
        if l == 0: 
            path = hhdir+'/%s/xi_cross_d%d_frac%d.txt'%(meanfile_name, A,frac)
            isFile = os.path.isfile(path)
            if isFile:
                data = np.loadtxt(path)
            else:
                data = np.loadtxt(hhdir+'/%s/xi_cross_L0_d%d_frac%d.txt'%(meanfile_name, A,frac))
        if l == 2: data = np.loadtxt(hhdir+'/%s/xi_cross_L2_d%d_frac%d.txt'%(meanfile_name, A,frac))
        xi0_read = np.interp(rp_list, data[:,0], data[:,1])
        return rp_list, xi0_read

    elif case == 0: #### no cross term
        return 0,0

def Legendre(mu, l):
    if l == 0:
        Legendre = np.ones(shape = mu.shape)
    elif l == 2:
        Legendre = 1/2*(3*mu**2-1)
    elif l == 4:
        Legendre = 1/8*(35*mu**4-30*mu**2+3)
    return Legendre

def xi0_model(theta, r, xi_mod_norm, r_cross, xi_cross):
    fi = theta[-1]
    alpha = theta[-2]
    b = theta[-3]
    xi0_auto = np.interp(r*alpha, r_list, xi_mod_norm) 
    if case == 8:
        xi0_cross = np.interp(r, r_cross, xi_cross)
        bc = 1/ratio_func(fi)
    elif case == 6:
        xi0_cross = np.interp(r*alpha, r_cross, xi_cross)
        bc = 1/ratio_func(fi)
    else:
        print('model not defined for this case')
        sys.exit()
    poly_a = polynomial_a(theta, r)
    
    xi_auto_final = b * xi0_auto + poly_a
    xi_cross_final = bc * xi0_cross
    xi0_eff = (1+2*fi**2-2*fi) * xi_auto_final + 2*fi*(1-fi) * xi_cross_final
    return [xi0_eff,xi_cross_final,xi_auto_final]

def xi0_model_nonfrac(theta, r, xi_mod_norm):
    alpha = theta[-1]
    b = theta[-2]
    xi0_auto = np.interp(r*alpha, r_list, xi_mod_norm) 
    poly_a = polynomial_a(theta, r)
    xi_auto_final = b * xi0_auto + poly_a
    return  xi_auto_final

def polynomial_a(theta, r):
    """_summary_
    Args:
        theta_quad (array): an array consisting of constants for monopole: 
            dof_a = 3: [a1,a2,a3]
            dof_a = 2: [a1,a2]
            dof_a = 1: [a1]
        r (array): the r array 

    Returns:
        array with len(r): the sum of all polynomials
    """
    vec_a = theta[:dof_a]
    poly_a = 0
    for i in range(1,dof_a+1):
        poly_a +=  vec_a[i-1] * r**(-(3-i))
    return poly_a
               
def log_likelihood(theta, r_data, xi_data, covmat, xi_mod_norm, r_cross, xi_cross, fix = False):
    if case == 0:
        model = xi0_model_nonfrac(theta, r_data, xi_mod_norm)
    else: 
        model = xi0_model(theta, r_data, xi_mod_norm, r_cross, xi_cross)[0]
    diff = xi_data - model
    covinv = np.linalg.inv(covmat)
    chi2 = (diff.dot(covinv)).dot(diff)
    ns = 1000
    nd = len(r_list)
    nt = len(theta)
    B = (ns-nd-2)/(ns-nd-1)/(ns-nd-4)
    m = nt + 2 + (ns - 1 + B*(nd-nt))/(1+B*(nd-nt))
    post = -m/2*np.log((1+chi2/(ns-1)))
    return post

def log_prior(theta):
    fi = theta[-1]
    alpha = theta[-2]
    b = theta[-3]
    vec_a = theta[:dof_a]
    if 0.6 < alpha < 1.4 and 0 < fi < 0.5: 
        return 0.0
    return -np.inf

def log_prior_nonfrac(theta):
    alpha = theta[-1]
    b = theta[-2]
    vec_a = theta[:dof_a]
    vec_c = theta[:]
    if 0.6 < alpha < 1.4:
        return 0.0
    return -np.inf

def log_probability(theta, r_data, xi_data, covmat, xi_mod_norm, r_cross, xi_cross, fix = False):
    if case != 0:
        lp = log_prior(theta)
    else:
        lp = log_prior_nonfrac(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, r_data, xi_data, covmat, xi_mod_norm, r_cross, xi_cross, fix)


def make_path(path):
    isFile = os.path.isdir(path)
    if not(isFile):
        print('path \n"%s"\n did not exist. Now creating...'%(path))
        os.mkdir(path)
    return 0