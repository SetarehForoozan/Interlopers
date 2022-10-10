import numpy as np
from spl import *
from default_config import *
from read_data import *
import mono_functions as func
from fitting_params import *

def ratio_func(f):
    """ _summary_
    This is the prefactor in front of the cross-correlation that comes from auto-correlation
    Returns: prefactor
        _type_: float
    """
    f = f
    return (1+2*f**2-2*f)

def xi_poles_obs(theta_mono, theta_quad, r_data, mu_data, weight, xi_fid, xi_cross, case):
    if case != 0: 
        xi_cross_L0, xi_cross_L2 = xi_cross[:,0], xi_cross[:,1]
        b_a,alpha,fi = theta_mono[-3:]
    if case == 0:
        b_a,alpha = theta_mono[-2:]
    xi0_fid, xi2_fid, xi4_fid = xi_fid[:,0], xi_fid[:,1], xi_fid[:,2]
    epsilon = theta_quad[-1]
    nr_fid = len(r_data)
    r_fid2 = r_data**2
    opep = 1+epsilon
    nmuo2 = int(np.size(mu_data)/2)
    mu_fid2 = mu_data[0:nmuo2]**2.
    weighto2 = weight[0:nmuo2]
    mu_obs2 = 1./(1. + (1./mu_fid2-1)/opep**6.)
    L2_obs = 0.5*(3.*mu_obs2-1.)
    L4_obs = 0.125*(35.*mu_obs2**2. - 30.*mu_obs2 + 3.)
    r_obs = np.zeros((nr_fid,nmuo2))
    xi0_obs = np.zeros((nr_fid,nmuo2))
    xi2_obs = np.zeros((nr_fid,nmuo2))
    xi4_obs = np.zeros((nr_fid,nmuo2))
    for i in range(0,nmuo2):
        r_obs[:,i] = r_data*alpha*np.sqrt(opep**4.*mu_fid2[i] + (1.-mu_fid2[i])/opep**2.)
    r_obs2 = r_obs**2.
    for i in range(0,nmuo2):
        xi0_obs[:,i] = np.interp(r_obs[:,i], r_data, xi0_fid*r_data**2)/r_obs[:,i]**2 #spllin(r_obs[:,i],r_fid2*xi0_fid,r_data) / r_obs2[:,i]
        xi2_obs[:,i] = np.interp(r_obs[:,i], r_data, xi2_fid*r_data**2)/r_obs[:,i]**2 * L2_obs[i] #spllin(r_obs[:,i],r_fid2*xi2_fid,r_data) / r_obs2[:,i] * L2_obs[i]
        xi4_obs[:,i] = np.interp(r_obs[:,i], r_data, xi4_fid*r_data**2)/r_obs[:,i]**2 * L4_obs[i] #spllin(r_obs[:,i],r_fid2*xi4_fid,r_data) / r_obs2[:,i] * L4_obs[i]
    xi_obs = xi0_obs + xi2_obs + xi4_obs
    xi_obs0 = xi_obs*0.0
    xi_obs2 = xi_obs*0.0
    L2_fid = 0.5*(3.*mu_fid2-1.)
    for i in range(0,nmuo2):
        xi_obs0[:,i] = xi_obs[:,i]*weighto2[i]
        xi_obs2[:,i] = xi_obs[:,i]*L2_fid[i]*weighto2[i]
    xi0 = np.sum(xi_obs0,1)
    xi2 = 5.*np.sum(xi_obs2,1)

    poly_q = polynomial_q(theta_quad, r_data)
    poly_a = polynomial_a(theta_mono, r_data)

    xi0 = xi0 * b_a + poly_a
    xi2 = xi2 * b_a + poly_q
    if case != 0:
        xi_cross_L0 = xi_cross_L0 / ratio_func(fi)
        xi_cross_L2 = xi_cross_L2 / ratio_func(fi)
        xi0_eff = (1+2*fi**2-2*fi) * xi0 + 2*fi*(1-fi) * xi_cross_L0
        xi2_eff = (1+2*fi**2-2*fi) * xi2 + 2*fi*(1-fi) * xi_cross_L2
        return [xi0_eff,xi2_eff, xi_cross_L0,xi_cross_L2, xi0, xi2]
    else:
        return [xi0, xi2]

def polynomial_q(theta_quad, r):
    """_summary_
    Args:
        theta_quad (array): an array consisting of constants for quadrupole: 
            dof_q = 3: [q1,q2,q3]
            dof_q = 2: [q1,q2]
            dof_q = 1: [q1]
        r (array): the r array 

    Returns:
        array with len(r): the sum of all polynomials
    """
    vec_a = theta_quad[:dof_q]
    poly_q = 0
    for i in range(1,dof_q+1):
        poly_q +=  vec_a[i-1] * r**(-(3-i))
    return poly_q

def polynomial_a(theta_mono, r):
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
    vec_a = theta_mono[:dof_a]
    poly_a = 0
    for i in range(1,dof_a+1):
        poly_a +=  vec_a[i-1] * r**(-(3-i))
    return poly_a

def chi2_func(theta_mono, theta_quad, r_data, mu_data, weight, xi_fid, xi_cross, xi_obs, covmat, case):
    xi0, xi2 = xi_poles_obs(theta_mono, theta_quad, r_data, mu_data, weight, xi_fid, xi_cross, case)[:2]
    model = np.concatenate((xi0, xi2), axis = 0)
    xi0_obs, xi2_obs = xi_obs[:,0], xi_obs[:,1]
    data = np.concatenate((xi0_obs, xi2_obs), axis = 0)
    #print('BE CAREFUL MODEL SHAPE SHOULD BE 1,2n',model.shape)
    diff = data - model
    covinv = np.linalg.inv(covmat)
    chi2 = (diff.dot(covinv)).dot(diff)
    return chi2

def log_likelihood(theta_mono, theta_quad, r_data, mu_data, weight, xi_fid, xi_cross, xi_obs, covmat, case, posterior = 'Percival'):
    """_summary_
    Args:
        theta_mono (_type_): _description_
        theta_quad (_type_): _description_
        r_data (_type_): _description_
        mu_data (_type_): _description_
        weight (_type_): _description_
        xi_fid (_type_): _description_
        xi_cross (_type_): _description_
        xi_obs (_type_): _description_
        covmat (_type_): _description_
        case (_type_): can be 0, 8 
        posterior (str, optional): Can be between 'chi2' or 'Percival'. Defaults to 'Percival'.
    Returns:
        float: posterior
    """
    chi2 = chi2_func(theta_mono, theta_quad, r_data, mu_data, weight, xi_fid, xi_cross, xi_obs, covmat, case)
    if posterior == 'Percival': #### given by https://arxiv.org/abs/2108.10402
        ns = 1000
        nd = len(r_list) * 2 #because we have mono+quad
        nt = len(theta_mono)+len(theta_quad)
        B = (ns-nd-2)/(ns-nd-1)/(ns-nd-4)
        m = nt + 2 + (ns - 1 + B*(nd-nt))/(1+B*(nd-nt))
        post = (-m/2)*np.log((1+chi2/(ns-1)))
    elif posterior == 'chi2':
        post = - 0.5 * chi2
    else:
        print('Posterior is not defined. Choose between Percival and chi2')
        sys.exit()
    return post

def log_prior(theta_mono, theta_quad, case):
    if case != 0:
        alpha,fi = theta_mono[-2:]
        fi_min = 0
        fi_max = 0.5
    else:
        alpha = theta_mono[-1]
    epsilon = theta_quad[-1]
    if case != 0 :
        if 0.6 < alpha < 1.4 and fi_min < fi < fi_max and -1<epsilon<1:
            return 0.0
    else: 
        if 0.6 < alpha < 1.4 and -1<epsilon<1:
            return 0.0
    return -np.inf
                    
def log_probability(theta_all, r_data, mu_data, weight, xi_fid, xi_cross, xi_obs, covmat, case, posterior):
    if case != 0: uuu = 3 
    else: uuu = 2
    theta_mono = theta_all[:dof_a+uuu]
    theta_quad = theta_all[dof_a+uuu:]
    lp = log_prior(theta_mono, theta_quad, case)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta_mono, theta_quad, r_data, mu_data, weight, xi_fid, xi_cross, xi_obs, covmat, case, posterior)

def set_initialparams(dof_a, dof_q, case):
    """_summary_

    Args:
        dof_a (int): Can be 1, 2, or 3. degrees of freedom for polynomials of monopole
        dof_q (int): Can be 1, 2, or 3. degrees of freedom for polynomials of monopole
        case (int): Can be 0 (no correction) or 8 (correction)

    Returns:
        array: initial parameters for MCMC code
    """
    #theta_mono = [a1,a2,a3,b,alpha,fi]
    #theta_quad = [q1,q2,q3,epsilon]    
    if dof_a == 3:
        theta_mono_0 = [10,1e-2,1e-4,1,1,0.01] 
        params = [r"$a_1$",r"$a_2$",r"$a_3$",r"$b$",r"$\alpha$", r"$f_{i}$"]
    elif dof_a == 2:
        theta_mono_0 = [10,1e-2,1,1,0.01] 
        params = [r"$a_1$",r"$a_2$",r"$b$",r"$\alpha$", r"$f_{i}$"]
    elif dof_a == 1:
        theta_mono_0 = [10,1,1,0.01] 
        params = [r"$a_1$",r"$b$",r"$\alpha$", r"$f_{i}$"]
    
    if case == 0:
        theta_mono_0 = theta_mono_0[:-1]        
        params = params[:-1]
        
    if dof_q == 3:    
        theta_quad_0 = [-50,2,5e-3,0.0]
        params += [r"$q_1$",r"$q_2$",r"$q_3$",r"$\epsilon$"]
    elif dof_q == 2:
        theta_quad_0 = [-50,2,0.0]
        params += [r"$q_1$",r"$q_2$",r"$\epsilon$"]
    elif dof_q == 1:
        theta_quad_0 = [-50,0.0]
        params += [r"$q_1$",r"$\epsilon$"]
    return theta_mono_0, theta_quad_0, params
    
def make_path(path):
    isFile = os.path.isdir(path)
    if not(isFile):
        print('path \n"%s"\n did not exist. Now creating...'%(path))
        os.mkdir(path)
    return 0