import numpy as np
import sys
#Arguments should be input as: 
#################  SIMS, CASE, A, FRAC, DOF_A, N  ###############

nsim  = int(sys.argv[1])
case  = int(sys.argv[2])
A     = float(sys.argv[3])
frac  = int(sys.argv[4])
dof_a = int(sys.argv[5])
dof_q = int(sys.argv[6])
sims  = np.arange(nsim) + int(sys.argv[7])
rsd   = str(sys.argv[8])

"""
#sims = np.arange(200) + 400
sims = np.arange(10)
nsim = 10
case = 8
A    = 85 #85,90,97.4 #displacement
frac = 10
dof_a = 3
dof_q = 3
rsd = "norsd"
"""
nsim_covmat = 1000

hhdir = './Interlopers/%s_%0.0f'%(rsd,int(A))
meanfile_name = 'mean%dto%d'%(np.min(sims), np.max(sims))
