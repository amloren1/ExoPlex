import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex as exo
import pdb

from params import *

import grid_search as gs

file_name = '1_ME_TESTER.dat'

########################################################################
'''
-Import data grids: CMF, Fe/Mg_bulk, Si/Mg_bulk, Radius

-find row in grid given the input parameters
'''
########################################################################

def data_grids(filename, cmf, femg_m, simg_m):

    dat = np.genfromtxt(filename, comments = '#')
    N_rows = int(len(cmf)*len(femg_m)*len(simg_m))

    R_p       = dat[:,3]
    FeMg_b    = dat[:,4]
    SiMg_b    = dat[:,5]

    #pdb.set_trace()
    return(R_p, FeMg_b, SiMg_b)

def search_grid(loc, CMF, FeMg_m, SiMg_m):

    cmf = loc[0]
    femg_m = loc[1]
    simg_m = loc[2]

    c = gs.bissect(CMF, cmf)[0]
    f = gs.bissect(FeMg_m, femg_m)[0]
    s = gs.bissect(SiMg_m, simg_m)[0]


    point = c*(len(FeMg_m)*len(SiMg_m))+f*len(SiMg_m)+s

    return(point)


#DEBUG
#print R[point]
#print FeMg_b[point]
#print SiMg_b[point]
#print point+2
#pdb.set_trace()


########################################################################
'''
Priors: (Fe/Mg)_mantle, (Si/Mg)_mantle, CMF (core mass fraction)

-because these parameters do not depend on one another, I can use uniform priors for
    all of the parameters
'''
########################################################################

def priors_uniform(FeMg_m, SiMg_m, CMF):

    n_femg = len(FeMg_m)
    n_simg = len(SiMg_m)
    n_cmf  = len(CMF)

    p_femg = 1/n_femg
    p_simg = 1/n_simg
    p_cmf  = 1/n_cmf

    joint_prior = p_femg*p_simg*p_cmf

    return(joint_prior)



########################################################################
'''
Likelihood function for data values:
R, (Fe/Mg)_P, (Si/Mg)_P

-each of the likelihoods will be a function of error between the model and data
    so this one function takes: d_i, sig_i, f_i (data, std_dev, model for parameter i)
'''
########################################################################


def likelihood_normal(d, sig, f):
    #variance
    var  = np.power(sig, 2)

    #model and data difference
    diff = d-f
    #coeficient to exponential in normal PDF equation
    coef = 1./np.sqrt(2*np.pi*var)


    arg_exp = (1./(2*var))*np.power(diff,2)
    exp  = np.exp(-arg_exp)

    return(exp*coef)


########################################################################
'''
Sampler
'''
########################################################################

def component_MH(n_samp, CMF, FeMg_m, SiMg_m, EP):
    #R, Fe/Mg_bulk, Si/Mg_bulk
    d   = [1.0, 0.8128, 0.9549 ]
    sig = [0.01, 0.1*0.8128, 0.1*0.9549]

    n_cmf  = len(CMF)
    n_femg = len(FeMg_m)
    n_simg = len(SiMg_m)

    cmf  = np.random.randint(0,n_cmf-1)
    femg = np.random.randint(0,n_femg-1)
    simg = np.random.randint(0,n_simg-1)

    params = [CMF[cmf], FeMg_m[femg], SiMg_m[simg]]
    params_test = [CMF[cmf], FeMg_m[femg], SiMg_m[simg]]

    params_test = [params[0], params[1], params[2]]
    #params_grid = np.transpose([CMF, FeMg_m, SiMg_m], axis = 0)
    params_grid = [CMF, FeMg_m, SiMg_m]

    #pdb.set_trace()
    for i in range(n_samp):
        for k in range(0,3):
            test_val = np.random.randint(params[k]-2, params[k]+2)

            if test_val < 0 or test_val > len(params_grid[k])-1:

                alpha = 0
                continue
            else:
                params_test[k] = params_grid[k][test_val]
                print test_val
                print params_test
                print params
                print 'k = {}'.format(k)
                print 'i = {}'.format(i)
                raw_input()

                point_test = search_grid(params_test, CMF, FeMg_m, SiMg_m)
                point_cur  = search_grid(params, CMF, FeMg_m, SiMg_m)
                #R, Fe/Mg_bulk, Si/Mg_bulk
                f_test = [EP[0][point_test], EP[1][point_test], EP[2][point_test]]
                f_cur  = [EP[0][point_cur], EP[1][point_cur], EP[2][point_cur]]

                prior_test    = priors_uniform(params_test[1], params_test[2], params_test[0])

                prior_current = priors_uniform(params[1], params[2], params[0])

                likely_test = likelihood_normal(d[0], sig[0], f_test[0])* \
                                likelihood_normal(d[1], sig[1], f_test[1])*\
                                  likelihood_normal(d[2], sig[2], f_test[2])


                likely_current = likelihood_normal(d[0], sig[0], f_cur[0])* \
                                    likelihood_normal(d[1], sig[1], f_cur[1])*\
                                      likelihood_normal(d[2], sig[2], f_cur[2])

                posterior_test = prior_test*likely_test
                posterior_current = prior_current*likely_current

                alpha = min([1, posterior_test/posterior_current])

                u = np.random.uniform(0,1)

                pdb.set_trace()
                if alpha > u:
                    params[k] = params_test[k]
                else:
                    continue

########################################################################
'''
Call sampler, test samples and accept/reject
'''
########################################################################













########################################################################
'''
Plot results on corner plot
'''
########################################################################



########################################################################
########################################################################
'''
RUN CODE HERE
'''

########################################################################
########################################################################

FeMg = np.arange(0.0,0.5,0.1)
SiMg = np.arange(0.1,0.5,0.1)
CMF  = np.arange(0.1,0.4,0.1)




R_p, FeMg_b, SiMg_b = data_grids(file_name, CMF, FeMg, SiMg)

EP = [R_p, FeMg_b, SiMg_b]

component_MH(1000, CMF, FeMg, SiMg, EP)

########################################################################








