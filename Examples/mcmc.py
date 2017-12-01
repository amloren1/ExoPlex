import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex as exo


from params import *




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

    return(p_femg, p_simg, p_cmf)



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














