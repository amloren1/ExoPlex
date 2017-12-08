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
import scipy.stats
import grid_search as gs

file_name = 'tester.dat'
plt.rc('font', family='serif')
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

    #pdb.set_trace()
    n_femg = len(FeMg_m)
    n_simg = len(SiMg_m)
    n_cmf  = len(CMF)

    p_femg = 1./n_femg
    p_simg = 1./n_simg
    p_cmf  = 1./n_cmf

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

    accept = 0

    n_cmf  = len(CMF)
    n_femg = len(FeMg_m)
    n_simg = len(SiMg_m)

    cmf_dist = np.empty(0)
    femg_dist = np.empty(0)
    simg_dist = np.empty(0)
    R_dist    = np.empty(0)
    femg_b_dist = np.empty(0)
    simg_b_dist = np.empty(0)

    cmf  = np.random.randint(0,n_cmf-1)
    femg = np.random.randint(0,n_femg-1)
    simg = np.random.randint(0,n_simg-1)

    params = [CMF[cmf], FeMg_m[femg], SiMg_m[simg]]
    params_test = [CMF[cmf], FeMg_m[femg], SiMg_m[simg]]

    params_test = [params[0], params[1], params[2]]
    #params_grid = np.transpose([CMF, FeMg_m, SiMg_m], axis = 0)
    params_grid = [CMF, FeMg_m, SiMg_m]

    sigma_p = 1.2
    #pdb.set_trace()
    for i in range(n_samp):
        for k in range(0,3):
            print i

            test_val = int(np.random.normal(params[k]*10, sigma_p))/10.
            #print 'test_val = {}'.format(test_val)
            #print 'curr_val = {}'.format(params[k])

            #raw_input()
            #continue
            if test_val < 0 or test_val > len(params_grid[k])-1:

                alpha = 0
                continue
            else:
                params_test[k] = test_val
                #print test_val
                #print params_test
                #print params
                #print 'k = {}'.format(k)
                #print 'i = {}'.format(i)
                #raw_input()

                q_test = scipy.stats.norm(params[k], sigma_p).pdf(test_val)

                q_cur =  scipy.stats.norm(test_val, sigma_p).pdf(params[k])

                point_test = search_grid(params_test, CMF, FeMg_m, SiMg_m)
                point_cur  = search_grid(params, CMF, FeMg_m, SiMg_m)

                #EP = R, Fe/Mg_bulk, Si/Mg_bulk
                f_test = [EP[0][point_test], EP[1][point_test], EP[2][point_test]]
                f_cur  = [EP[0][point_cur], EP[1][point_cur], EP[2][point_cur]]

                prior_test    = priors_uniform(params_grid[1], params_grid[2], params_grid[0])
                prior_current = priors_uniform(params_grid[1], params_grid[2], params_grid[0])

                likely_test = likelihood_normal(d[0], sig[0], f_test[0])* \
                                likelihood_normal(d[1], sig[1], f_test[1])*\
                                  likelihood_normal(d[2], sig[2], f_test[2])

                likely_current = likelihood_normal(d[0], sig[0], f_cur[0])* \
                                    likelihood_normal(d[1], sig[1], f_cur[1])*\
                                      likelihood_normal(d[2], sig[2], f_cur[2])

                posterior_test = prior_test*likely_test*q_cur
                posterior_current = prior_current*likely_current*q_test

                alpha = min([1, np.exp(np.log(posterior_test)-np.log(posterior_current))])

                u = np.random.uniform(0,1)

                #pdb.set_trace()
                if alpha > u:
                    accept +=1
                    params[k] = params_test[k]

                    if accept > 0:
                        print 'here'
                        cmf_dist = np.append(cmf_dist, params[0])
                        femg_dist = np.append(femg_dist, params[1])
                        simg_dist = np.append(simg_dist, params[2])
                        R_dist    = np.append(R_dist, f_cur[0])
                        femg_b_dist = np.append(femg_b_dist, f_cur[1])
                        simg_b_dist = np.append(simg_b_dist, f_cur[2])
                        print 'i = {}'.format(i)
                    continue
                else:
                    #print 'rejected!'
                    params_test[k] = params[k]
                    continue

    fs =22
    ls =16

    fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize = (10,8))

    ax1.hist(cmf_dist,  bins = len(CMF))
    ax2.hist(femg_dist, bins = len(FeMg_m))
    ax3.hist(simg_dist, bins = len(SiMg_m))

    ax1.set_ylabel('count', fontsize = fs)

    ax1.set_xlabel('CMF', fontsize = fs)
    ax2.set_xlabel('Fe/Mg_m', fontsize = fs)
    ax3.set_xlabel('Si/Mg_m', fontsize = fs)

    ax1.tick_params(axis='both', size=1, labelsize=ls)

    plt.tight_layout()


    fig, (ax4, ax5, ax6) = plt.subplots(3,1 ,figsize = (12,8))

    ax4.hist(R_dist)
    ax5.hist(femg_b_dist,range = (min(femg_b_dist), max(femg_b_dist)))
    ax6.hist(simg_b_dist, bins = len(SiMg_m))

    ax6.set_ylabel('count', fontsize = fs)

    ax4.set_xlabel('R', fontsize = fs)
    ax5.set_xlabel('Fe/Mg_bulk', fontsize = fs)
    ax6.set_xlabel('Si/Mg_bulk', fontsize = fs)

    ax4.tick_params(axis='both', size=1, labelsize=ls)
    ax5.tick_params(axis='both', size=1, labelsize=ls)
    ax6.tick_params(axis='both', size=1, labelsize=ls)

    plt.tight_layout()
    plt.show()

    data = np.stack((T_chain, v_chain), axis = -1)
#figure = corner.corner(data)




    pdb.set_trace()

########################################################################
'''
Call sampler, test samples and accept/reject
'''
########################################################################



def cornering(data):

    figure = corner.corner(data, labels = ['T (mK)', r'$\nu_0$'], show_titles = True,\
                       title_kwargs = {'fontsize': 18}, label_kwargs = {'fontsize': 18}, plot_contours = True)










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

FeMg = np.arange(0, 2.1,0.1)
SiMg = np.arange(0.1,2.1,0.1)
CMF  = np.arange(0.1,0.9,0.1)



R_p, FeMg_b, SiMg_b = data_grids(file_name, CMF, FeMg, SiMg)

EP = [R_p, FeMg_b, SiMg_b]

component_MH(50000, CMF, FeMg, SiMg, EP)

########################################################################








