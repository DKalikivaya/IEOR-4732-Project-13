# -*- coding: utf-8 -*-
"""
Created on Tue Dec 8 10:47:35 2018

@author: p12de
"""


import numpy as np
import pandas as pd
import cmath
import math
from scipy.optimize import fmin
from scipy.optimize import brute
from scipy.optimize import fmin_bfgs
from calibration import eValue


def myRange(start, finish, increment):
    while (start <= finish):
        yield start
        start += increment
#Model Parameters Starting Points        
def generic_model(model):
    if (model == 'Heston'):
#         kappa = 2.3
#         theta = 0.046
#         sig = 0.0825
#         rho = -0.53
#         v0 = 0.054
        kappa = 1.3
        theta = 0.056
        sig = 0.0625
        rho = -0.43
        v0 = 0.064

        params = []
        params.append(kappa)
        params.append(theta)
        params.append(sig)
        params.append(rho)
        params.append(v0)
        print("Parameters of the model",params)
        
    elif (model == 'VGSA'):
        sig = 0.5
        nu = 0.13
        theta = -0.13
        kappa = 10
        eta = 2
        lda = 10
        
        params = []
        params.append(sig)
        params.append(nu)
        params.append(theta)
        params.append(kappa)
        params.append(eta)
        params.append(lda)
        print("Parameters of the model",params)
        
    elif (model == 'VGSSD'):
        sig = 0.25
        nu = 0.4
        theta = -0.3
        gamma = 0.05
        
        params = []
        params.append(sig)
        params.append(nu)
        params.append(theta)
        params.append(gamma)
        print("Parameters of the model",params)
        
    return params


#Optimal Parameters
def generic_algo(algo, model, params, arg,  marketPrices, maturities_years, strikes, r, q, S0, alpha, eta_global, n):
    
    #BFGS
    if (algo == 'BFGS'):
        num_iter=1
        [xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflg] = fmin_bfgs(eValue,params,args=arg, fprime=None, 
                                                                              callback=None,
                                                                          maxiter=10, full_output=True, retall=False)
        params2 = xopt
        print(params2)

       
    #Nelder Mead
    elif (algo == 'Nelder_Mead'):
        num_iter=1
        t = fmin(eValue, params, args=arg, xtol=1e-4, ftol=1e-4,maxiter=10,maxfun=400,callback=None,disp=True,retall=False,full_output=True)
        print("t:", t)
        params2 = t[0]
        

    elif (algo == 'Grid_Search'):
        rmseMin = 1.0e6
        num_iter = 1
        if (model == 'Heston'):
            for kappa in myRange(1.3,1.8,0.5):
                for theta in myRange(0.036,0.056,0.01):
                    for sig in myRange(0.0725,0.0925,0.1):
                        for rho in myRange(-0.63,-0.53,0.1):
                            for v0 in myRange(0.044,0.054,0.01):
                                params = []
                                params.append(kappa)
                                params.append(theta)
                                params.append(sig)
                                params.append(rho)
                                params.append(v0)
                                #print('i = ' + str(num_iter))
                                #print("params:",params)
                                num_iter += 1
                                rmse = eValue(params, marketPrices, maturities_years, strikes, r, q, S0, alpha, eta_global, n, model)
                                #print("rmse:", rmse)
                                if (rmse < rmseMin):
                                    rmseMin = rmse
                                    params2 = params
                                    print('\nnew min found')
                                    print(rmseMin)
                                    #print(params2)
                                    print('')
        
        elif (model == 'VGSA'):
            for sig in myRange(0.1,0.2,0.1):
                for nu in myRange(0.2,0.3,0.1):
                    for theta in myRange(-0.4,-0.1,0.1):
                        for kappa in myRange(6.0,7.0,1.0):
                            for eta in myRange(0.3,0.4,0.1):
                                for lda in myRange(0.3,0.4,0.1):
                                    params = []
                                    params.append(sig)
                                    params.append(nu)
                                    params.append(theta)
                                    params.append(kappa)
                                    params.append(eta)
                                    params.append(lda)
                                    num_iter += 1
                                    rmse = eValue(params, marketPrices, maturities_years, strikes, r, q, S0, alpha, eta_global, n, model)
                                    if (rmse < rmseMin):
                                        rmseMin = rmse
                                        params2 = params
                                        print('\nnew min found')
                                        print(rmseMin)
                                        #print(params2)
                                        print('')
        elif (model == 'VGSSD'):
            for sig in myRange(0.2,0.4,0.1):
                for nu in myRange(0.05,0.15,0.05):
                    for theta in myRange(-0.05,0.05,0.05):
                        for gamma in myRange(0.2,0.5,0.1):
                            params = []
                            params.append(sig)
                            params.append(nu)
                            params.append(theta)
                            params.append(gamma)
                            num_iter += 1
                            rmse = eValue(params, marketPrices, maturities_years, strikes, r, q, S0, alpha, eta_global, n, model)
                            if (rmse < rmseMin):
                                rmseMin = rmse
                                params2 = params
                                print('\nnew min found')
                                print(rmseMin)
                                #print(params2)
                                print('')

        print('\nSolution of grid search:')                        
        print(params2)
        print('Optimal rmse = ' + str(rmseMin))
    return params2
