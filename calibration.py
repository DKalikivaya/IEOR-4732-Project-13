# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 10:50:40 2018

@author: p12de
"""

import numpy as np
import math
import cmath

# Periodic Linear Extension
def paramMapping(x, c, d):
    
    if ((x>=c) & (x<=d)):       
        y = x
    else:       
        range = d-c
        n = math.floor((x-c)/range)
        if (n%2 == 0):
            y = x - n*range;
        else:
            y = d + n*range - (x-c)            
    return y

#VG Characteristic Function
def logCharFuncVG(u, sig, nu, theta):
    return -np.log(1.0 - 1j*nu*theta*u + sig*sig*nu*u*u/2.0)/nu;

def charFuncCIR(u, t, y, kappa, eta, lda):
    gm=np.sqrt(kappa**2-2*(lda**2)*u*1j)
    bbn=2*(np.exp(gm*t)-1)*(1j*u)
    bbd=gm-kappa+np.exp(gm*t)*(gm+kappa)
    bb=bbn/bbd    
    aau=2*(gm)*np.exp((gm+kappa)*t/2)
    aad=gm-kappa+np.exp(gm*t)*(gm+kappa)
    aa=(kappa*eta*2/(lda**2))*np.log(aau/aad)
    z = np.exp(aa+bb*y)
    return z
 
def generic_CF(u, params, S0, r, q, T, model):
    
    if (model == 'GBM'):      
        sig = params[0]
        mu = np.log(S0) + (r-q-sig**2/2)*T
        a = sig*np.sqrt(T)
        phi = np.exp(1j*mu*u-(a*u)**2/2)  
        
    elif(model == 'Heston'):
        
        kappa  = params[0]
        theta  = params[1]
        sigma  = params[2]
        rho    = params[3]
        v0     = params[4]
        
        kappa = paramMapping(kappa,0.1, 20)
        theta = paramMapping(theta,0.001, 0.4)
        sigma = paramMapping(sigma,0.01, 0.6)
        rho   = paramMapping(rho  ,-1.0, 1.0)
        v0    = paramMapping(v0   ,0.005, 0.25)
        
        tmp = (kappa-1j*rho*sigma*u)
        g = np.sqrt((sigma**2)*(u**2+1j*u)+tmp**2)
        
        pow1 = 2*kappa*theta/(sigma**2)
        
        numer1 = (kappa*theta*T*tmp)/(sigma**2) + 1j*u*T*r + 1j*u*math.log(S0)
        log_denum1 = pow1 * np.log(np.cosh(g*T/2)+(tmp/g)*np.sinh(g*T/2))
        tmp2 = ((u*u+1j*u)*v0)/(g/np.tanh(g*T/2)+tmp)
        log_phi = numer1 - log_denum1 - tmp2
        phi = np.exp(log_phi)
        

    elif (model == 'VG'):
        
        sigma  = params[0];
        nu     = params[1];
        theta  = params[2];

        if (nu == 0):
            mu = math.log(S0) + (r-q - theta -0.5*sigma**2)*T;
            phi  = math.exp(1j*u*mu) * math.exp((1j*theta*u-0.5*sigma**2*u**2)*T);
        else:
            mu  = math.log(S0) + (r-q + math.log(1-theta*nu-0.5*sigma**2*nu)/nu)*T;
            phi = cmath.exp(1j*u*mu)*((1-1j*nu*theta*u+0.5*nu*sigma**2*u**2)**(-T/nu));
            
    elif(model =='VGSA'):
        sig   = params[0]
        nu    = params[1]
        theta = params[2]
        kappa = params[3]
        eta   = params[4]
        lda   = params[5]

        tmp =  1j*(np.log(S0)+(r-q)*T)*u
        u1 = -1j*logCharFuncVG(u, sig, nu, theta)
        u2 = -1j*logCharFuncVG(-1j, sig, nu, theta)
        numer = charFuncCIR(u1, T, 1.0/nu, kappa, eta, lda)
        denom = charFuncCIR(u2, T, 1.0/nu, kappa, eta, lda)
        phi = np.exp(tmp)*numer/np.exp(denom, 1j*u)
            
    elif(model == 'VGSSD'):
        sig   = params[0]
        nu    = params[1]
        theta = params[2]
        gamma = params[3]
        
        first_term = np.exp(1j*u*(np.log(S0)+(r-q)*T))
        si_u = pow((1.0 - 1j*u*pow(T,gamma)*theta*nu + pow(sig, 2.0)*pow(u*pow(T,gamma), 2.0)*nu / 2.0),-1/nu)
        si_complx = pow((1.0 - 1j*(-1j)*pow(T,gamma)*theta*nu + pow(sig, 2.0)*pow((-1j)*pow(T,gamma), 2.0)*nu / 2.0),-1/nu)
        phi = first_term*si_u/si_complx
        
    return phi

#SAME module from mudules for Calb
def generic_FFT(params, S0, K, r, q, T, alpha, eta_global, n, model):
    
    N = 2**n
    
    # step-size in log strike space
    lda = (2*np.pi/N)/eta_global
    
    #Choice of beta
    #beta = np.log(S0)-N*lda/2
    beta = np.log(K)
    
    # forming vector x and strikes km for m=1,...,N
    km = np.zeros((N))
    xX = np.zeros((N))
    
    # discount factor
    df = math.exp(-r*T)
    
    nuJ = np.arange(N)*eta_global
    psi_nuJ = generic_CF(nuJ-(alpha+1)*1j, params, S0, r, q, T, model)/((alpha + 1j*nuJ)*(alpha+1+1j*nuJ))
    
    for j in range(N):  
        km[j] = beta+j*lda
        if j == 0:
            wJ = (eta_global/2)
        else:
            wJ = eta_global
        xX[j] = cmath.exp(-1j*beta*nuJ[j])*df*psi_nuJ[j]*wJ
     
    yY = np.fft.fft(xX)
    cT_km = np.zeros((N))  
    for i in range(N):
        multiplier = math.exp(-alpha*km[i])/math.pi
        cT_km[i] = multiplier*np.real(yY[i])
    
    return km, cT_km

#def eValue(params, marketPrices, maturities, strikes, r, q, S0, alpha, eta_global, n, model):

def eValue(params, *args):
    marketPrices = args[0]
    maturities = args[1]
    strikes = args[2]
    r = args[3]
    q = args[4]
    S0 = args[5]
    alpha = args[6]
    eta_global = args[7]
    n = args[8]
    model = args[9]

    lenT = len(maturities)
    lenK = len(strikes)
    
    modelPrices = np.zeros((lenT, lenK))
    #print(marketPrices.shape)

    count = 0
    mae = 0
    for i in range(lenT):
        for j in range(lenK):
            count  = count+1
            T = maturities[i]
            K = strikes[j]
            [km, cT_km] = generic_FFT(params, S0, K, r, q, T, alpha, eta_global, n, model)
            modelPrices[i,j] = cT_km[0]
            tmp = marketPrices[i,j]-modelPrices[i,j]
            mae += tmp**2
    
    rmse = math.sqrt(mae/count)
    return rmse