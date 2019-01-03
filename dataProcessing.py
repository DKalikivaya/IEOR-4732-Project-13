# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 11:04:46 2018

@author: p12de
"""
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm

def readNPlot(excel_file, option):
        
    # read data into a data frame
    df = pd.read_excel(excel_file)
    # create the 'Mid' variable
    df['Mid'] = df[['Bid','Ask']].mean(axis=1)
     

    # define strikes and maturities
    #strikes = np.arange(170., 250. + 5.0, 5.0)
    strikes = np.sort(df.Strike.unique())
    maturities = np.sort(df.Maturity_days.unique())
    
    # define a grid for the surface
    X, Y = np.meshgrid(strikes, maturities)
    optionPrices = np.empty([len(maturities), len(strikes)])
    
    if(option == 'Call'):
        
        df_calls = df[df['Option_type'] == 'Call'][['Maturity_days', 'Strike', 'Mid']]
        #print(df_calls.head())

        # we use linear interpolation for missing strikes
        for i in range(len(maturities)):
            s = df_calls[df_calls.Maturity_days == maturities[i]]['Strike']
            price = df_calls[df_calls.Maturity_days == maturities[i]]['Mid']
            f = interpolate.interp1d(s, price, bounds_error=False, fill_value="extrapolate")
            optionPrices[i, :] = f(strikes) 

    elif(option == 'Put'):
        df_puts = df[df['Option_type'] == 'Put'][['Maturity_days', 'Strike', 'Mid']]
        #print(df_puts.head())

        for i in range(len(maturities)):
            s = df_puts[df_puts.Maturity_days == maturities[i]]['Strike']
            price = df_puts[df_puts.Maturity_days == maturities[i]]['Mid']
            f = interpolate.interp1d(s, price, bounds_error=False, fill_value="extrapolate")
            optionPrices[i, :] = f(strikes) 

    #plot the surface
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    #ax.plot_wireframe(X, Y, callPrices, rstride=1, cstride=1)
    ax.plot_surface(X, Y, optionPrices, cmap=cm.coolwarm)
    ax.set_ylabel('Maturity (days)') 
    ax.set_xlabel('Strike') 
    #plt.save
    #fig('appleCallSurface.png')
    plt.show()
    
    return maturities, strikes, optionPrices