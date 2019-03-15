# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:13:51 2019

@author: Rick
"""


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               Imports
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

from math import *
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from datareader import CL_CD_series1
from conversion_helpers import lbs_to_kg



#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               Weight calculations
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------


def weight(weight_zero, CL_CD_series1):
    # Input: initial total weight; elapsed time; fuel mass used so far
    
    weight = []
    weight.append(weight_zero)
    
    # Import fuel used and convert from lbs to kg
    f_used = lbs_to_kg*CL_CD_series1[:,-1]
    
    i = 0
    for i in range(len(f_used)):
        # Subtract the fuel mass that has been burned during the time interval
        burn = weight[i] - f_used[i]
        weight.append(burn)
        
    # Output: array with weight at each time interval
    return weight


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               CL calculations
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------


def CL(Vtas, rho, weight, S):
    # Input: true airspeed, air density, weight, wing surface area
    
    CL = []
    i = 0
    while i < len(weight):
        # Calculate CL for each time interval
        CL.append(weight[i]/(1./2*Vtas[i]**2*S*rho[i]))
        i+=1
        
    # Output: array with CL values at each time interval
    return CL


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               CD calculations
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------



def CD(CL, Tr, A, rho, Vtas, S):
    # Input: CL; thrust; aspect ratio; air density; true airspeed; 
    # wing surface area
    
    CD =[]
    x = []
    i = 0
    while i < len(rho):
        # Calculate CD for each time interval, from given values for thrust
        CD.append(2*Tr/(rho[i]*S*Vtas[i]**2))
        x.append(CL[i]**2)
        i+=1
        
    # Obtain slope of CD CL^2 diagram to find the Oswald factor
    slope = np.polyfit(x,CD,1,full=False)[0]
    oswald_factor = 1./(pi*A*slope)
    
    # Get CD_zero from the intersection with the y-axis
    CD_zero = np.polyfit(x,CD,1,full=False)[1]

    # Output: array with CD values at each time interval; CD_zero; oswald factor 
    return CD, CD_zero, oswald_factor


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               Plots
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------



def Plots(CL, CD, alpha):
    # Plot CL against angle of attack
    plt.figure()
    plt.plot(alpha,CL)
    plt.title('CL-alpha')
    plt.xlabel('Angle of attack [degrees]')
    plt.ylabel('CL [-]')
    plt.show()

    # Plot CD against angle of attack
    plt.figure()
    plt.plot(alpha,CD)
    plt.title('CD-alpha')
    plt.xlabel('Angle of attack [degrees]')
    plt.ylabel('CD [-]')
    plt.show()
    
    # Plot CD against CL
    plt.figure()
    plt.plot(CD,CL)
    plt.title('CD-CL')
    plt.xlabel('CL [-]')
    plt.ylabel('CD [-]')
    plt.show()
    return

