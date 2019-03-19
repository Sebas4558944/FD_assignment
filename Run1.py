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
from datareader import *
from Cit_par import S, g
from conversion_helpers import lbs_to_kg , celsius_to_kelvin ,kts_to_ms
from reduced_condition_calculator import Conditions

CL_CD_series1 = importExcelData('Post_Flight_Datasheet_13_03_V2.csv')[9]
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               Weight calculations
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
weight_zero = (9165. + 2800. + 89. + 82. + 70. + 62. + 74. + 65. + 80. + 82. + 80.)

def weight(weight_zero, CL_CD_series1):
    # Input: initial total weight; elapsed time; fuel mass used so far
    
    weight = []
    weight.append(weight_zero)
    # Import fuel used and convert from lbs to kg
    f_used = lbs_to_kg*CL_CD_series1[:,-1]
    
    i = 0
    for i in range(len(f_used)):
        # Subtract the fuel mass that has been burned during the time interval
        burn = weight_zero - f_used[i]
        weight.append(burn)
        
    # Output: array with weight at each time interval in kg
    return weight


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               CL calculations
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------


def CL(CL_CD_series1, weight, S):
    # Input: datareader file [no adjustments required], mass [kg], wing surface area [m^2]
    V_c = CL_CD_series1[:,3]*kts_to_ms
    T_mCelsius = CL_CD_series1[:,-1]
    T_m = T_mCelsius + celsius_to_kelvin
    h = CL_CD_series1[:,2]*ft_to_m
    weight = weight(weight_zero, CL_CD_series1)
    CL = []
    i = 0
    while i < len(h):
        # Find the true airspeed [m/s]
        reduced_calculator = Conditions(h)
        rho = reduced_calculator.calc_density()
        mach = reduced_calculator.calc_mach(V_c[i])
        temp = reduced_calculator.calc_temperature(T_m[i],mach)
        Vtas = reduced_calculator.calc_V_t(temp-celsius_to_kelvin,mach)
        # Calculate CL for each time interval (and convert mass [kg] to weight [N])
        CL.append(weight[i]*g/(1./2*Vtas[i]**2*S*rho[i]))
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
    reduced_calculator = Conditions(h)
    rho = reduced_calculator.calc_density()
    mach = reduced_calculator.calc_m
    v_tas = reduced_calculator.calc_V_t()    
    
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

alpha = CL_CD_series1[:,4]

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

