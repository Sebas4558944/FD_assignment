# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:13:51 2019

@author: Rick
"""

from math import *
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from datareader import CL_CD_series1
from conversion_helpers import lbs_to_kg

#df = pd.read_excel(r'C:/Users/Rick/Documents/Python/FD_assignment/REFERENCE_Post_Flight_Datasheet_Flight.xlsx')
#dat = pd.DataFrame(df, columns = ['Unnamed: 2','Unnamed: 3','Unnamed: 4','Unnamed: 5','Unnamed: 6','Unnamed: 7','Unnamed: 8','Unnamed: 9','Unnamed: 10','Unnamed: 11','Unnamed: 12'])

#xl_workbook = pd.ExcelFile(r'C:/Users/Rick/Documents/Python/FD_assignment/REFERENCE_Post_Flight_Datasheet_Flight.xlsx')  # Load the excel workbook
#df = xl_workbook.parse("Sheet1")  # Parse the sheet into a dataframe
#aList = df.columns('Unnamed: 3')

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               Data
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------


weight_zero = 60500/9.80665     # from N to kg

#Time = [0,19*60+17, 21*60+37, 23*60+46, 26*60+4, 29*60+47, 32*60]
#time = np.array(Time)/3600.
#F_used = [0,360, 412, 447, 478, 532, 570]
#f_used = np.array(F_used)*lbs_to_kg


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               Weight calculations
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#print 2*CL_CD_series1[:,-1]
f_used = list(CL_CD_series1[:,-1])
for i in f_used:
    f_sed[i] = float(f_used[i])
print f_sed

def weight(weight_zero, CL_CD_series1):
    # Input: initial total weight; elapsed time; fuel mass used so far
    
    weight = []
    weight.append(weight_zero)
    
    print f_used
    
    i = 1
    for i in range(len(f_used)):
        # Subtract the fuel mass that has been burned during the time interval
        burn = weight[i-1] - f_used[i]
        weight.append(burn)     # Check this!!
        
    # Output: array with weight at each time interval
    return weight

answer = weight(15000, CL_CD_series1)
print answer
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

