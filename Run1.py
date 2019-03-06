# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:13:51 2019

@author: Rick
"""

from math import *
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


#df = pd.read_excel(r'C:/Users/Rick/Documents/Python/FD_assignment/REFERENCE_Post_Flight_Datasheet_Flight.xlsx')
#dat = pd.DataFrame(df, columns = ['Unnamed: 2','Unnamed: 3','Unnamed: 4','Unnamed: 5','Unnamed: 6','Unnamed: 7','Unnamed: 8','Unnamed: 9','Unnamed: 10','Unnamed: 11','Unnamed: 12'])

#xl_workbook = pd.ExcelFile(r'C:/Users/Rick/Documents/Python/FD_assignment/REFERENCE_Post_Flight_Datasheet_Flight.xlsx')  # Load the excel workbook
#df = xl_workbook.parse("Sheet1")  # Parse the sheet into a dataframe
#aList = df.columns('Unnamed: 3')

lbs_to_kg = 0.45359237
weight_zero = 60500/9.80665     # from N to kg

Time = [0,19*60+17, 21*60+37, 23*60+46, 26*60+4, 29*60+47, 32*60]
time = np.array(Time)/3600.
F_used = [0,360, 412, 447, 478, 532, 570]
f_used = np.array(F_used)*lbs_to_kg

def weight(weight_zero, time, f_used):
    # Input: initial total weight; elapsed time; fuel mass used so far
    
    weight = []
    weight.append(weight_zero)
    i = 1
    for i in range(1,len(time)):
        # Subtract the fuel mass that has been burned during the time interval
        burn = weight[i-1] - f_used[i]*(time[i]-time[i-1])
        weight.append(burn)     # Check this!!

    # Output: array with weight at each time interval
    return weight






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






def CD(CL, Tr, aspect_ratio, rho, Vtas, S):
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
    oswald_factor = 1./(pi*aspect_ratio*slope)
    
    # Get CD_zero from the intersection with the y-axis
    CD_zero = np.polyfit(x,CD,1,full=False)[1]

    # Output: array with CD values at each time interval; CD_zero; oswald factor 
    return CD, CD_zero, oswald_factor






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

