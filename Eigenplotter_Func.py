# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:12:34 2019

@author: daanv
"""

from datareader import getFDValues,importExcelData,convertToSec,convertToTimeStr,fixList
import matplotlib.pyplot as plt
import numpy as np
import scipy.fftpack

def noiseFilter(y,dt=0.1,threshold=0.5, plotFourier=False):
    #needs fixed list as input!
    ylen=len(y)
    yf = scipy.fftpack.fft(y)    
    xf = np.linspace(0.0, 1.0/(2.0*dt), ylen/2)
    if plotFourier:
        plt.figure()
        #only using positive half of fourier
        plt.plot(xf, 2.0/ylen * np.abs(yf[:ylen/2]))
        plt.show()
    
    #filter out noise (defined of percentage of the threshold)
    absTreshold=threshold*max(np.abs(yf[1:ylen/2]))    
    for i in range(len(yf)):
        if np.abs(yf[i])<absTreshold:
            yf[i]=0
    return scipy.fftpack.ifft(yf)


def time_stamps(n, eigenmotions, lengths, time_list):
    #n is which eigenmode is to be plotted
    time_start = convertToSec(eigenmotions[n])
    time_end = time_start+lengths[n]
    
    indices = []
    for i in range(len(time_list)):
        if time_start<time_list[i]<time_end:
            indices.append(i)
    return indices




def GetLists():   
    keylist,desclist,unitlist,newDict = getFDValues('FlightData.mat')
    date_of_flight, flight_number, TO_time, LND_time, passengerMass, passengerNames\
    , passengerPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2, ACC_Trim,\
     El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, Cg_shift, eigenmotions\
     = importExcelData('Post_Flight_Datasheet_13_03_V2.csv')
    time_list = newDict.get("time") 
    velocity_list = newDict.get("Dadc1_tas")
    altitude_list = newDict.get("Dadc1_alt")
    alpha_list = newDict.get("vane_AOA")
    elevator_list = newDict.get("delta_e")
    roll_list = newDict.get("Ahrs1_Roll")
    yaw_list = newDict.get("Ahrs1_Pitch")
    aileron_list = newDict.get("delta_a")
    rudder_list = newDict.get("delta_r")
    return time_list,velocity_list,altitude_list,alpha_list,elevator_list,roll_list,yaw_list,aileron_list,rudder_list, eigenmotions


def getEigenmotions():
    lengths = [200.,60.,45.,30.,60.,200.]
    modes = ["phugoid", "short period", "dutch roll","dutch roll yd", "aperiodic roll", "spiral" ]
    time_list,velocity_list,altitude_list,alpha_list,elevator_list,roll_list,yaw_list,aileron_list,rudder_list, eigenmotions=GetLists()
    times=[]
    velocities=[]
    altitudes=[]
    alphas=[]
    elevators=[]
    rolls=[]
    yaws=[]
    ailerons=[]
    rudders=[]
    for n in range(len(modes)):
        indices=time_stamps(n, eigenmotions, lengths, time_list)
        times.append(time_list[indices[0]:indices[-1]])
        velocities.append(velocity_list[indices[0]:indices[-1]])
        altitudes.append(altitude_list[indices[0]:indices[-1]])
        alphas.append(alpha_list[indices[0]:indices[-1]])
        elevators.append(elevator_list[indices[0]:indices[-1]])
        rolls.append(roll_list[indices[0]:indices[-1]])
        yaws.append(yaw_list[indices[0]:indices[-1]])
        ailerons.append(aileron_list[indices[0]:indices[-1]])
        rudders.append(rudder_list[indices[0]:indices[-1]])
    return times, velocities, alphas, elevators, rolls, yaws, ailerons, rudders

######phugoid (n=0): plotting speed, altitude and angle of attack against time
######short period (n=1): plotting speed, altitude and angle of attack against time 
#####dutch roll (n=2): plot yaw angle, roll angle, altitude and true airspeed
#####dutch roll YD (n=3): plot yaw angle, roll angle, altitude and true airspeed
#####aperiodic roll (n=4): plot yaw angle, roll angle, altitude and true airspeed
#####spiral (n=5): plot yaw angle, roll angle, altitude and true airspeed    
    
