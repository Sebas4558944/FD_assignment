# -*- coding: utf-8 -*-
"""
Created on Thu May 04 09:31:53 2017

@author: daanv
"""
from math import e, log, exp

def calcISA(hft=0,hm=0,P=0):
    #ALL UNITS IN SI STANDARD
    #Inputs:
    #Modes of calculation:
    #1-geopotential altitude in feet as input
    #2-geopotential altitude in m as input
    #3-pressure as input
    #Returns all properties of the local atmosphere

    #Only input one!!
    #P in Pa
    #hm is geopotential altitude in m
    #hft is geopotential altitude in ft
    #constants

    g=9.80665
    R=287.0
    #Altitude and temp gradients
    hl = [0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0, 100000.0]
    al = [-6.5e-3, 0.0e-3, 1.0e-3, 2.8e-3, 0.0e-3, -2.8e-3, -2.0e-3, 0.0]
    #Sea level values
    Tl = [288.15]
    Pl = [101325.0]
    rhol = [1.225]

    #Initialize lists
    i=0
    for i in range(len(al)):
        Ti = Tl[i] + al[i]*(hl[i+1]-hl[i])
        Tl.append(Ti)
        if al[i] == 0:
            Pi = Pl[i]* (e**-((g/(R*Tl[i]))*(hl[i+1]-hl[i])))
            rhoi = rhol[i]* (e**-((g/(R*Tl[i]))*(hl[i+1]-hl[i])))
        else:        
            Pi = Pl[i] * ((Tl[i+1]/Tl[i]) ** (-g / (R * al[i])))
            rhoi = rhol[i] * (Tl[i+1]/Tl[i]) ** ((-g / (R * al[i])) - 1)
        Pl.append(Pi)
        rhol.append(rhoi)   

    #If pressure is the input, calculate height, density and temperature
    if P!=0:
        i = 0
        while P<Pl[i+1]:
            i = i+1
        if al[i] == 0:
            T = Tl[i]
            h1 = hl[i] - log(P/Pl[i])*(R*T/g)
            rho = rhol[i] * ( e **( (g/(R*T))*(h1 - hl[i]))) 
        else: 
            T = Tl[i] *((P / Pl[i]) ** (-(al[i] * R)/ g))
            h1 = hl[i] +((T - Tl[i]) / al[i])
            rho = rhol[i] *((T / Tl[i])**((-g /(R * al[i]))-1))

    #If height is the input, calculate the pressure, density and temperature.
    elif hm!=0 or hft!=0:
        if hm!=0:
            h1=hm
        elif hft!=0:
            h1=hft*0.3048 
        i = 0
        while h1>hl[i+1]:
            i = i+1
        T= Tl[i] + al[i] * (h1-hl[i])
        if al[i] == 0:
            rho = rhol[i] * exp( -(g / (R * T) * (h1 - hl[i])))
            P = Pl[i]*exp(-(g/(R*T)*(h1 - hl[i])))
        else: 
            P = Pl[i]*((T / Tl[i])**(-g/(R*al[i])))
            rho = rhol[i]*((T / Tl[i])**((-g/(R*al[i]))-1))
    
    return P,T,rho,h1 #Pressure [Pa], T [K], Rho [kg/m^3], Geometric alt [m]
            
        
        
        
    
    
