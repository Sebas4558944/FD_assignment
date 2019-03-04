# -*- coding: utf-8 -*-
"""
Created on Thu May 04 09:31:53 2017

@author: daanv
"""
from math import *
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
#Set up lists
i=0
while i<len(al)-1:
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
    i=i+1

while True:
    print"ISA Calculations"
    print"What do you want to do?"
    print""
    print"   1) enter an altitude in feet"
    print"   2) enter an altitude in meters"
    print"   3) enter a pressure"
    print"   4) convert geopotential to geometric altitude" 
    print"   5) convert geometric to geopotential altitude" 
    print"   9) quit"
    choice = input("your choice:")
    
    if choice == 1 or choice == 2: 
        if choice == 1:
            hfeet = input("Geopotential altitude h in feet = ")
            h1 = hfeet*0.3048
        else: 
            h1 = input("Geopotential altitude h in meters = ")
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
    if choice == 3:
        P = input("Pressure in Pa:") 
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
    if choice == 4:
        E = 6371e3
        H = input("Input geopotential altitude in m:")
        Z = (E*H)/(E-H)
        print"geometric altitude =", round(Z,3), "m"
        print""
        dummy_input = raw_input("Press enter to continue")
    if choice == 5:
        E = 6371e3
        Z = input("Input geometric altitude in m:")
        H = (E*Z)/(E+Z)
        print"geopotential altitude =", round(H,3), "m"
        print""
        dummy_input = raw_input("Press enter to continue")
    if choice == 9:
        break
    if choice in range(1,4):
        Prel=(P/Pl[0])*100
        rhorel=(rho/rhol[0])*100
        hfeet=h1/0.3048
        FL=hfeet/100
        Tc=T-273.15
        print""
        print""
        print"*** International Standard Atmosphere calculations ***"
        print""
        print"Geopotential altitude h in feet (1ft = 0.3048m) =", round(hfeet,2), "(", round(h1,2),"meters)"
        print"FL", int(round(FL, 0))
        print""
        print"   T =", round(T, 3), "K (", round(Tc, 3), "C)"
        print"   p =", P, "Pa"
        print"   rho =", rho, "kg/m3"
        print""
        print""
        print"Relative to Sea level:"
        print"   p=", Prel, "%"
        print"   rho=", rhorel, "%" 
        dummy_input = raw_input("Press enter to continue")
    
    
    
    
    
