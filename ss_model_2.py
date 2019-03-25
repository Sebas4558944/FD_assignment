# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 14:59:51 2019

@author: msjor
"""
import matplotlib.pyplot as plt
import control.matlab as co
import numpy as np
from Eigenplotter_Func import getEigenmotions

def getStateSpace(alpha0,V0,th0):
# combining state vectors gives:
#x = [u,alhpa,theta,q,h,beta,psi,p,r,phi]
#u = [de,da,dr]
    from Cit_par import muc,c,Cmadot,KY2,CXu,CXa,CZa,CX0,CZq,Cmu,Cma,KX2,Cmq,mub,\
    CYr,KXZ,b,Clr,Cnr,Clp,Cnp,CZadot,CZ0,CXq,CZu,CXde,CZde,Cmde,CYbdot,Cnbdot,KZ2,\
    CYb,CL,CYp,Clb,Cnb,CYda,CYdr,Clda,Cldr,Cnda,Cndr
    C1S=np.matrix([[-2.*muc*c/(V0**2), 0., 0., 0., 0., 0. ,0., 0., 0., 0.],\
                   [0., (CZadot - 2. * muc) * c/V0, 0., 0. ,0., 0., 0., 0., 0., 0.],\
                   [0., 0., -c/V0, 0.,  0., 0. ,0., 0., 0., 0.],\
                   [0., Cmadot * c/V0, 0., -2.*muc*KY2*(c/V0)**2., 0., 0., 0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  1., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.]])
    
    C2S=np.matrix([[CXu/V0, CXa, CZ0 ,CXq*c/V0,  0., 0. ,0., 0., 0., 0.],\
                   [CZu/V0, CZa, -CX0, (CZq+2*muc) * (c/V0),  0., 0. ,0., 0., 0., 0.],\
                   [0. , 0. , 0.,  (c/V0),  0., 0. ,0., 0., 0., 0.],\
                   [Cmu/V0, Cma, 0, Cmq*c/V0,  0., 0. ,0., 0., 0., 0.],\
                   [alpha0-th0, V0*(1.+alpha0*th0), -V0*(1.+alpha0*th0),0.,  0., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.],\
                   [ 0. ,0., 0., 0.,  0., 0. ,0., 0., 0., 0.]])
    
    C1A=np.matrix([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., (CYbdot-2.*mub)*b/V0, 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., -0.5*b/V0, 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., 0., -2.*mub*KX2*(b/V0)**2., 2.*mub*KXZ*(b/V0)**2.,0.],\
                   [0., 0., 0., 0., 0., Cnbdot*b/V0, 0., 2.*mub*KXZ*(b/V0)**2., -2.*mub*KZ2*(b/V0)**2.,0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]])
    
    C2A=np.matrix([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],\
                   [0., 0., 0., 0., 0., CYb, CL, CYp*b/(2.*V0), (CYr-4.*mub)*b/(2.*V0),0.],\
                   [0., 0., 0., 0., 0., 0., 0., b/(2.*V0), 0., 0.],\
                   [0., 0., 0., 0., 0., Clb, 0., Clp*b/(2.*V0), Clr*b/(2.*V0), 0.],\
                   [0., 0., 0., 0., 0., Cnb, 0., Cnp*b/(2.*V0), Cnr*b/(2.*V0), 0.],\
                   [0., 0., 0., 0., 0., 0., 0., 0., 1., 0.]])
    
    ###COMBINING MATRICES
    C3=np.matrix([[CXde, 0., 0.],\
                   [CZde, 0., 0.],\
                   [0., 0., 0.],\
                   [Cmde, 0., 0.],\
                   [0., 0., 0.],\
                   [0., CYda, CYdr],\
                   [0., 0., 0.],\
                   [0., Clda, Cldr],\
                   [0., Cnda, Cndr],\
                   [0., 0., 0.]])
    C1=C1S+C1A     #np.add(C1S,C1A)
    C2=C2S+C2A     #np.add(C2S,C2A)
    
    A = np.linalg.inv(C1)*-C2
    B = np.linalg.inv(C1)*-C3
    #Outputs
    #desired outputs:
    #y = [u,h,theta,psi,phi]  (all of which are states)
    #
    CS = np.matrix([[1,0,0,0,0,0,0,0,0,0],\
                    [0,0,0,0,1,0,0,0,0,0],\
                    [0,0,1,0,0,0,0,0,0,0],\
                    [0,0,0,0,0,0,1,0,0,0],\
                    [0,0,0,0,0,0,0,0,0,1]])
    DS = np.zeros((5,3))
    #
    ##init ss
    SS=co.ss(A,B,CS,DS)
    
    return SS

#######plotting vars
label_font = 20
title_font = 25
reallabel='real'
fakelabel='fake'
modes = ["phugoid", "short period", "dutch roll","dutch roll yd", "aperiodic roll", "spiral" ]

plotting=[True, True, False, False, False, False]

indices, times, altitudes, velocities, alphas,  pitches, rolls, yaws, ailerons, rudders, elevators = getEigenmotions()
######phugoid (n=0): plotting speed, altitude and angle of attack against time          =
######short period (n=1): plotting speed, altitude and angle of attack against time     =
#####dutch roll (n=2): plot yaw angle, roll angle, altitude and true airspeed           =
#####dutch roll YD (n=3): plot yaw angle, roll angle, altitude and true airspeed        =
#####aperiodic roll (n=4): plot yaw angle, roll angle, altitude and true airspeed       =
#####spiral (n=5): plot yaw angle, roll angle, altitude and true airspeed               =


#======================================================================
#======                     phugoid                             =======
#======================================================================
#Phugoid is a pulse on the elevator
for n in range(6):
    if plotting[n]:
        #Get Flight Data
        ind=indices[n]
        time=times[n]
        h=altitudes[n]
        V=velocities[n]
        A=alphas[n]
        th=pitches[n]
        roll=rolls[n]
        yaw=yaws[n]
        dA=ailerons[n]
        dR=rudders[n]
        dE=elevators[n]
        
        #set init values
        V0=V[0]
        alpha0=A[0]
        th0=th[0]
        hp0=h[0]
        
        SS=getStateSpace(alpha0,V0,th0)
        if n==0 or n==1: #Sym
            #x = [u,alhpa,theta,q,h,beta,psi,p,r,phi]
            #u = [de,da,dr]
            #Rudders and ailerons set to zero since this shit is symmetrical
            u_input=[]
            for i in range(len(dE)):
                u_input.append([-dE[i],0.,0.])
                
            response, T, state = co.lsim(SS, T = time,U = u_input)
            
            #plotting u,h,theta,psi(roll),phi(yaw)
            speed_out = response[:,0]
            h_out = response[:,1]
            theta_out = response[:,2]
            psi_out = response[:,3]
            phi_out = response[:,4]
            
            for i in range(len(speed_out)):
                speed_out[i] = speed_out[i]+V0
                h_out[i] = h_out[i]+hp0
                theta_out[i] = theta_out[i]+th0
                
            plt.figure()
            plt.subplot(221)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Velocity [m/s]", fontsize = label_font)
            plt.plot(T,speed_out, label=fakelabel)
            plt.plot(time, V, label=reallabel)
            plt.legend()
            
            plt.subplot(222)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Altitude [m]", fontsize = label_font)
            plt.plot(T,h_out)
            plt.plot(time, h)
            
            plt.subplot(223)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Pitch angle [rad]", fontsize = label_font)
            plt.plot(T,theta_out)
            plt.plot(time, th)
            
            plt.subplot(224)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Elevator [rad]", fontsize = label_font)
            plt.plot(time, dE)
            
            plt.suptitle(modes[n], fontsize = title_font)  
            plt.show()
            
        if n==2 or n==3 or n==4 or n==5: #Asym
            #x = [u,alhpa,theta,q,h,beta,psi (roll),p,r,phi]
            #u = [de,da,dr]
            #Elevator set to zero since this shit is asymmetrical
            u_input=[]
            for i in range(len(dE)):
                u_input.append([0.,-dA[i],-dR[i]])
                
            response, T, state = co.lsim(SS, T = time,U = u_input)
            
            #plotting u,h,theta,psi(roll),phi(yaw)
            speed_out = response[:,0]
            h_out = response[:,1]
            theta_out = response[:,2]
            psi_out = response[:,3]
            phi_out = response[:,4]
            
            for i in range(len(speed_out)):
                speed_out[i] = speed_out[i]+V0
                h_out[i] = h_out[i]+hp0
                theta_out[i] = theta_out[i]+th0
                
            plt.figure()
            plt.subplot(231)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Velocity [m/s]", fontsize = label_font)
            plt.plot(T,speed_out, label=fakelabel)
            plt.plot(time, V, label=reallabel)
            plt.legend()
            
            plt.subplot(232)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Altitude [m]", fontsize = label_font)
            plt.plot(T,h_out)
            plt.plot(time, h)
            
            plt.subplot(233)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Yaw angle [rad]", fontsize = label_font)
            plt.plot(T,phi_out)
            plt.plot(time, yaw)
            
            plt.subplot(234)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Roll Angle [rad]", fontsize = label_font)
            plt.plot(T,psi_out)
            plt.plot(time, roll)
            
            plt.subplot(235)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Rudder Angle [rad]", fontsize = label_font)
            plt.plot(time, dR)
            
            plt.subplot(236)
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Aileron Angle [rad]", fontsize = label_font)
            plt.plot(time, dA)
            
            plt.suptitle(modes[n], fontsize = title_font)  
            plt.show()