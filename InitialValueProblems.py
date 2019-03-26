# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 16:18:56 2019

@author: daanv
"""

import matplotlib.pyplot as plt
import control.matlab as co
import numpy as np
from Eigenplotter_Func import getEigenmotions

def getStateSpace(alpha0,V0,th0,changed):
    # combining state vectors gives:
    #x = [u,alhpa,theta,q,h,beta,psi,p,r,phi]
    #u = [de,da,dr]
    #Changed indicates whether the default or the changed cit_par is used
    if changed==1:
        from Cit_par import muc,c,Cmadot,KY2,CXu,CXa,CZa,CX0,CZq,Cmu,Cma,KX2,Cmq,mub,\
        CYr,KXZ,b,Clr,Cnr,Clp,Cnp,CZadot,CZ0,CXq,CZu,CXde,CZde,Cmde,CYbdot,Cnbdot,KZ2,\
        CYb,CL,CYp,Clb,Cnb,CYda,CYdr,Clda,Cldr,Cnda,Cndr, alpha0, th0
    else:
        from Cit_par_default import muc,c,Cmadot,KY2,CXu,CXa,CZa,CX0,CZq,Cmu,Cma,KX2,Cmq,mub,\
        CYr,KXZ,b,Clr,Cnr,Clp,Cnp,CZadot,CZ0,CXq,CZu,CXde,CZde,Cmde,CYbdot,Cnbdot,KZ2,\
        CYb,CL,CYp,Clb,Cnb,CYda,CYdr,Clda,Cldr,Cnda,Cndr, alpha0, th0

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
                    [0,0,0,0,0,0,0,0,1,0],\
                    [0,0,0,0,0,0,0,1,0,0],\
                    [0,1,0,0,0,0,0,0,0,0],\
                    [0,0,0,0,0,0,0,0,0,1]])
    DS = np.zeros((8,3))
    
    #
    ##init ss
    SS=co.ss(A,B,CS,DS)
    #print eigenvalues A matrix dimonsionalized
    eigenvals, eigenvectors = np.linalg.eig(A)

    return SS

#######plotting vars
label_font = 20
title_font = 25
reallabel='FlightData'
fakelabel='Model Data with Original Parameters'
fixedfakelabel='Model Data with Adapted Parameters'



#Get Flight Data        
#set init values
alpha0=0.
V0=180.
th0=0.
hp0=5000.

SS=getStateSpace(alpha0,V0,th0,0) #Default Values State Space
SSF=getStateSpace(alpha0,V0,th0,1) #Fixed State Space 

steps = 1000*4+1
tmax = 100.*4 
T = np.linspace(0,tmax,steps)

plotting=[1,1,1,1,1,1,1,1,1,1]
modes=['Velocity','alpha','theta','pitch rate','height','Sideslip','psi','p','r','phi' ]
DX0 = np.matrix([[1],\
                [0.05],\
                [0.05],\
                [0.05],\
                [1],\
                [0.05],\
                [0.05],\
                [0.05],\
                [0.05],\
                [0.05]])

    
for n in range(10): #Sym
    if plotting[n]:
        #x = [u,alhpa,theta,q,h,beta,psi,p,r,phi]
        #u = [de,da,dr]
        #Rudders and ailerons set to zero since this shit is symmetrical
        X0 = np.matrix([[0],\
                        [0],\
                        [0],\
                        [0],\
                        [0],\
                        [0],\
                        [0],\
                        [0],\
                        [0],\
                        [0]])
        X0[n]=DX0[n]
        
        response, T = co.initial(SS, T = T,X0 = X0)
        responseF, TF = co.initial(SSF, T = T,X0 = X0)
    
        #plotting u,h,theta,psi(roll),phi(yaw)
        speed_out = response[:,0]
        h_out = response[:,1]
        theta_out = response[:,2]
        psi_out = response[:,3]
        phi_out = response[:,4]
        rollr_out = response[:,5]
        AoA_out = response[:,6]
        
        speed_outF = responseF[:,0]
        h_outF = responseF[:,1]
        theta_outF = responseF[:,2]
        psi_outF = responseF[:,3]
        phi_outF = responseF[:,4]
        rollr_outF = responseF[:,5]
        AoA_outF = responseF[:,6]
        
        for i in range(len(speed_out)):
            speed_out[i] = speed_out[i]+V0
            h_out[i] = h_out[i]+hp0
            theta_out[i] = theta_out[i]+th0
            AoA_out[i] = AoA_out[i]+alpha0
            
            speed_outF[i] = speed_outF[i]+V0
            h_outF[i] = h_outF[i]+hp0
            theta_outF[i] = theta_outF[i]+th0
            AoA_outF[i] = AoA_outF[i]+alpha0
        
        plt.figure()
        plt.subplot(221)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Velocity [m/s]", fontsize = label_font)
        plt.plot(T,speed_out, label=fakelabel)
        plt.plot(TF,speed_outF, label=fixedfakelabel)
        plt.legend()
        
        plt.subplot(222)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Altitude [m]", fontsize = label_font)
        plt.plot(T,h_out, label=fakelabel)
        plt.plot(TF,h_outF, label=fixedfakelabel)
        
        plt.subplot(223)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Pitch angle [rad]", fontsize = label_font)
        plt.plot(T,theta_out, label=fakelabel)
        plt.plot(TF,theta_outF, label=fixedfakelabel)
        
        plt.subplot(224)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("AoA [rad]", fontsize = label_font)
        plt.plot(T,AoA_out, label=fakelabel)
        plt.plot(TF,AoA_outF, label=fixedfakelabel)
        plt.suptitle(modes[n], fontsize = title_font)  
        plt.show()
    
        plt.suptitle(modes[n], fontsize = title_font)  
        plt.show()
        
        plt.figure()
        plt.subplot(221)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Velocity [m/s]", fontsize = label_font)
        plt.plot(T,speed_out, label=fakelabel)
        plt.plot(TF,speed_outF, label=fixedfakelabel)
        plt.legend()
        
        plt.subplot(222)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Roll rate [rad/s]", fontsize = label_font)
        plt.plot(T,rollr_out, label=fakelabel)
        plt.plot(TF,rollr_outF, label=fixedfakelabel)
        
        plt.subplot(223)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Yaw rate [rad/s]", fontsize = label_font)
        plt.plot(T,phi_out, label=fakelabel)
        plt.plot(TF,phi_outF, label=fixedfakelabel)
        
        plt.subplot(224)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Roll Angle [rad]", fontsize = label_font)
        plt.plot(T,psi_out, label=fakelabel)
        plt.plot(TF,psi_outF, label=fixedfakelabel)
        
        plt.suptitle(modes[n], fontsize = title_font)  
        plt.show()