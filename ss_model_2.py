# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 14:59:51 2019

@author: msjor
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
                    [0,0,0,1,0,0,0,0,0,0]])
    DS = np.zeros((8,3))
    #
    ##init ss
    SS=co.ss(A,B,CS,DS)
    #print eigenvalues A matrix dimonsionalized
    eigenvals, eigenvectors = np.linalg.eig(A)
    print eigenvals, changed
    return SS

#######plotting vars
label_font = 20
title_font = 25
reallabel='FlightData'
fakelabel='Model Data with Original Parameters'
fixedfakelabel='Model Data with Adapted Parameters'
modes = ["Phugoid", "Short Period", "Dutch Roll","Dutch Roll Yd", "Aperiodic Roll", "Spiral" ]

plotting=[0,0,1,0,0,0]
printRMS=1
indices, times, altitudes, velocities, alphas,  pitches, rolls, yaws, ailerons, rudders, elevators, rollrates,pitchrates = getEigenmotions()
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
        rollr = rollrates[n]
        pitchr = pitchrates[n]
        
        #set init values
        hp0=h[0]
        V0=V[0]
        alpha0=A[0]
        th0=th[0]
        r0=roll[0]
        y0=yaw[0]
        rr0 = rollr[0]        
        pr0 = pitchr[0]
        
        SS=getStateSpace(alpha0,V0,th0,0) #Default Values State Space
        SSF=getStateSpace(alpha0,V0,th0,1) #Fixed State Space 
        
        if n==0 or n==1: #Sym
            #x = [u,alhpa,theta,q,h,beta,psi,p,r,phi]
            #u = [de,da,dr]
            #Rudders and ailerons set to zero since this shit is symmetrical
            u_input=[]
            dE0=dE[0]
            for i in range(len(dE)):
                u_input.append([(dE[i]-dE0),0.,0.])
                
            response, T, state = co.lsim(SS, T = time,U = u_input)
            responseF, TF, stateF = co.lsim(SSF, T = time,U = u_input)
            
            #plotting u,h,theta,psi(roll),phi(yaw)
            speed_out = response[:,0]
            h_out = response[:,1]
            theta_out = response[:,2]
            psi_out = response[:,3]
            phi_out = response[:,4]
            rollr_out = response[:,5]
            AoA_out = response[:,6]
            pitchr_out = response[:,7]
            
            speed_outF = responseF[:,0]
            h_outF = responseF[:,1]
            theta_outF = responseF[:,2]
            psi_outF = responseF[:,3]
            phi_outF = responseF[:,4]
            rollr_outF = responseF[:,5]
            AoA_outF = responseF[:,6]
            pitchr_outF = responseF[:,7]
            
            speed_delta=0
            theta_delta=0
            h_delta=0
            
            q=0
            for i in range(len(speed_out)):
                q+=1
                speed_out[i] = speed_out[i]+V0
                h_out[i] = h_out[i]+hp0
                theta_out[i] = theta_out[i]+th0
                AoA_out[i] = AoA_out[i]+alpha0
                pitchr_out[i] = pitchr_out[i]+pr0                
                
                speed_outF[i] = speed_outF[i]+V0
                h_outF[i] = h_outF[i]+hp0
                theta_outF[i] = theta_outF[i]+th0
                AoA_outF[i] = AoA_outF[i]+alpha0
                pitchr_outF[i] = pitchr_outF[i]+pr0 
                
                speed_delta+=(V[i]-speed_outF[i])**2              
                theta_delta+=(th[i]-theta_outF[i])**2
                h_delta+=(h[i]-h_outF[i])**2
            
            speed_delta=np.sqrt(speed_delta/(q-1))
            theta_delta=np.sqrt(theta_delta/(q-1))
            h_delta=np.sqrt(h_delta/(q-1))
            
            if printRMS:
                print ''
                print modes[n]
                print "V rms=",speed_delta
                print "theta rms=",theta_delta
                print "h rms=",h_delta
            
            plt.figure()
            plt.subplot(231)
            plt.grid()
            plt.xlabel("t [sec]", fontsize = label_font)
            plt.ylabel("V [m/s]", fontsize = label_font)
            plt.plot(T,speed_out, label=fakelabel)
            plt.plot(TF,speed_outF, label=fixedfakelabel)
            plt.plot(time, V, label=reallabel)
            plt.legend(loc = "lower left")
            
            plt.subplot(232)
            plt.grid()
            plt.xlabel("t [sec]", fontsize = label_font)
            plt.ylabel("h [m]", fontsize = label_font)
            plt.plot(T,h_out, label=fakelabel)
            plt.plot(TF,h_outF, label=fixedfakelabel)
            plt.plot(time, h, label=reallabel)
            
            plt.subplot(233)
            plt.grid()
            plt.xlabel("t [sec]", fontsize = label_font)
            plt.ylabel(r"$\theta$ [rad]", fontsize = label_font)
            plt.plot(T,theta_out, label=fakelabel)
            plt.plot(TF,theta_outF, label=fixedfakelabel)
            plt.plot(time, th, label=reallabel)
            
            plt.subplot(234)
            plt.grid()
            plt.xlabel("t [sec]", fontsize = label_font)
            plt.ylabel(r"$\alpha$ [rad]", fontsize = label_font)
            plt.plot(T,AoA_out, label=fakelabel)
            plt.plot(TF,AoA_outF, label=fixedfakelabel)
            plt.plot(time, A, label=reallabel)
            plt.suptitle(modes[n], fontsize = title_font)  
            plt.show()
            
            plt.subplot(235)
            plt.grid()
            plt.xlabel("t [sec]", fontsize = label_font)
            plt.ylabel(r"$\delta_e$ [rad]", fontsize = label_font)
            plt.plot(time, dE, label=reallabel)

#            plt.subplot(236)
#            plt.grid()
#            plt.xlabel("Time [sec]", fontsize = label_font)
#            plt.ylabel("Pitch rate [rad/s]", fontsize = label_font)
#            plt.plot(T,pitchr_out, label=fakelabel)
#            plt.plot(TF,pitchr_outF, label=fixedfakelabel)
#            plt.plot(time, th, label=reallabel)            
            
            plt.suptitle(modes[n], fontsize = title_font)  
            plt.show()
            
        if n==2 or n==3 or n==4 or n==5: #Asym
            #x = [u,alhpa,theta,q,h,beta,psi (roll),p,r,phi]
            #u = [de,da,dr]
            #Elevator set to zero since this shit is asymmetrical
            u_input=[]
            dA0=dA[0]
            dR0=dR[0]
            for i in range(len(dE)):
                u_input.append([0.,(dA[i]-dA0),(dR[i]-dR0)])
                
            response, T, state = co.lsim(SS, T = time,U = u_input)
            responseF, TF, stateF = co.lsim(SSF, T = time,U = u_input)
            
            #plotting u,h,theta,psi(roll),phi(yaw)
            speed_out = response[:,0]
            h_out = response[:,1]
            theta_out = response[:,2]
            psi_out = response[:,3]*-1.
            phi_out = response[:,4]*-1.
            rollr_out = response[:,5]*-1            
            
            speed_outF = responseF[:,0]
            h_outF = responseF[:,1]
            theta_outF = responseF[:,2]
            psi_outF = responseF[:,3]*-1.
            phi_outF = responseF[:,4]*-1.
            rollr_outF = responseF[:,5]*-1.
            
            psi_delta=0
            phi_delta=0
            rollr_delta=0
            
            q=0
            
            for i in range(len(speed_out)):
                q+=1
                
                speed_out[i] = speed_out[i]+V0
                h_out[i] = h_out[i]+hp0
                theta_out[i] = theta_out[i]+th0
                psi_out[i] = psi_out[i]+r0
                phi_out[i] = phi_out[i]+y0
                rollr_out[i] = rollr_out[i]+rr0
                
                speed_outF[i] = speed_outF[i]+V0
                h_outF[i] = h_outF[i]+hp0
                theta_outF[i] = theta_outF[i]+th0
                psi_outF[i] = psi_outF[i]+r0
                phi_outF[i] = phi_outF[i]+y0
                rollr_outF[i] = rollr_outF[i]+rr0
                           
                psi_delta+=(roll[i]-psi_outF[i])**2
                phi_delta+=(yaw[i]-phi_outF[i])**2
                rollr_delta+=(rollr[i]-rollr_outF[i])**2
                
            psi_delta=np.sqrt(psi_delta/(q-1))
            phi_delta=np.sqrt(phi_delta/(q-1))
            rollr_delta=np.sqrt(rollr_delta/(q-1))
            
            if printRMS:
                print ''
                print modes[n]
                print 'Psi/roll rms=', psi_delta
                print 'Phi/yaw rms=', phi_delta
                print 'psidot/rollrate rms=', rollr_delta
            
            plt.figure()
#            plt.subplot(231)
#            plt.grid()
#            plt.xlabel("Time [sec]", fontsize = label_font)
#            plt.ylabel("Velocity [m/s]", fontsize = label_font)
#            plt.plot(T,speed_out, label=fakelabel)
#            plt.plot(TF,speed_outF, label=fixedfakelabel)
#            plt.plot(time, V, label=reallabel)

            
            plt.subplot(231)
            plt.grid()
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Roll rate [rad/s]", fontsize = label_font)
#            plt.plot(T,rollr_out, label=fakelabel)
            plt.plot(TF,rollr_outF, label=fixedfakelabel)
            plt.plot(time, rollr, label=reallabel)
            plt.legend(loc = "lower left")
            
            plt.subplot(232)
            plt.grid()
            plt.xlabel("Time [sec]", fontsize = label_font)
            plt.ylabel("Yaw rate [rad/s]", fontsize = label_font)
#            plt.plot(T,phi_out, label=fakelabel)
            plt.plot(TF,phi_outF, label=fixedfakelabel)
            plt.plot(time, yaw, label=reallabel)
            
            plt.subplot(233)
            plt.grid()
            plt.xlabel("t [sec]", fontsize = label_font)
            plt.ylabel(r"$\psi$ [rad]", fontsize = label_font)
#            plt.plot(T,psi_out, label=fakelabel)
            plt.plot(TF,psi_outF, label=fixedfakelabel)
            plt.plot(time, roll, label=reallabel)
            
            plt.subplot(234)
            plt.grid()
            plt.xlabel("t [sec]", fontsize = label_font)
            plt.ylabel(r"$\delta_r$ [rad]", fontsize = label_font)
            plt.plot(time, dR, label=reallabel)
            
            plt.subplot(235)
            plt.grid()
            plt.xlabel("t [sec]", fontsize = label_font)
            plt.ylabel(r"$\delta_a$ [rad]", fontsize = label_font)
            plt.plot(time, dA, label=reallabel)
            
            plt.suptitle(modes[n], fontsize = title_font)  
            plt.show()