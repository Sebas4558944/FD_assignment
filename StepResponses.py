# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 18:47:44 2019

@author: daanv
"""



import matplotlib.pyplot as plt
import control.matlab as co
import numpy as np

from Cit_par_default import muc,c,alpha0, V0,th0,Cmadot,KY2,CXu,CXa,CZa,CX0,CZq,Cmu,Cma,KX2,Cmq,mub,\
        CYr,KXZ,b,Clr,Cnr,Clp,Cnp,CZadot,CZ0,CXq,CZu,CXde,CZde,Cmde,CYbdot,Cnbdot,KZ2,\
        CYb,CL,CYp,Clb,Cnb,CYda,CYdr,Clda,Cldr,Cnda,Cndr
#Gets spiral stability
Clb = Clb*1.3
Clr = Clr*0.7

Cnb = Cnb*0.7
Cnr = Cnr*1.3

def getStateSpace():
    # combining state vectors gives:
    #x = [u,alhpa,theta,q,h,beta,psi,p,r,phi]
    #u = [de,da,dr]
    #Changed indicates whether the default or the changed cit_par is used
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
#print CXu, CXa, CXde, CXq
#print Cma,Cmadot,Cmq,Cmde,Cmu
#print CZa, CZq, CZadot, CZ0, CZu, CZde
#
#print CYr, CYp, CYdr, CYda, CYb
#print CL, Clr, Cldr, Clda, Clb, Clp, Cnp
#print Cnr, Cnb, Cnda, Cndr 
niks=1
sym='CXu,CXa,CXde,CXq,Cma,Cmadot,Cmq,Cmde,Cmu,CZa,CZq,CZadot,CZ0,CZu,CZde,niks' #List for plot titles
symList=[CXu,CXa,CXde,CXq,Cma,Cmadot,Cmq,Cmde,Cmu,CZa,CZq,CZadot,CZ0,CZu,CZde,niks] #List to be changed
symNamesList=sym.split(',')

asym='CYr,CYp,CYdr,CYda,CYb,CL,Clr,Cldr,Clda,Clb,Clp,Cnp,Cnr,Cnb,Cnda,Cndr,niks'
asymList=[CYr,CYp,CYdr,CYda,CYb,CL,Clr,Cldr,Clda,Clb,Clp,Cnp,Cnr,Cnb,Cnda,Cndr,niks]
asymNamesList=asym.split(',')

label_font = 20
title_font = 25

plotSym=0
plotAsym=1

if plotSym:
    for i in range(len(symList)):
    #    print symNamesList[i]+'=-'+str(symList[i])   #Check string
        #Invert command
        exec(symNamesList[i]+'=-'+str(symList[i]))
    #    exec('print '+str(symNamesList[i]))  #Check inversion
        SS=getStateSpace()
        co.step(SS)
        
        #############step input from t=0 to t=tstep ###################
        steps = 1000*4
        tmax = 100.*4
        tstep = 10. #10 sec step input
        nstep = tstep/(tmax/float(steps))
        T = np.linspace(0,tmax,steps)
        
        #create impulse vector for t = 0 
        u_input = []
        #[de,da,dr]
        u_val = [0.01,0.0,0.0]
        
        #move forcing to 0 for anything past the initial input
        for o in range(steps):
            if o<=nstep:
                u_input.append(u_val)
            elif o>nstep:
                u_input.append([0.,0.,0.])
        u_input = np.array(u_input)
        
        #calculate response
        response, T, state = co.lsim(SS, T = T,U = u_input)
        speed_out = response[:,0]
        h_out = response[:,1]
        theta_out = response[:,2]
        psi_out = response[:,3]
        phi_out = response[:,4]
        
        plt.figure()
            
        plt.subplot(221)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Velocity [m/s]", fontsize = label_font)
        plt.plot(T,speed_out)
        
        plt.subplot(222)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Altitude [m]", fontsize = label_font)
        plt.plot(T,h_out)
        
        plt.subplot(223)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Pitch angle [rad]", fontsize = label_font)
        plt.plot(T,theta_out)
        
        plt.subplot(224)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Elevator [rad]", fontsize = label_font)
        plt.plot(T, u_input[:,0])
        
        plt.suptitle(('Variable changed: '+symNamesList[i]), fontsize = title_font)  
        plt.show()
    
            #Revert to normal
        exec(symNamesList[i]+'='+str(symList[i]))
    #    exec('print '+str(symNamesList[i]))  #Check revert

if plotAsym:
    for i in range(len(asymList)):
    #    print asymNamesList[i]+'=-'+str(asymList[i])   #Check string
        #Invert command
        exec(asymNamesList[i]+'=-'+str(asymList[i]))
    #    exec('print '+str(asymNamesList[i]))  #Check inversion
        SS=getStateSpace()
        co.step(SS)
        
        #############step input from t=0 to t=tstep ###################
        steps = 1000*5
        tmax = 100.*5. 
        tstep = 10. #10 sec step input
        nstep = tstep/(tmax/float(steps))
        T = np.linspace(0,tmax,steps)
        
        #create impulse vector for t = 0 
        u_input = []
        #[de,da,dr]
        u_val = [0.0,0.01,0.01]
        
        #move forcing to 0 for anything past the initial input
        for o in range(steps):
            if o<=nstep:
                u_input.append(u_val)
            elif o>nstep:
                u_input.append([0.,0.,0.])
        u_input = np.array(u_input)
        
        #calculate response
        response, T, state = co.lsim(SS, T = T,U = u_input)
        speed_out = response[:,0]
        h_out = response[:,1]
        theta_out = response[:,2]
        psi_out = response[:,3]
        phi_out = response[:,4]
        
        plt.figure()
        plt.subplot(231)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Velocity [m/s]", fontsize = label_font)
        plt.plot(T,speed_out)
        
        plt.subplot(232)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Altitude [m]", fontsize = label_font)
        plt.plot(T,h_out)
        
        plt.subplot(233)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Yaw angle [rad]", fontsize = label_font)
        plt.plot(T,phi_out)
        
        plt.subplot(234)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Roll Angle [rad]", fontsize = label_font)
        plt.plot(T,psi_out)
        
        plt.subplot(235)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Rudder Angle [rad]", fontsize = label_font)
        plt.plot(T, u_input[:,2])
        
        plt.subplot(236)
        plt.grid()
        plt.xlabel("Time [sec]", fontsize = label_font)
        plt.ylabel("Aileron Angle [rad]", fontsize = label_font)
        plt.plot(T, u_input[:,1])
        
        plt.suptitle(('Variable changed: '+asymNamesList[i]), fontsize = title_font)  
        plt.show()
    
        #Revert to normal
        exec(asymNamesList[i]+'='+str(asymList[i]))
        #exec('print '+str(asymNamesList[i]))  #Check revert
