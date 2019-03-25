# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 14:59:51 2019

@author: msjor
"""
import matplotlib.pyplot as plt
import control.matlab as co
import numpy as np
from Cit_par import muc,c,hp0, V0,Cmadot,KY2,CXu,CXa,CZa,CX0,CZq,Cmu,Cma,KX2,Cmq,mub,\
CYr,KXZ,b,Clr,Cnr,Clp,Cnp,CZadot,CZ0,CXq,CZu,CXde,CZde,Cmde,CYbdot,Cnbdot,KZ2,\
CYb,CL,CYp,Clb,Cnb,CYda,CYdr,Clda,Cldr,Cnda,Cndr, alpha0, th0


#\theta=angle between X and horizontal
#Dc,Db=Differnetial operator
#D_c=c/V0 #*d/dt
#D_b=b/V0 #*d/dt


#Symmetric flight
#C1*xdot+C2*x+C3*u=0
#C1  [(du/dt), (da/dt), (dθ/dt),(dq/dt)]
#-2*muc*c, 0, 0, 0
#0, (Czad - 2 muc) * c/V0, 0, 0
#0, 0, c/V0, 0
#0, Cmad * c/V0, 0, -2muc*ky2* (c/V0)**2
#
#C2   * [u, a, θ, q]
#Cxu * V0, Cxa, Cz0 ,Cxq * c/V0
#Czu*V0, Cza, -Cx0, (Czq+2muc) * (c/V0) 
#0 , 0 , 0,  c/V0
#Cmu*V0, Cma, 0, Cmq*c/V0
#
#C3 [de]
#Cxde
#Czde
#0
#Cmde

#Asymmetric flight
# x = [beta,phi,p roll, r yaw]
# u = [da,dr] 
#C1
#(Cybd-2mub)*b/V0, 0, 0, 0
#0, -0.5*b/V0, 0, 0
#0,0,-2*mub*KX2*(b/V0)^2, 2*mub*KXZ2*(b/V0)^2 
#Cnbd*b/V0,0, 2*mub*KXZ2*(b/V0)^2, -2*mub*KZ2*(b/V0)^2 
#C2
#Cyb, CL, Cyp*b/2V0, (Cyr-4mub)*b/2V0
#0, 0, b/2V0, 0
#Clb, 0, Clp*b/2V0, Clr*b/2V0
#Cnb, 0, Cnp*b/2V0, Cnr*b/2V0
#C3
#Cyda,Cydr
#0,0
#Clda,Cldr
#Cnda,Cndr


#xdot=(C1)^-1 * -C2*x - (C1)^-1 * -C3*u
#A=C1^-1*-C2
#B=C1^-1*-C2

# combining state vectors gives:
#x = [u,alhpa,theta,q,h,beta,psi,p,r,phi]
#u = [de,da,dr]
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

C3S=np.matrix([[CXde],\
                [CZde],\
                [0.],\
                [Cmde],\
                [0.],\
                [0.],\
                [0.],\
                [0.],\
                [0.],\
                [0.]])

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

C3A=np.matrix([[0,0],\
               [0,0],\
               [0,0],\
               [0,0],\
               [0,0],\
               [0,0],\
               [CYda,CYdr],\
               [0,0],\
               [Clda,Cldr],\
               [Cnda,Cndr],\
               [0,0]])

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

#h = co.tf(SS)

#print eigenvalues A matrix
eigenvals, eigenvectors = np.linalg.eig(A)
print eigenvals
#######plotting responses
label_font = 20
title_font = 25

#############################inital value problem ################################
#T = np.linspace(0,100,1000)
#X0 = np.matrix([[0],\
#                [0],\
#                [0],\
#                [0],\
#                [0],\
#                [0],\
#                [0],\
#                [0],\
#                [0],\
#                [0]])
#
#response, T = co.initial(SS, T = T,X0 = X0)
##plotting u,h,theta,psi,phi
#u = response[:,0]
#h = response[:,1]
#theta = response[:,2]
#psi = response[:,3]
#phi = response[:,4]
#
#plt.figure()
#plt.subplot(231)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Velocity [m/s]", fontsize = label_font)
#plt.plot(T,u)
#
#plt.subplot(232)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Altitude [m]", fontsize = label_font)
#plt.plot(T,h)
#
#plt.subplot(233)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Pitch angle [rad]", fontsize = label_font)
#plt.plot(T,theta)
#
#plt.subplot(234)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Roll angle [rad]", fontsize = label_font)
#plt.plot(T,psi)
#
#plt.subplot(235)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Yaw angle [rad]", fontsize = label_font)
#plt.plot(T,phi)
#plt.show()

#################response to impulse on control vector#####################
#steps = 1000
#tmax = 100 
#T = np.linspace(0,tmax,steps)
#
##create impulse vector for t = 0 
#u = []
##[de,da,dr]
#u_impulse = [0.1,0.,0.]
#u.append(u_impulse)

#move forcing to 0 for anything past the initial input
#for i in range(steps-1):
#    u.append([0.,0.,0.])
#u = np.array(u)
#
##calculate response
#response, T, state = co.lsim(SS, T = T,U = u)
#
##plotting u,h,theta,psi,phi
#u = response[:,0]
#h = response[:,1]
#theta = response[:,2]
#psi = response[:,3]
#phi = response[:,4]
#
#plt.figure()
#plt.subplot(231)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Velocity [m/s]", fontsize = label_font)
#plt.plot(T,u)
#
#plt.subplot(232)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Altitude [m]", fontsize = label_font)
#plt.plot(T,h)
#
#plt.subplot(233)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Pitch angle [rad]", fontsize = label_font)
#plt.plot(T,theta)
#
#plt.subplot(234)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Roll angle [rad]", fontsize = label_font)
#plt.plot(T,psi)
#
#plt.subplot(235)
#plt.xlabel("Time [sec]", fontsize = label_font)
#plt.ylabel("Yaw angle [rad]", fontsize = label_font)
#plt.plot(T,phi)
#plt.show()


#############step input from t=0 to t=tstep ###################
steps = 1000
tmax = 20. 
tstep = 1.
nstep = tstep/(tmax/float(steps))
T = np.linspace(0,tmax,steps)

#create impulse vector for t = 0 
u_input = []
#[de,da,dr]
u_val = [0.0,0.0,0.025]

#move forcing to 0 for anything past the initial input
for i in range(steps): 
    if i<=nstep:
        u_input.append(u_val)
    elif i>nstep:
        u_input.append([0.,0.,0.0])
u_input = np.array(u_input)

#calculate response
response, T, state = co.lsim(SS, T = T,U = u_input)

#plotting u,h,theta,psi,phi
u = response[:,0]
h = response[:,1]
theta = response[:,2]
psi = response[:,3]
phi = response[:,4]

for i in range(len(u)):
    u[i] = u[i]+V0
    h[i] = h[i]+hp0
    
#plt.figure()
plt.subplot(231)
plt.xlabel("Time [sec]", fontsize = label_font)
plt.ylabel("Velocity [m/s]", fontsize = label_font)
plt.plot(T,u)

plt.subplot(232)
plt.xlabel("Time [sec]", fontsize = label_font)
plt.ylabel("Altitude [m]", fontsize = label_font)
plt.plot(T,h)

plt.subplot(233)
plt.xlabel("Time [sec]", fontsize = label_font)
plt.ylabel("Pitch angle [rad]", fontsize = label_font)
plt.plot(T,theta)

plt.subplot(234)
plt.xlabel("Time [sec]", fontsize = label_font)
plt.ylabel("Roll angle [rad]", fontsize = label_font)
plt.plot(T,psi)

plt.subplot(235)
plt.xlabel("Time [sec]", fontsize = label_font)
plt.ylabel("Yaw angle [rad]", fontsize = label_font)
plt.plot(T,phi)
plt.show()
