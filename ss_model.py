# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 14:59:51 2019

@author: msjor
"""

import control.matlab as co
import numpy as np
from Cit_par import muc,c,V0,Cmadot,KY2,CXu,CXa,CZa,CX0,CZq,Cmu,Cma,KX2,Cmq,mub,\
CYr,KXZ,b,Clr,Cnr,Clp,Cnp,CZadot,CZ0,CXq,CZu,CXde,CZde,Cmde,CYbdot,Cnbdot,KZ2,\
CYb,CL,CYp,Clb,Cnb,CYda,CYdr,Clda,Cldr,Cnda,Cndr


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


C1S=np.matrix([[-2*muc*c, 0, 0, 0],\
                [0, (CZadot - 2 * muc) * c/V0, 0, 0],
                [0, 0, c/V0, 0,],
                [0, Cmadot * c/V0, 0, -2*muc*KY2*(c/V0)**2]])

C2S=np.matrix([[CXu*V0, CXa, CZ0 ,CXq * c/V0],\
              [CZu*V0, CZa, -CX0, (CZq+2*muc) * (c/V0)],\
              [0 , 0 , 0,  c/V0],\
              [Cmu*V0, Cma, 0, Cmq*c/V0]])

C3S=np.matrix([[CXde],\
                [CZde],\
                [0],\
                [Cmde]])

C1A=np.matrix([[(CYbdot-2*mub)*b/V0, 0, 0, 0,0],\
                [0, -0.5*b/V0, 0, 0,0],\
                [0, 0, -2*mub*KX2*(b/V0)**2, 2*mub*KXZ*(b/V0)**2,0],\
                [Cnbdot*b/V0, 0, 2*mub*KXZ*(b/V0)**2, -2*mub*KZ2*(b/V0)**2,0],\
                [0,0,0,0,-0.5*b/V0]])

C2A=np.matrix([[CYb, CL, CYp*b/(2*V0), (CYr-4*mub)*b/(2*V0),0],\
                [0, 0, b/(2*V0), 0,0],\
                [Clb, 0, Clp*b/(2*V0), Clr*b/(2*V0),0],\
                [Cnb, 0, Cnp*b/(2*V0), Cnr*b/(2*V0),0],\
                [0,0,0,b/(2*V0),0]])

C3A=np.matrix([[CYda,CYdr],\
                [0,0],\
                [Clda,Cldr],\
                [Cnda,Cndr],\
                [0,0]])

AS=np.linalg.inv(C1S)*-C2S
BS=np.linalg.inv(C1S)*-C3S

AA=np.linalg.inv(C1A)*-C2A
BA=np.linalg.inv(C1A)*-C3A



#Outputs
#Symmetric outputs -> u, theta
#x=[u,a,theta,p]T
#y=[u,theta]T
#asymetric outputs -> psi, phi
#x=[beta, psi, p,r,phi]
#y=[psi, phi]

CS=np.matrix([[1,0,0,0],\
              [0,0,1,0]])

DS=np.zeros((2,1))


CA=np.matrix(([0,1,0,0,0],\
              [0,0,0,0,1]))

DA=np.zeros((2,2))

#init ss
SSS=co.ss(AS,BS,CS,DS)
SSA=co.ss(AA,BA,CA,DA)



