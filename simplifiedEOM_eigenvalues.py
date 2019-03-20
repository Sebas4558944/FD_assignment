# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:46:19 2019

@author: msjor
"""

import numpy as np
from Cit_par_testing import muc,c,V0,Cmadot,KY2,CXu,CXa,CZa,CX0,CZq,Cmu,Cma,KX2,Cmq,mub,\
CYr,KXZ,b,Clr,Cnr,Clp,Cnp,CZadot,CZ0,CXq,CZu,CXde,CZde,Cmde,CYbdot,Cnbdot,KZ2,\
CYb,CL,CYp,Clb,Cnb,CYda,CYdr,Clda,Cldr,Cnda,Cndr, alpha0, th0
import scipy as sp

#from eigenvalues to damping ratio and natural frequency
def freq_damp_period(eigen):
    xi = np.real(eigen)
    eta = np.imag(eigen)
    freq = np.sqrt(xi**2+eta**2)
    damp = -xi/freq
    period = 2*np.pi/freq
    return freq, damp, period
    
#short period motion
def short_period():
    #V = constant and gamma = 0 which removes theta equation
    #CZadot and CZq << muc
    A = 2*muc*KY2*(2*muc-CZadot)
    B = -2*muc*KY2*CZa-(2*muc+CZq)*Cmadot-(2*muc-CZadot)*Cmq
    C = CZa*Cmq-(2*muc+CZq)*Cma
    eig1 =(-B+np.sqrt(B-4*A*C+0j))/2/A
    eig2 = (-B-np.sqrt(B-4*A*C+0j))/2/A
    return eig1, eig2

def phugoid():
    #no change in pitch rate and no in angle of attack
    #neglect CZq and CX0
    A = 2*muc*(CZa*Cmq-2*muc*Cma)
    B = 2*muc*(CXu*Cma-Cma*CXa)+Cmq*(CZa*CXa-CXu*CZa)
    C = CZ0*(Cmu*CZa-CZu*Cma)
    eig1 =(-B+np.sqrt(B-4*A*C+0j))/2/A
    eig2 = (-B-np.sqrt(B-4*A*C+0j))/2/A
    return eig1, eig2
    
def aperiodic_roll():
    #only a real eigenvalue as aperiodic
    #only rolls around longitudinal axis
    #only roll so no beta nd r 
    #direct result from rolling moment equation
    eig = Clp/(4*mub*KX2)
    return eig
    
def dutch_roll():
    #CYbdot and Cnbdot are neglected
    #CYr << mub
    A = 8*mub**2*KZ2
    B = -2*mub*(Cnr+2*KZ2*CYb)
    C = 4*mub*Cnb+CYb*Cnr
    eig1 =(-B+np.sqrt(B-4*A*C+0j))/2/A
    eig2 = (-B-np.sqrt(B-4*A*C+0j))/2/A
    return eig1, eig2
    
def spiral():
    #only a real eigenvalue as not aperiodic motion
    #CYr and CYp neglected 
    eig = 2*CL*(Clb*Cnr-Cnb*Clr)/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb))
    return eig


eigenlist = [["short period", short_period()],\
            ["Phugoid", phugoid()],\
            ["Aperiodic roll", aperiodic_roll()],\
            ["dutch_roll", dutch_roll()],\
            ["Spiral", spiral()]]
            
eigenlist_dimensional = [["short period", (short_period()[0]*V0/c, short_period()[1]*V0/c)],\
            ["Phugoid",  (phugoid()[0]*V0/c, phugoid()[1]*V0/c)],\
            ["Aperiodic roll", aperiodic_roll()*V0/b],\
            ["dutch_roll", (dutch_roll()[0]*V0/b, dutch_roll()[1]*V0/b)],\
            ["Spiral", spiral()*V0/b]]
            
def sym():
    A = 4*(muc**2)*KY2*(CZadot-2*muc)
    B = Cmadot*2*muc*(CZq+2*muc)-Cmq*2*muc*(CZadot-2*muc)-2*muc*KY2*(CXu*(CZadot-2*muc)-2*muc*CZa)
    C = Cma*2*muc*(CZq+2*muc)-Cmadot*(2*muc*CX0+CXu*(CZq+2*muc))+Cmq*(CXu*(CZadot-2*muc)-2*muc*CZa)+2*muc*KY2*(CXa*CZu-CZa*CXu)
    D = Cmu*(CXa*(CZq+2*muc)-CZ0*(CZadot-2*muc))-Cma*(2*muc*CX0+CXu*(CZq+2*muc))+Cmadot*(CX0*CXu-CZ0*CZu)+Cmq*(CXu*CZa-CZu*CXa)
    E = -Cmu*(CX0*CXa+CZ0*CZa)+Cma*(CX0*CXu+CZ0*CZu)
    R=B*C*D-A*D**2-B**2*E
    eigs = sp.roots(np.array([A,B,C,D,E]))
    return eigs,A,B,C,D,E, R



def asym():
    A=16*mub**3*(KX2*KZ2-KXZ**2)
    B=-4*mub**2*(2*CYb*(KX2*KZ2-KXZ**2)+Cnr*KX2+Clp*KZ2+(Clr+Cnp)*KXZ)
    C=2*mub*((CYb*Cnr-CYr*Cnb)*KX2+(CYb*Clp-Clb*CYp)*KZ2+((CYb*Cnp-Cnb*CYp)+(CYb*Clr-Clb*CYr))*KXZ+4*mub*Cnb*KX2+4*mub*Clb*KXZ+0.5*(Clp*Cnr-Cnp*Clr))
    D=-4*mub*CL*(Clb*KZ2+Cnb*KXZ)+2*mub*(Clb*Cnp-Cnb*Clp)+0.5*CYb*(Clr*Cnp-Cnr*Clp)+0.5*CYp*(Clb*Cnr-Cnb*Clr)+0.5*CYr*(Clp*Cnb-Cnp*Clb)
    E=CL*(Clb*Cnr-Cnb*Clr)
    R=B*C*D-A*D**2-B**2*E
    eigs = sp.roots(np.array([A,B,C,D,E]))
    return eigs,A,B,C,D,E, R

