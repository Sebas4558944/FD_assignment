# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:46:19 2019

@author: msjor
"""

import numpy as np
from Cit_par import muc,c,V0,Cmadot,KY2,CXu,CXa,CZa,CX0,CZq,Cmu,Cma,KX2,Cmq,mub,\
CYr,KXZ,b,Clr,Cnr,Clp,Cnp,CZadot,CZ0,CXq,CZu,CXde,CZde,Cmde,CYbdot,Cnbdot,KZ2,\
CYb,CL,CYp,Clb,Cnb,CYda,CYdr,Clda,Cldr,Cnda,Cndr, alpha0, th0


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