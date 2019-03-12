# -*- coding: utf-8 -*-
"""
Created on Wed Mar 06 09:47:34 2019

@author: msjor
"""
from conversion_helpers import *
from Cit_par import *
import numpy as np


##input for me (Martijn to get my values to solve array by hand)
# hp = 5790.*ft_to_m  #pressure altitude
# Vc = 161.*kts_to_ms #indicated airspeed
# Tm = 5.             #measured airspeed in celsius

def true_airspeed(hp, Vc, Tm):
    Tm = Tm + 273.15
    Prat = 1 / ((1 + (lambdas * hp) / Temp0) ** (-g / (lambdas * R)))  # P0/P in ISA
    M = np.sqrt((2 / (gamma - 1)) * (
                (1 + Prat * ((1 + (gamma - 1) / (2 * gamma) * rho0 / P0 * Vc ** 2) ** (gamma / (gamma - 1)) - 1)) ** (
                    (gamma - 1) / gamma) - 1))

    Vt = M * (np.sqrt(gamma * R * Tm / (1 + (gamma - 1) / 2 * M ** 2)))
    return Vt
