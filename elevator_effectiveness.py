from Cit_par import *
from reduced_condition_calculator import Conditions
from datareader import getFDValues

import numpy as np
import matplotlib.pyplot as plt

reduction_calculator = Conditions(5000)

Cm0 = 0.0297
CNa = CNwa + CNha


def calc_CN(weight, density, v, s):
    return weight / (0.5 * density * (v ** 2) * s)


def calc_DCm(x_cg_2, x_cg_1, C_N, chord):
    return C_N * ((x_cg_2 - x_cg_1) / chord)


def calc_Cm_delta(diff_delta, diff_Cm):
    return -(diff_Cm / diff_delta)


def calc_delta_reduced(delta_measured, Cm_delta, Tcs, Tc):
    CmTc = - 0.0064
    delta_reduced = delta_measured - (1 / Cm_delta) * CmTc * (Tcs - Tc)

    return delta_reduced


key_list, description_list, unit_list, newDict = getFDValues('reference.mat')

alpha = newDict.get('vane_AOA', {})
delta_e = newDict.get('delta_e', {})

plt.plot()

print alpha
print delta_e