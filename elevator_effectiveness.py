from Cit_par import *
from reduced_condition_calculator import Conditions
from datareader import importExcelData

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


def calc_derivative(alpha_1, delta_1, alpha_2, delta_2):
        return (delta_2-delta_1)/(alpha_2 - alpha_1)


f='Reference_Datasheet.csv'
date_of_flight, flight_number, TO_time, LND_time, passengerMass, passengerNames\
, passengerPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2, ACC_Trim,\
 El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, Cg_shift, eigenmotions \
 = importExcelData(f)

alpha = []
delta = []
for i in range(len(El_Trim_Curve)):
    alpha.append(El_Trim_Curve[i][4])
    delta.append(El_Trim_Curve[i][5])

plt.scatter(alpha, delta)
plt.show()
