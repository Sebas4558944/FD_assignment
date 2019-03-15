from Cit_par import *
from reduced_condition_calculator import Conditions
from datareader import importExcelData

import numpy as np
import matplotlib.pyplot as plt
import pylab

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


def calc_de_dalpha(trim_curve):
    alpha = []
    delta = []

    for i in range(len(trim_curve)):
        alpha.append(trim_curve[i][4])
        delta.append(trim_curve[i][5])

    plt.scatter(alpha, delta)

    z = np.polyfit(alpha, delta, 1)
    p = np.poly1d(z)
    der = p.deriv()

    # print "the derivative :" + str(der)

    pylab.plot(alpha, p(alpha), "b")
    plt.show()

    return der


f = 'Reference_Datasheet.csv'
date_of_flight, flight_number, TO_time, LND_time, passengerMass, passengerNames \
    , passengerPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2, ACC_Trim, \
El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, Cg_shift, eigenmotions \
    = importExcelData(f)

de_dalpha = calc_de_dalpha(El_Trim_Curve)
print Cg_shift
