from Cit_par import *
from reduced_condition_calculator import Conditions
from datareader import importExcelData, ThrustingAllDayEveryday
from conversion_helpers import *

import numpy as np
import matplotlib.pyplot as plt
import pylab

Cm0 = 0.0297
CNa = CNwa + CNha


def calc_CN(weight, density, v, s):
    C_N = weight / (0.5 * density * (v ** 2) * s)

    if C_N < 0:
        print "Warning: The aircraft is unstable, C_N = " + str(C_N)

    return C_N


def calc_DCm(x_cg_2, x_cg_1, C_N, chord):
    return C_N * ((x_cg_2 - x_cg_1) / chord)


def calc_Cm_delta(diff_delta, diff_Cm):
    Cm_delta = -(diff_Cm / diff_delta)

    if Cm_delta > 0:
        print "Warning: The aircraft is unstable, Cm_delta = " + str(Cm_delta)
    return Cm_delta


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

    # plt.scatter(alpha, delta)

    z = np.polyfit(alpha, delta, 1)
    p = np.poly1d(z)
    der = p.deriv()

    # print "the derivative :" + str(der)

    # pylab.plot(alpha, p(alpha), "b")
    # plt.show()

    if der < 0:
        print "Warning: The aircraft is unstable, de_dalpha = " + str(der)

    return der


def calc_Cm_alpha(derivative, Cm_delta):
    Cm_alpha = -1 * derivative * Cm_delta

    if Cm_alpha < 0:
        print "Warning: The aircraft is unstable, Cm_alpha = " + str(Cm_alpha)

    return Cm_alpha


def calc_thrust_coefficient(thrust, s, v, density):
    return thrust / (0.5 * density * (v ** 2) * s)


def calc_reduced_stick_force(w, w_s, force):
    return force*(w_s/w)


f = 'Reference_Datasheet.csv'
f2 = 'Post_Flight_Datasheet_13_03_V2.csv'
date_of_flight, flight_number, TO_time, LND_time, passengerMass, passengerNames \
    , passengerPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2, ACC_Trim, \
El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, Cg_shift, eigenmotions \
    = importExcelData(f2)

altitude = []
elevator = []
velocity = []
temperature = []
massflows = []
for i in range(len(Cg_shift)):
    altitude.append(Cg_shift[i][2] * ft_to_m)
    elevator.append(Cg_shift[i][5])
    velocity.append(Cg_shift[i][3] * kts_to_ms)
    temperature.append(Cg_shift[i][11])
    massflows.append([Cg_shift[i][8] * lbs_per_hour_to_kg_per_s, Cg_shift[i][9] * lbs_per_hour_to_kg_per_s])

reduced_calculator = Conditions(altitude[0] * ft_to_m)

m_flow_standard = 0.048
C_mtc = -0.0064

dens = reduced_calculator.calc_density()
mach = reduced_calculator.calc_mach(velocity[0])
temp = reduced_calculator.calc_temperature(temperature[0], mach)
speed_true = reduced_calculator.calc_V_t(temp, mach)

thrust_standard = [altitude[0], mach, temp, m_flow_standard, m_flow_standard]
T_standard = sum(ThrustingAllDayEveryday(thrust_standard))
Tcs = calc_thrust_coefficient(T_standard, S, speed_true, dens)

thrust_non_standard = [altitude[0], mach, temp, massflows[0][0], massflows[0][1]]
T_non_standard = sum(ThrustingAllDayEveryday(thrust_non_standard))
Tc = calc_thrust_coefficient(T_non_standard, S, speed_true, dens)

delta_diff = elevator[1] - elevator[0]
CN = calc_CN(7500 * g, dens, velocity[0], S)
dcm = calc_DCm(7.5, 8, CN, c)
cm_delta = calc_Cm_delta(delta_diff, dcm)

# ------------------------------------------------------------------------------------------
delta_reduced = []
speed_reduced = []
alpha = []
stick_force = []
print El_Trim_Curve
for i in range(len(El_Trim_Curve)):
    elevator_reduced = Conditions(El_Trim_Curve[i][2] * ft_to_m)
    delta_e = El_Trim_Curve[i][5]

    dens = elevator_reduced.calc_density()
    mach = elevator_reduced.calc_mach(El_Trim_Curve[i][3] * kts_to_ms)
    temp = elevator_reduced.calc_temperature(El_Trim_Curve[i][11], mach)
    speed_true = elevator_reduced.calc_V_t(temp, mach)
    speed_equivalent = elevator_reduced.calc_final(El_Trim_Curve[i][3], El_Trim_Curve[i][11], 7500*g)

    reduced_stick_force = calc_reduced_stick_force(7500*g, 60500, El_Trim_Curve[i][7])

    m_flow_l = El_Trim_Curve[i][8] * lbs_per_hour_to_kg_per_s
    m_flow_r = El_Trim_Curve[i][9] * lbs_per_hour_to_kg_per_s

    thrust_standard = [El_Trim_Curve[i][2], mach, temp, m_flow_standard, m_flow_standard]
    T_standard = sum(ThrustingAllDayEveryday(thrust_standard))
    Tcs = calc_thrust_coefficient(T_standard, S, speed_true, dens)

    thrust_non_standard = [altitude[0], mach, temp, m_flow_l, m_flow_r]
    T_non_standard = sum(ThrustingAllDayEveryday(thrust_non_standard))
    Tc = calc_thrust_coefficient(T_non_standard, S, speed_true, dens)

    stick_force.append(reduced_stick_force)
    alpha.append(El_Trim_Curve[i][4])
    speed_reduced.append(speed_equivalent)
    delta_reduced.append((delta_e - (1 / cm_delta) * C_mtc * (Tcs - Tc)))

slope = calc_de_dalpha(El_Trim_Curve)
cm_alpha = calc_Cm_alpha(slope, cm_delta)

print "Flying at an altitude of : " + str(altitude[0])
print "Flying at a mach number of : " + str(mach)
print "Flying with a speed of : " + str(speed_true)

print "Having a standard thrust of : " + str(T_standard)
print "Thus Tcs equals : " + str(Tcs)

print "Having a non-standard thrust of : " + str(T_non_standard)
print "Thus Tc equals : " + str(Tc)

print
print
print

print "Difference in elevator deflection : " + str(delta_diff)
print "CN equals : " + str(CN)
print "DCm equals : " + str(dcm)
print "Cm_delta equals : " + str(cm_delta)
print "Cm_alpha equals : " + str(cm_alpha)

z1 = np.polyfit(speed_reduced, delta_reduced, 2)
p1 = np.poly1d(z1)

z2 = np.polyfit(alpha, delta_reduced, 1)
p2 = np.poly1d(z2)

z3 = np.polyfit(speed_reduced, stick_force, 2)
p3 = np.poly1d(z3)

speed_reduced.sort()
alpha.sort()

plt.figure(1)
pylab.plot(speed_reduced, p1(speed_reduced), "b")
plt.gca().invert_yaxis()
plt.xlabel("Reduced velocity")
plt.ylabel("Reduced elevator deflection")
plt.title("Elevator-trim curve")

plt.figure(2)
pylab.plot(alpha, p2(alpha), "b")
plt.gca().invert_yaxis()
plt.xlabel("Angle of attack")
plt.ylabel("Reduced elevator deflection")
plt.title("Angle plot")

plt.figure(3)
pylab.plot(speed_reduced, p3(speed_reduced), "b")
plt.gca().invert_yaxis()
plt.xlabel("Reduced velocity")
plt.ylabel("Reduced stick-force")
plt.title("Control-force curve")

plt.show()

