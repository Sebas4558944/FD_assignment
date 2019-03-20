from Cit_par import *
from reduced_condition_calculator import Conditions
from datareader import importExcelData, ThrustingAllDayEveryday
from conversion_helpers import *
from Run1 import weight, f_used_trim

import numpy as np
import matplotlib.pyplot as plt
import pylab


class Elevator(Conditions):

    def __init__(self, h_p, v_c, t_m, w):
        Conditions.__init__(self, h_p)
        # Aircraft parameters
        self.s = S
        self.chord = c
        self.c_m_tc = -0.0064
        self.mass_flow_standard = 0.048
        self.weight = w

        # Calculate conditions
        self.calc_pressure()
        self.calc_density()
        self.calc_mach(v_c)
        self.calc_temperature(t_m)
        self.calc_true_airspeed()
        self.calc_equivalent_airspeed()
        self.calc_reduced_airspeed(self.weight)

        # Initializing parameters to calculate
        self.c_n = 0.0
        self.c_m_delta = 0.0
        self.c_m_alpha = 0.0

    def get_reduced_equivalent_airspeed(self):
        return self.reduced_airspeed

    def calc_temperature_difference(self):
        return self.temperature - self.T0 + (self.lambdas * self.altitude)

    def calc_c_n(self):
        self.c_n = self.weight / (0.5 * self.density * (self.true_airspeed ** 2) * self.s)

        if self.c_n <= 0:
            print "Warning: The aircraft is unstable, c_n = " + str(self.c_n)

        return self.c_n

    def calc_d_c_m(self, x_cg_2, x_cg_1):
        return self.c_n * ((x_cg_2 - x_cg_1) / self.chord)

    def calc_c_m_delta(self, delta_difference, c_m_difference):
        self.c_m_delta = -(c_m_difference / delta_difference)

        if self.c_m_delta >= 0:
            print "Warning: The aircraft is unstable, cm_delta = " + str(self.c_m_delta)
        return self.c_m_delta

    def calc_delta_reduced(self, delta_measured, t_cs, t_c):
        delta_reduced = delta_measured - (1 / self.c_m_delta) * self.c_m_tc * (t_cs - t_c)

        return delta_reduced

    def calc_c_m_alpha(self, derivative):
        self.c_m_alpha = -1 * derivative * self.c_m_delta

        if self.c_m_alpha <= 0:
            print "Warning: The aircraft is unstable, cm_alpha = " + str(self.c_m_alpha)

        return self.c_m_alpha

    def calc_thrust_coefficients(self, mass_flow_l, mass_flow_r):
        delta_temperature = self.calc_temperature_difference()

        standard_input = [self.altitude, self.mach, delta_temperature, self.mass_flow_standard, self.mass_flow_standard]
        thrust_standard = sum(ThrustingAllDayEveryday(standard_input))

        t_cs = thrust_standard / (0.5 * self.density * (self.true_airspeed ** 2) * self.s)

        non_standard_input = [self.altitude, self.mach, delta_temperature, mass_flow_l, mass_flow_r]
        thrust_non_standard = sum(ThrustingAllDayEveryday(non_standard_input))
        t_c = thrust_non_standard / (0.5 * self.density * (self.true_airspeed ** 2) * self.s)

        return t_cs, t_c

    def calc_reduced_delta(self, delta_measured, mass_flow_l, mass_flow_r):
        thrust_cs, thrust_c = self.calc_thrust_coefficients(mass_flow_l, mass_flow_r)

        reduced_delta = delta_measured - (1 / self.c_m_delta) * self.c_m_tc * (thrust_cs - thrust_c)

        return reduced_delta

    def calc_reduced_stick_force(self, force):
        return force * (self.standard_weight / self.weight)


def calc_d_e_d_alpha(trim_curve_data):
    angle_of_attack = []
    elevator_angle = []

    for i in range(len(trim_curve_data)):
        angle_of_attack.append(trim_curve_data[i][4])
        elevator_angle.append(trim_curve_data[i][5])

    z = np.polyfit(angle_of_attack, elevator_angle, 1)
    p = np.poly1d(z)
    derivative = p.deriv()

    # print "the derivative :" + str(der)

    # pylab.plot(alpha, p(alpha), "b")
    # plt.show()

    if derivative < 0:
        print "Warning: The aircraft is unstable, d_e_d_alpha = " + str(derivative)

    return derivative


# ----------------------------------------------------------------------------------------------------------------------

f_1 = 'Reference_Datasheet.csv'
f_2 = 'Post_Flight_Datasheet_13_03_V2.csv'

date_of_flight, flight_number, TO_time, LND_time, passengerMass, passengerNames, passengerPos, blockfuel, ACC_CLCD, \
CL_CD_series1, CL_CD_series2, ACC_Trim, El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, Cg_shift, \
eigenmotions = importExcelData(f_1)

# ----------------------------------------------------------------------------------------------------------------------
cg_weight = weight(15000, [1053])
elevator_effectiveness = Elevator(Cg_shift[0][2] * ft_to_m, Cg_shift[0][3] * kts_to_ms,
                                  Cg_shift[0][11] + celsius_to_kelvin, cg_weight[0])

elevator_difference = Cg_shift[1][5] - Cg_shift[0][5]

c_n = elevator_effectiveness.calc_c_n()
difference_cm = elevator_effectiveness.calc_d_c_m(float(pos_shifted), float(newpos_shifted))
cm_delta = elevator_effectiveness.calc_c_m_delta(elevator_difference, difference_cm)

slope = calc_d_e_d_alpha(El_Trim_Curve)
cm_alpha = elevator_effectiveness.calc_c_m_alpha(slope)

# ----------------------------------------------------------------------------------------------------------------------

delta = []
speed = []
alpha = []
stick_force = []
weights = weight(15000, f_used_trim)

for j in range(len(El_Trim_Curve)):
    m_flow_l = El_Trim_Curve[j][8] * lbs_per_hour_to_kg_per_s
    m_flow_r = El_Trim_Curve[j][9] * lbs_per_hour_to_kg_per_s
    delta_r = elevator_effectiveness.calc_reduced_delta(El_Trim_Curve[j][5], m_flow_l, m_flow_r)

    reducing_elevator = Elevator(El_Trim_Curve[j][2] * ft_to_m, El_Trim_Curve[j][3] * kts_to_ms,
                                 El_Trim_Curve[j][11] + celsius_to_kelvin, 7500*g)  # weights[j])

    airspeed = reducing_elevator.get_reduced_equivalent_airspeed()
    stick_force_reduced = reducing_elevator.calc_reduced_stick_force(El_Trim_Curve[j][7])

    alpha.append(El_Trim_Curve[j][4])
    speed.append(airspeed)
    stick_force.append(stick_force_reduced)
    delta.append(delta_r)

# ----------------------------------------------------------------------------------------------------------------------

print "c_n equals : " + str(c_n)
print "cm_delta equals : " + str(cm_delta)
print "cm_alpha equals : " + str(cm_alpha)
print "The weights are : " + str(weights)
print "The airspeed is : " + str(speed)
print "The angle of attack is : " + str(alpha)
print "The elevator angle is : " + str(newpos_shifted)
print "The stick force is : " + str(pos_shifted)

z1 = np.polyfit(speed, delta, 2)
p1 = np.poly1d(z1)

z2 = np.polyfit(alpha, delta, 1)
p2 = np.poly1d(z2)

z3 = np.polyfit(speed, stick_force, 3)
p3 = np.poly1d(z3)

speed.sort()
alpha.sort()

plt.figure(1)
pylab.plot(speed, p1(speed), "b")
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
pylab.plot(speed, p3(speed), "b")
plt.gca().invert_yaxis()
plt.xlabel("Reduced velocity")
plt.ylabel("Reduced stick-force")
plt.title("Control-force curve")

plt.show()
