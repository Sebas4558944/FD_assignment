from Cit_par import *

import numpy as np
from conversion_helpers import *


class Conditions:

    def __init__(self, h_p):
        # Set standard conditions
        self.altitude = h_p
        self.p0 = P0
        self.T0 = Temp0
        self.rho0 = rho0
        self.g0 = g
        self.lambdas = lambdas
        self.R = R
        self.gamma = gamma
        self.standard_weight = 60500.0

        # Set iterable variables
        self.pressure = 0.0
        self.density = 0.0
        self.mach = 0.0
        self.temperature = 0.0
        self.true_airspeed = 0.0
        self.equivalent_airspeed = 0.0
        self.reduced_airspeed = 0.0

    def calc_pressure(self):
        self.pressure = self.p0 * np.power((1 + (self.lambdas * self.altitude) / self.T0),
                                           -self.g0 / (self.lambdas * self.R))
        return self.pressure

    def calc_density(self):
        self.density = self.rho0 * np.power((1 + (self.lambdas * self.altitude / self.T0)),
                                            (-(self.g0 / (self.lambdas * self.R) + 1)))

        return self.density

    def calc_mach(self, v_c):
        self.calc_pressure()
        inner = np.power((1 + ((self.gamma - 1) / (2 * self.gamma)) * (self.rho0 / self.p0) * (v_c ** 2)),
                         (self.gamma / (self.gamma - 1))) - 1

        self.mach = np.sqrt(
            (2 / (self.gamma - 1)) * (
                        np.power((1 + (self.p0 / self.pressure) * inner), ((self.gamma - 1) / self.gamma)) - 1))

        return self.mach

    def calc_temperature(self, t_m):
        self.temperature = t_m / (1 + ((self.gamma - 1) / 2) * (self.mach ** 2))
        return self.temperature

    def calc_true_airspeed(self):
        a = np.sqrt((self.gamma * self.R * self.temperature))

        self.true_airspeed = self.mach * a
        return self.true_airspeed

    def calc_equivalent_airspeed(self):
        self.equivalent_airspeed = self.true_airspeed * np.sqrt((self.density / self.rho0))
        return self.equivalent_airspeed

    def calc_reduced_airspeed(self, w):
        self.reduced_airspeed = self.equivalent_airspeed * np.sqrt((self.standard_weight / w))
        return self.reduced_airspeed

# reducend_conditions = Conditions(5000)
# pressure = reducend_conditions.calc_pressure()
# print "the pressure equals : " + str(pressure)
#
# mach = reducend_conditions.calc_mach(120)
# print "the mach number equals : " + str(mach)
# 
# temp = reducend_conditions.calc_temperature(250, mach)
# print "the temperature equals : " + str(temp)
#
# true_speed = reducend_conditions.calc_V_t(temp, mach)
# print "the true airspeed equals : " + str(true_speed)
