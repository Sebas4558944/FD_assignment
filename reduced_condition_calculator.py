from Cit_par import *

import numpy as np
from conversion_helpers import *


class Conditions:

    def __init__(self, h_p):
        self.h_p = h_p
        self.p0 = P0
        self.T0 = Temp0
        self.lambdas = lambdas
        self.R = R
        self.rho0 = rho0
        self.g0 = g
        self.W_S = 60500
        self.gamma = gamma

    def calc_pressure(self):
        return self.p0 * np.power((1 + (self.lambdas * self.h_p) / self.T0), -self.g0 / (self.lambdas * self.R))

    def calc_density(self):
        return self.rho0 * np.power((1 + (self.lambdas * self.h_p / self.T0)),
                                    (-(self.g0 / (self.lambdas * self.R) + 1)))

    def calc_temperature(self, T_m, mach):
        return T_m / (1 + ((self.gamma - 1) / 2) * (mach ** 2))

    def calc_mach(self, V_c):
        inner = np.power((1 + ((self.gamma - 1) / (2 * self.gamma)) * (self.rho0 / self.p0) * (V_c ** 2)),
                         (self.gamma / (self.gamma - 1))) - 1
        p = self.calc_pressure()

        return np.sqrt(
            (2 / (self.gamma - 1)) * (np.power((1 + (self.p0 / p) * inner), ((self.gamma - 1) / self.gamma)) - 1))

    def calc_V_t(self, temp, mach):
        temp += celsius_to_kelvin
        a = np.sqrt((self.gamma * self.R * temp))

        return mach * a

    def calc_V_e(self, V_t):
        density = self.calc_density()

        return V_t * np.sqrt((density/self.rho0))

    def calc_V_e_reduced(self, V_e, W):
        return V_e * np.sqrt((self.W_S/W))

    def calc_final(self, V_c, T_m):
        mach = self.calc_mach(V_c)
        temp = self.calc_temperature(T_m, mach)

        V_t = self.calc_V_t(temp, mach)
        V_e = self.calc_V_e(V_t)
        V_e_reduced = self.calc_V_e_reduced(V_e)

        return V_e_reduced


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




