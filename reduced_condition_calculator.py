from Cit_par import *

import numpy as np


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
        a = np.sqrt((self.gamma * self.R * temp))
        print "the speed of sound equals : " + str(a)

        return mach * a

    def calc_V_e(self, V_t):
        density = self.calc_density()

        return V_t * np.sqrt((density/self.rho0))

    def calc_V_e_reduced(self, V_e, W):
        return V_e * np.sqrt((self.W_S/W))


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




