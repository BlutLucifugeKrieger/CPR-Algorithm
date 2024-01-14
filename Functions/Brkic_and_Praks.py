import math
from math import sqrt, log10

from scipy.stats import qmc

#variables#
fseed = 0.001
#variables#

#lists#
laminar_brkic_p = []
turbulent_brkic_p = []
brkic_praks_list = []
relative_error_brkic_p= []
vector_message_regime_brkic_p = []
message_vector_length_brkic_p = []
r_min_value_brkic_p = []
e_min_value_brkic_p = []
r_max_value_brkic_p = []
e_max_value_brkic_p = []


# lists#

class Brkic_praks:

    def __init__(self, i, j, x, ii, jj, s, o):

        self.lower_bounds = [i, ii]
        self.upper_bounds = [j, jj]
        self.sampler = qmc.Sobol(d=2, scramble=s, optimization=o)
        self.sample = self.sampler.random_base2(m=x)
        self.reynolds_and_relative_roughness_list = qmc.scale(self.sample, self.lower_bounds, self.upper_bounds)
        self.reynolds_list_brkic_p = [sub_list[0] for sub_list in self.reynolds_and_relative_roughness_list]
        self.relative_roughness_list_brkic_p = [sub_list[1] for sub_list in self.reynolds_and_relative_roughness_list]

    def message(self,r):
        if (r <= 2000):
            vector_message_regime_brkic_p.append("laminar")
        elif (r >= 4000 and r <= 10 ** 8):
            vector_message_regime_brkic_p.append("turbulent")

    def colebrook_equation(self,r, e, x0):
        g = -2 * log10((e / 3.7) + 2.51 * x0 / r)
        return g

    def colebrook_derivative_equation(self,r, e, x0):
        g_derivative = (-2 / math.log(10)) * (2.51 / (r * ((e / 3.7) + 2.51 * x0 / r)))
        return g_derivative

    def brkic_praks_list_equation(self,r, e):
        y1 = 1 - 1048 / ((4.489 / 10 ** 20) * (r ** 6) * (0.148 * r - 2.306 * r / (0.003133 * r + 9.646)) + 1050)
        y2 = 1.012 - 1 / (0.02521 * r * e + 2.202)
        y3 = 1 - 1 / (0.000389 * (r ** 2) * (e ** 2) + 0.0000239 * r + 1.61)
        f1 = (64 / r) * (1 - y1) + (0.316 / r ** 0.25) * (y1 - y3) + (0.25 / ((log10(e / 3.71)) ** 2)) * y2
        return f1

    def CRP(self):
        x0 = 1 / sqrt(fseed)
        for e in self.relative_roughness_list_brkic_p:
            for r in self.reynolds_list_brkic_p:
                iteration = 0
                max_iteration = 7

                if (r <= 2000):
                    f = 64 / r
                    laminar_brkic_p.append(f)

                elif (r >= 4000 and r <= 10 ** 8):
                    while iteration < max_iteration:
                        iteration = iteration + 1
                        g = self.colebrook_equation(r, e, x0)
                        g_derivative = self.colebrook_derivative_equation(r, e, x0)
                        x1 = x0 - (g - x0) / (g_derivative - 1)
                        if abs(x1 - x0) == 0:
                            break
                        else:
                            x0 = x1
                    f = 1 / x1 ** 2
                    turbulent_brkic_p.append(f)

                c = self.brkic_praks_list_equation(r, e)
                brkic_praks_list.append(c)
                self.message(r)


    def laminar_relative_error(self):
        for x, y in zip(laminar_brkic_p, brkic_praks_list):
            relative_error_brkic_p.append((abs(x - y) / x) * 100)

    def turbulent_relative_error(self):
        for x, y in zip(turbulent_brkic_p, brkic_praks_list):
            relative_error_brkic_p.append((abs(x - y) / x) * 100)

    def min_relative_error_value(self):
        val = min(relative_error_brkic_p)
        return val

    def avg_relative_error(self):
        val = sum(relative_error_brkic_p) / len(relative_error_brkic_p)
        return val

    def max_relative_error_value(self):
        val = max(relative_error_brkic_p)
        return val

    def lenght_values(self):
        for r in self.reynolds_list_brkic_p:
            if (r <= 2000):
                return laminar_brkic_p
            elif (r >= 4000 and r <= 10 ** 8):
                return turbulent_brkic_p



