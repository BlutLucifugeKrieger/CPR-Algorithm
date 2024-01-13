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
relative_roughness_list_brkic_p = []
reynolds_list_brkic_p = []
relative_error_brkic_p= []
vector_message_regime_brkic_p = []
message_vector_length_brkic_p = []
r_min_value_brkic_p = []
e_min_value_brkic_p = []
r_max_value_brkic_p = []
e_max_value_brkic_p = []


# lists#

class Brkic_praks:

    def sobol(self,i,j,x,ii,jj):

        lower_bounds = [i, ii]
        upper_bounds = [j, jj]
        sampler = qmc.Sobol(d=2, scramble=False)
        sample = sampler.random_base2(m=x)
        reynolds_and_relative_roughness_list = qmc.scale(sample, lower_bounds, upper_bounds)
        reynolds_list_brkic_p.extend([sub_list[0] for sub_list in reynolds_and_relative_roughness_list if sub_list[0] != 0 and sub_list[0] != 4000])
        relative_roughness_list_brkic_p.extend([sub_list[1] for sub_list in reynolds_and_relative_roughness_list if sub_list[1] != 0])

    def message(self,r):
        if (r <= 2300):
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
        for e in relative_roughness_list_brkic_p:
            for r in reynolds_list_brkic_p:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
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

                if (r != 0):
                    c = self.brkic_praks_list_equation(r, e)
                    brkic_praks_list.append(c)
                    self.message(r)

    def laminar_brkic_p_relative_error_brkic_p(self):
        for x, y in zip(laminar_brkic_p, brkic_praks_list):
            relative_error_brkic_p.append((abs(x - y) / x) * 100)

    def turbulent_brkic_p_relative_error_brkic_p(self):
        for x, y in zip(turbulent_brkic_p, brkic_praks_list):
            relative_error_brkic_p.append((abs(x - y) / x) * 100)

    def min_relative_error_brkic_p_value(self):
        val = min(relative_error_brkic_p)
        return val

    def avg_relative_error_brkic_p(self):
        val = sum(relative_error_brkic_p) / len(relative_error_brkic_p)
        return val

    def max_relative_error_brkic_p_value(self):
        val = max(relative_error_brkic_p)
        return val

    def value_r_e_min_relative_error_brkic_p(self):
        x0 = 1 / sqrt(fseed)
        for e in relative_roughness_list_brkic_p:
            for r in reynolds_list_brkic_p:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
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

                if (r != 0):
                    f1 = self.brkic_praks_list_equation(r, e)
                    relative_error_brkic_p = (abs(f - f1) / f) * 100

                    if (relative_error_brkic_p == self.min_relative_error_brkic_p_value()):
                        r_min_value_brkic_p.append(r)
                        e_min_value_brkic_p.append(e)
                        return print("Relative roughness value: ", e, ",", "Reynolds value: ", r, ",",
                                     "minimum relative error: ", self.min_relative_error_brkic_p_value(), "%")

    def value_r_e_max_relative_error_brkic_p(self):
        x0 = 1 / sqrt(fseed)
        for e in relative_roughness_list_brkic_p:
            for r in reynolds_list_brkic_p:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
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

                if (r != 0):
                    f1 = self.brkic_praks_list_equation(r, e)
                    relative_error_brkic_p = (abs(f - f1) / f) * 100

                    if (relative_error_brkic_p == self.max_relative_error_brkic_p_value()):
                        r_max_value_brkic_p.append(r)
                        e_max_value_brkic_p.append(e)
                        return print("Relative roughness value: ", e, ",", "Reynolds value: ", r, ",",
                                     "maximum relative error: ", self.max_relative_error_brkic_p_value(), "%")

    def lenght_values(self):
        for r in reynolds_list_brkic_p:
            if (r <= 2000):
                return laminar_brkic_p
            elif (r >= 4000 and r <= 10 ** 8):
                return turbulent_brkic_p



