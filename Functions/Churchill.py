import math
from math import sqrt, log10
import matplotlib.pyplot as plt
from scipy.stats import qmc

# variables#
fseed = 0.001
# variables#

# lists#
laminar_churchill = []
turbulent_churchill = []
churchill_list = []
relative_roughness_list_churchill = []
reynolds_list_churchill = []
relative_error_churchill = []
vector_message_regime_churchill = []
message_vector_length_churchill = []
r_min_value_churchill = []
e_min_value_churchill = []
r_max_value_churchill = []
e_max_value_churchill = []
y_axis_churchill = []


class Churchill:


    def sobol(self,i,j,x,ii,jj):

        lower_bounds = [i, ii]
        upper_bounds = [j, jj]
        sampler = qmc.Sobol(d=2, scramble=False)
        sample = sampler.random_base2(m=x)
        reynolds_and_relative_roughness_list = qmc.scale(sample, lower_bounds, upper_bounds)
        reynolds_list_churchill.extend([sub_list[0] for sub_list in reynolds_and_relative_roughness_list if sub_list[0] != 0 and sub_list[0] != 4000])
        relative_roughness_list_churchill.extend([sub_list[1] for sub_list in reynolds_and_relative_roughness_list if sub_list[1] != 0])

    def message(self, r):
        if (r <= 2000):
            vector_message_regime_churchill.append("laminar")
        elif (r >= 4000 and r <= 10 ** 8):
            vector_message_regime_churchill.append("turbulent")

    def colebrook_equation(self, r, e, x0):
        g = -2 * log10((e / 3.7) + 2.51 * x0 / r)
        return g

    def colebrook_derivative_equation(self, r, e, x0):
        g_derivative = (-2 / math.log(10)) * (2.51 / (r * ((e / 3.7) + 2.51 * x0 / r)))
        return g_derivative

    def churchill_equation(self, r, e):
        A = (2.457 * math.log(1 / ((7 / r) ** 0.9 + 0.27 * e))) ** 16
        B = (37530 / r) ** 16
        f1 = 8 * ((8 / r) ** 12 + (A + B) ** (-3 / 2)) ** (1 / 12)
        return f1

    def CRP(self):

        x0 = 1 / sqrt(fseed)
        for e in relative_roughness_list_churchill:
            for r in reynolds_list_churchill:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
                    f = 64 / r
                    laminar_churchill.append(f)

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
                    turbulent_churchill.append(f)

                if (r != 0):
                    c = self.churchill_equation(r, e)
                    churchill_list.append(c)
                    self.message(r)

    def laminar_relative_error(self):
        for x, y in zip(laminar_churchill, churchill_list):
            relative_error_churchill.append((abs(x - y) / x) * 100)

    def turbulent_relative_error(self):
        for x, y in zip(turbulent_churchill, churchill_list):
            relative_error_churchill.append((abs(x - y) / x) * 100)

    def min_relative_error_value(self):
        val = min(relative_error_churchill)
        return val

    def avg_relative_error(self):
        val = sum(relative_error_churchill) / len(relative_error_churchill)
        return val

    def max_relative_error_value(self):
        val = max(relative_error_churchill)
        return val

    def value_r_e_min_relative_error(self):
        x0 = 1 / sqrt(fseed)
        for e in relative_roughness_list_churchill:
            for r in reynolds_list_churchill:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
                    f = 64 / r
                    laminar_churchill.append(f)


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
                    turbulent_churchill.append(f)

                if (r != 0):
                    f1 = self.churchill_equation(r, e)
                    relative_error_churchill = (abs(f - f1) / f) * 100

                    if (relative_error_churchill == self.min_relative_error_value()):
                        r_min_value_churchill.append(r)
                        e_min_value_churchill.append(e)
                        return print("Relative roughness value: ", e, ",", "Reynolds value: ", r, ",",
                                     "minimum relative error: ", self.min_relative_error_value(), "%")

    def value_r_e_max_relative_error(self):
        x0 = 1 / sqrt(fseed)
        for e in relative_roughness_list_churchill:
            for r in reynolds_list_churchill:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
                    f = 64 / r
                    laminar_churchill.append(f)

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
                    turbulent_churchill.append(f)

                if (r != 0):
                    f1 = self.churchill_equation(r, e)
                    relative_error_churchill = (abs(f - f1) / f) * 100

                    if (relative_error_churchill == self.max_relative_error_value()):
                        r_max_value_churchill.append(r)
                        e_max_value_churchill.append(e)
                        return print("Relative roughness value: ", e, ",", "Reynolds value: ", r, ",",
                                     "maximum relative error: ", self.max_relative_error_value(), "%")

    def lenght_values(self):
        for r in reynolds_list_churchill:
            if (r <= 2000):
                return laminar_churchill
            elif (r >= 4000 and r <= 10 ** 8):
                return turbulent_churchill


