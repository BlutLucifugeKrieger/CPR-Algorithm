import math
from math import sqrt, log10

from scipy.stats import qmc

# variables#
fseed = 0.001
# variables#

# lists#
laminar_russian = []
turbulent_russian = []
russian_list = []
relative_error_russian= []
vector_message_regime_russian = []
message_vector_length_russian = []
r_min_value_russian = []
e_min_value_russian = []
r_max_value_russian = []
e_max_value_russian = []




class Russian:

    def __init__(self, i, j, x, ii, jj, s, o):

        self.lower_bounds = [i, ii]
        self.upper_bounds = [j, jj]
        self.sampler = qmc.Sobol(d=2, scramble=s, optimization=o)
        self.sample = self.sampler.random_base2(m=x)
        self.reynolds_and_relative_roughness_list = qmc.scale(self.sample, self.lower_bounds, self.upper_bounds)
        self.reynolds_list_russian = [sub_list[0] for sub_list in self.reynolds_and_relative_roughness_list]
        self.relative_roughness_list_russian = [sub_list[1] for sub_list in self.reynolds_and_relative_roughness_list]

    def message(self,r):
        if (r <= 2000):
            vector_message_regime_russian.append("laminar")
        elif (r >= 4000 and r <= 10 ** 8):
            vector_message_regime_russian.append("turbulent")

    def colebrook_equation(self,r, e, x0):
        g = -2 * log10((e / 3.7) + 2.51 * x0 / r)
        return g

    def colebrook_derivative_equation(self,r, e, x0):
        g_derivative = (-2 / math.log(10)) * (2.51 / (r * ((e / 3.7) + 2.51 * x0 / r)))
        return g_derivative
    def russian_equation(self,r, e):
        alfa = 68 / r
        x = (28 * alfa) ** 10
        f1 = 0.11 * ((alfa + e + x ** 1.4) / (115 * x + 1)) ** 0.25
        return f1

    def CRP(self):
        x0 = 1 / sqrt(fseed)
        for e in self.relative_roughness_list_russian:
            for r in self.reynolds_list_russian:
                iteration = 0
                max_iteration = 7

                if (r <= 2000):
                    f = 64 / r
                    laminar_russian.append(f)

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
                    turbulent_russian.append(f)

                c = self.russian_equation(r, e)
                russian_list.append(c)
                self.message(r)


    def laminar_relative_error(self):
        for x, y in zip(laminar_russian, russian_list):
            relative_error_russian.append((abs(x - y) / x) * 100)

    def turbulent_relative_error(self):
        for x, y in zip(turbulent_russian, russian_list):
            relative_error_russian.append((abs(x - y) / x) * 100)

    def min_relative_error_value(self):
        val = min(relative_error_russian)
        return val

    def avg_relative_error(self):
        val = sum(relative_error_russian) / len(relative_error_russian)
        return val

    def max_relative_error_value(self):
        val = max(relative_error_russian)
        return val

    def lenght_values(self):
        for r in self.reynolds_list_russian:
            if (r <= 2000):
                return laminar_russian
            elif (r >= 4000 and r <= 10 ** 8):
                return turbulent_russian
