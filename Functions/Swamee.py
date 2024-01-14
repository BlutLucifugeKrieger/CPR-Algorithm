import math
from math import sqrt, log10

from scipy.stats import qmc

# variables#
fseed = 0.001
# variables#

# lists#
laminar_swamee = []
turbulent_swamee = []
swamee_list = []
relative_error_swamee = []
vector_message_regime_swamee = []
message_vector_length_swamee = []
r_min_value_swamee = []
e_min_value_swamee = []
r_max_value_swamee = []
e_max_value_swamee = []


# lists#

class Swamee:

    def __init__(self, i, j, x, ii, jj, s, o):

        self.lower_bounds = [i, ii]
        self.upper_bounds = [j, jj]
        self.sampler = qmc.Sobol(d=2, scramble=s, optimization=o)
        self.sample = self.sampler.random_base2(m=x)
        self.reynolds_and_relative_roughness_list = qmc.scale(self.sample, self.lower_bounds, self.upper_bounds)
        self.reynolds_list_swamee = [sub_list[0] for sub_list in self.reynolds_and_relative_roughness_list]
        self.relative_roughness_list_swamee = [sub_list[1] for sub_list in self.reynolds_and_relative_roughness_list]

    # reynolds list#

    def message(self, r):
        if (r <= 2000):
            vector_message_regime_swamee.append("laminar")
        elif (r >= 4000 and r <= 10 ** 8):
            vector_message_regime_swamee.append("turbulent")

    def colebrook_equation(self, r, e, x0):
        g = -2 * log10((e / 3.7) + 2.51 * x0 / r)
        return g

    def colebrook_derivative_equation(self, r, e, x0):
        g_derivative = (-2 / math.log(10)) * (2.51 / (r * ((e / 3.7) + 2.51 * x0 / r)))
        return g_derivative

    def swamee_equation(self,r, e):
        f1 = ((64 / r) ** 8 + 9.5 * (math.log((e / 3.7) + 5.74 / (r ** 0.9)) - (2500 / r) ** 6) ** -16) ** 0.125
        return f1


    def CRP(self):
        x0 = 1 / sqrt(fseed)
        for e in self.relative_roughness_list_swamee:
            for r in self.reynolds_list_swamee:
                iteration = 0
                max_iteration = 7

                if (r <= 2000):
                    f = 64 / r
                    laminar_swamee.append(f)

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
                    turbulent_swamee.append(f)

                c = self.swamee_equation(r, e)
                swamee_list.append(c)
                self.message(r)


    def laminar_relative_error(self):
        for x, y in zip(laminar_swamee, swamee_list):
            relative_error_swamee.append((abs(x - y) / x) * 100)


    def turbulent_relative_error(self):
        for x, y in zip(turbulent_swamee, swamee_list):
            relative_error_swamee.append((abs(x - y) / x) * 100)


    def min_relative_error_value(self):
        val = min(relative_error_swamee)
        return val


    def avg_relative_error(self):
        val = sum(relative_error_swamee) / len(relative_error_swamee)
        return val


    def max_relative_error_value(self):
        val = max(relative_error_swamee)
        return val


    def lenght_values(self):
        for r in self.reynolds_list_swamee:
            if (r <= 2000):
                return laminar_swamee
            elif (r >= 4000 and r <= 10 ** 8):
                return turbulent_swamee


