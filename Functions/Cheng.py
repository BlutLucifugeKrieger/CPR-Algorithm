import math
from math import sqrt, log10

from scipy.stats import qmc

#variables#
fseed = 0.001
#variables#

#lists#
laminar_cheng = []
turbulent_cheng = []
cheng_list = []
relative_error_cheng= []
vector_message_regime_cheng = []
message_vector_length_cheng = []
r_min_value_cheng = []
e_min_value_cheng = []
r_max_value_cheng = []
e_max_value_cheng = []
#lists#



class Cheng:

    def __init__(self, i, j, x, ii, jj, s, o):

        self.lower_bounds = [i, ii]
        self.upper_bounds = [j, jj]
        self.sampler = qmc.Sobol(d=2, scramble=s, optimization=o)
        self.sample = self.sampler.random_base2(m=x)
        self.reynolds_and_relative_roughness_list = qmc.scale(self.sample, self.lower_bounds, self.upper_bounds)
        self.reynolds_list_cheng = [sub_list[0] for sub_list in self.reynolds_and_relative_roughness_list]
        self.relative_roughness_list_cheng = [sub_list[1] for sub_list in self.reynolds_and_relative_roughness_list]

    def message(self,r):
        if (r <= 2000):
            vector_message_regime_cheng.append("laminar")
        elif (r >= 4000 and r <= 10 ** 8):
            vector_message_regime_cheng.append("turbulent")

    def colebrook_equation(self,r, e, x0):
        g = -2 * log10((e / 3.7) + 2.51 * x0 / r)
        return g

    def colebrook_derivative_equation(self,r, e, x0):
        g_derivative = (-2 / math.log(10)) * (2.51 / (r * ((e / 3.7) + 2.51 * x0 / r)))
        return g_derivative

    def cheng_equation(self,r,e):
        alfa = 1/(1 + (r/2720)**9)
        beta = 1/(1 + (r*e/320)**2)
        f1 = 1/(((r/64)**alfa)*((1.8*log10(r/6.8))**(2*beta*(1-alfa)))*((2*log10(3.7/e))**(2*(1-alfa)*(1-beta))))
        return f1

    def CRP(self):
        x0 = 1 / sqrt(fseed)
        for e in self.relative_roughness_list_cheng:
            for r in self.reynolds_list_cheng:
                iteration = 0
                max_iteration = 7

                if (r <= 2000):
                    f = 64 / r
                    laminar_cheng.append(f)

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
                    turbulent_cheng.append(f)

                c = self.cheng_equation(r, e)
                cheng_list.append(c)
                self.message(r)

    def laminar_relative_error(self):
        for x, y in zip(laminar_cheng, cheng_list):
            relative_error_cheng.append((abs(x - y) / x) * 100)

    def turbulent_relative_error(self):
        for x, y in zip(turbulent_cheng, cheng_list):
            relative_error_cheng.append((abs(x - y) / x) * 100)

    def min_relative_error_value(self):
        val = min(relative_error_cheng)
        return val

    def avg_relative_error(self):
        val = sum(relative_error_cheng) / len(relative_error_cheng)
        return val

    def max_relative_error_value(self):
        val = max(relative_error_cheng)
        return val

    def lenght_values(self):
        for r in self.reynolds_list_cheng:
            if (r <= 2000):
                return laminar_cheng
            elif (r >= 4000 and r <= 10 ** 8):
                return turbulent_cheng

