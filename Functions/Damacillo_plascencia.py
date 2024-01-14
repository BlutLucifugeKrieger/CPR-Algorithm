import math
from math import sqrt, log10
from scipy.stats import qmc
import numpy as np
#variables#
fseed = 0.001
#variables#

#lists#
laminar_damacillo_plascencia = []
turbulent_damacillo_plascencia = []
damacillo_plascencia_list = []
relative_error_damacillo_plascencia= []
vector_message_regime_damacillo_plascencia = []
message_vector_length_damacillo_plascencia = []
r_min_value_damacillo_plascencia = []
e_min_value_damacillo_plascencia = []
r_max_value_damacillo_plascencia = []
e_max_value_damacillo_plascencia = []
#lists#

class Damacillo_plasc:

    def __init__(self, i, j, x, ii, jj, s, o):

        self.lower_bounds = [i, ii]
        self.upper_bounds = [j, jj]
        self.sampler = qmc.Sobol(d=2, scramble=s, optimization=o)
        self.sample = self.sampler.random_base2(m=x)
        self.reynolds_and_relative_roughness_list = qmc.scale(self.sample, self.lower_bounds, self.upper_bounds)
        self.reynolds_list_damacillo_plascencia = [sub_list[0] for sub_list in self.reynolds_and_relative_roughness_list]
        self.relative_roughness_list_damacillo_plascencia = [sub_list[1] for sub_list in self.reynolds_and_relative_roughness_list]

    def message(self,r):
        if (r <= 2000):
            vector_message_regime_damacillo_plascencia.append("laminar")
        elif (r >= 4000 and r <= 10 ** 8):
            vector_message_regime_damacillo_plascencia.append("turbulent")

    def colebrook_equation(self,r, e, x0):
        g = -2 * log10((e / 3.7) + 2.51 * x0 / r)
        return g

    def colebrook_derivative_equation(self,r, e, x0):
        g_derivative = (-2 / math.log(10)) * (2.51 / (r * ((e / 3.7) + 2.51 * x0 / r)))
        return g_derivative

    def damacillo_plascencia_equation(self,r,e):
        lambda1 = 0.02
        tau1 = 3000
        lambda2 = abs(lambda1-(1/(-2*log10(e/3.71)))**2)
        tau2 = 0.77505/e**2 - 10.984/e + 7953.8
        f1 = 64/r + lambda1/(1+np.exp((tau1-r)/100)) + lambda2/(1+np.exp(((tau2-r)/150)*e))
        return f1

    def CRP(self):
        x0 = 1 / sqrt(fseed)
        for e in self.relative_roughness_list_damacillo_plascencia:
            for r in self.reynolds_list_damacillo_plascencia:
                iteration = 0
                max_iteration = 7

                if (r <= 2000):
                    f = 64 / r
                    laminar_damacillo_plascencia.append(f)

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
                    turbulent_damacillo_plascencia.append(f)

                c = self.damacillo_plascencia_equation(r, e)
                damacillo_plascencia_list.append(c)
                self.message(r)

    def laminar_relative_error(self):
        for x, y in zip(laminar_damacillo_plascencia, damacillo_plascencia_list):
            relative_error_damacillo_plascencia.append((abs(x - y) / x) * 100)

    def turbulent_relative_error(self):
        for x, y in zip(turbulent_damacillo_plascencia, damacillo_plascencia_list):
            relative_error_damacillo_plascencia.append((abs(x - y) / x) * 100)

    def min_relative_error_value(self):
        val = min(relative_error_damacillo_plascencia)
        return val

    def avg_relative_error(self):
        val = sum(relative_error_damacillo_plascencia) / len(relative_error_damacillo_plascencia)
        return val

    def max_relative_error_value(self):
        val = max(relative_error_damacillo_plascencia)
        return val

    def lenght_values(self):
        for r in self.reynolds_list_damacillo_plascencia:
            if (r <= 2000):
                return laminar_damacillo_plascencia
            elif (r >= 4000 and r <= 10 ** 8):
                return turbulent_damacillo_plascencia

