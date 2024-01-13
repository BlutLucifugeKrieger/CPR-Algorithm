import math
from math import sqrt, log10
from scipy.stats import qmc
import numpy as np
#variables#
fseed = 0.001
#variables#

#lists#
laminar_mbpls = []
turbulent_mbpls = []
modified_mbpls = []
relative_roughness_list_mbpls = []
reynolds_list_mbpls = []
relative_error_mbpls = []
vector_message_regime_mbpls = []
message_vector_length_mbpls = []
r_min_value_mbpls = []
e_min_value_mbpls = []
r_max_value_mbpls = []
e_max_value_mbpls = []
#lists#


class Mbpls:

    def sobol(self,i,j,x,ii,jj):

        lower_bounds = [i, ii]
        upper_bounds = [j, jj]
        sampler = qmc.Sobol(d=2, scramble=False)
        sample = sampler.random_base2(m=x)
        reynolds_and_relative_roughness_list = qmc.scale(sample, lower_bounds, upper_bounds)
        reynolds_list_mbpls.extend([sub_list[0] for sub_list in reynolds_and_relative_roughness_list if sub_list[0] != 0 and sub_list[0] != 4000])
        relative_roughness_list_mbpls.extend([sub_list[1] for sub_list in reynolds_and_relative_roughness_list if sub_list[1] != 0])

    def message(self,r):
        if(r <= 2000):
            vector_message_regime_mbpls.append("laminar")
        elif (r >= 4000 and r <= 10**8):
            vector_message_regime_mbpls.append("turbulent")

    def colebrook_equation(self,r,e,x0):
        g = -2*log10((e/3.7)+2.51*x0/r)
        return g

    def colebrook_derivative_equation(self,r,e,x0):
        g_derivative =(-2/math.log(10))*(2.51/(r*((e/3.7)+2.51*x0/r)))
        return g_derivative

    def mbpls_equation(self,r,e):
        f1 = 61.395/r + (0.024444+0.60915*e)/(np.exp(8188400/(r**2)))
        return f1

    def CRP(self):
        x0 = 1/sqrt(fseed)
        for e in relative_roughness_list_mbpls:
            for r in reynolds_list_mbpls:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
                    f = 64/r
                    laminar_mbpls.append(f)


                elif (r >= 4000 and r <= 10**8):
                    while iteration < max_iteration:
                        iteration = iteration + 1
                        g = self.colebrook_equation(r,e,x0)
                        g_derivative = self.colebrook_derivative_equation(r,e,x0)
                        x1 = x0 - (g-x0)/(g_derivative-1)
                        if abs(x1 - x0) == 0:
                            break
                        else:
                            x0 = x1
                    f = 1/x1**2
                    turbulent_mbpls.append(f)

                if (r != 0):
                    c = self.mbpls_equation(r,e)
                    modified_mbpls.append(c)
                    self.message(r)

    def laminar_relative_error(self):
        for x,y in zip(laminar_mbpls,modified_mbpls):
            relative_error_mbpls.append((abs(x-y)/x)*100)

    def turbulent_relative_error(self):
        for x,y in zip(turbulent_mbpls,modified_mbpls):
            relative_error_mbpls.append((abs(x-y)/x)*100)

    def min_relative_error_value(self):
        val = min(relative_error_mbpls)
        return val

    def avg_relative_error(self):
        val = sum(relative_error_mbpls)/len(relative_error_mbpls)
        return val

    def max_relative_error_value(self):
        val = max(relative_error_mbpls)
        return val

    def value_r_e_min_relative_error(self):
        x0 = 1/sqrt(fseed)
        for e in relative_roughness_list_mbpls:
            for r in reynolds_list_mbpls:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
                    f = 64/r
                    laminar_mbpls.append(f)

                elif (r >= 4000 and r <= 10 ** 8):
                    while iteration < max_iteration:
                        iteration = iteration + 1
                        g = self.colebrook_equation(r,e,x0)
                        g_derivative = self.colebrook_derivative_equation(r,e,x0)
                        x1 = x0 - (g - x0) / (g_derivative - 1)
                        if abs(x1 - x0) == 0:
                            break
                        else:
                            x0 = x1
                    f = 1/x1**2
                    turbulent_mbpls.append(f)

                if (r != 0):
                    f1 = self.mbpls_equation(r,e)
                    relative_error_mbpls = (abs(f-f1)/f)*100

                    if (relative_error_mbpls == self.min_relative_error_value()):
                        r_min_value_mbpls.append(r)
                        e_min_value_mbpls.append(e)
                        return print("Relative roughness value: ",e,",","Reynolds value: ",r,",","minimum relative error: ",self.min_relative_error_value(),"%")

    def value_r_e_max_relative_error(self):
        x0 = 1/sqrt(fseed)
        for e in relative_roughness_list_mbpls:
            for r in reynolds_list_mbpls:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
                    f = 64/r
                    laminar_mbpls.append(f)

                elif (r >= 4000 and r <= 10 ** 8):
                    while iteration < max_iteration:
                        iteration = iteration + 1
                        g = self.colebrook_equation(r,e,x0)
                        g_derivative = self.colebrook_derivative_equation(r,e,x0)
                        x1 = x0 - (g - x0) / (g_derivative - 1)
                        if abs(x1 - x0) == 0:
                            break
                        else:
                            x0 = x1
                    f = 1/x1**2
                    turbulent_mbpls.append(f)

                if (r != 0):
                    f1 = self.mbpls_equation(r,e)
                    relative_error_mbpls = (abs(f - f1) / f) * 100

                    if (relative_error_mbpls == self.max_relative_error_value()):
                        r_max_value_mbpls.append(r)
                        e_max_value_mbpls.append(e)
                        return print("Relative roughness value: ",e,",","Reynolds value: ",r,",","maximum relative error: ",self.max_relative_error_value(),"%")

    def lenght_values(self):
        for r in reynolds_list_mbpls:
            if (r <= 2000):
                return laminar_mbpls
            elif (r >= 4000 and r <= 10**8):
                return turbulent_mbpls


