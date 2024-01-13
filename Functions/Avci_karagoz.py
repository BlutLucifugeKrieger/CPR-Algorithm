import math
from math import sqrt, log10

from scipy.stats import qmc

#variables#
fseed = 0.001
#variables#

#lists#
laminar_avci_karagoz = []
turbulent_avci_karagoz = []
avci_karagoz_list = []
relative_roughness_list_avci_karagoz = []
reynolds_list_avci_karagoz = []
relative_error_avci_karagoz= []
vector_message_regime_avci_karagoz = []
message_vector_length_avci_karagoz = []
r_min_value_avci_karagoz = []
e_min_value_avci_karagoz = []
r_max_value_avci_karagoz = []
e_max_value_avci_karagoz = []
#lists#


class Avci_karagoz:
    def sobol(self, i, j, x,ii,jj):

        lower_bounds = [i, ii]
        upper_bounds = [j, jj]
        sampler = qmc.Sobol(d=2, scramble=False)
        sample = sampler.random_base2(m=x)
        reynolds_and_relative_roughness_list = qmc.scale(sample, lower_bounds, upper_bounds)
        reynolds_list_avci_karagoz.extend([sub_list[0] for sub_list in reynolds_and_relative_roughness_list if sub_list[0] != 0 and sub_list[0] != 4000])
        relative_roughness_list_avci_karagoz.extend([sub_list[1] for sub_list in reynolds_and_relative_roughness_list if sub_list[1] != 0])
    def message(self,r):
        if(r <= 2000):
            vector_message_regime_avci_karagoz.append("laminar")
        elif (r >= 4000 and r <= 10**8):
            vector_message_regime_avci_karagoz.append("turbulent")

    def colebrook_equation(self,r,e,x0):
        g = -2*log10((e/3.7)+2.51*x0/r)
        return g

    def colebrook_derivative_equation(self,r,e,x0):
        g_derivative =(-2/math.log(10))*(2.51/(r*((e/3.7)+2.51*x0/r)))
        return g_derivative

    def avci_karagoz_equation(self,r,e):
        cm = 1 + e + (e*sqrt(e))/(1+(225*e**3)) + 500*e**4
        ft = 6.4/abs((math.log((1/r) + 0.01*e*(1 + ((10*sqrt(e))/(1+(225*(e**2)))) + 5000*(e**3))))**2.4)
        f1 = ft + ((64/r) - ft)*math.exp(-((cm*r)/2560)**8)
        return f1

    def CRP(self):
        x0 = 1/sqrt(fseed)
        for e in relative_roughness_list_avci_karagoz:
            for r in reynolds_list_avci_karagoz:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
                    f = 64/r
                    laminar_avci_karagoz.append(f)

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
                    turbulent_avci_karagoz.append(f)

                if (r != 0):
                    c = self.avci_karagoz_equation(r,e)
                    avci_karagoz_list.append(c)
                    self.message(r)

    def laminar_relative_error(self):
        for x,y in zip(laminar_avci_karagoz,avci_karagoz_list):
            relative_error_avci_karagoz.append((abs(x-y)/x)*100)

    def turbulent_relative_error(self):
        for x,y in zip(turbulent_avci_karagoz,avci_karagoz_list):
            relative_error_avci_karagoz.append((abs(x-y)/x)*100)

    def min_relative_error_value(self):
        val = min(relative_error_avci_karagoz)
        return val

    def avg_relative_error(self):
        val = sum(relative_error_avci_karagoz)/len(relative_error_avci_karagoz)
        return val

    def max_relative_error_value(self):
        val = max(relative_error_avci_karagoz)
        return val

    def value_r_e_min_relative_error(self):
        x0 = 1/sqrt(fseed)
        for e in relative_roughness_list_avci_karagoz:
            for r in reynolds_list_avci_karagoz:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
                    f = 64/r
                    laminar_avci_karagoz.append(f)


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
                    turbulent_avci_karagoz.append(f)

                if (r != 0):
                    f1 = self.avci_karagoz_equation(r,e)
                    relative_error_avci_karagoz = (abs(f-f1)/f)*100

                    if (relative_error_avci_karagoz == self.min_relative_error_value()):
                        r_min_value_avci_karagoz.append(r)
                        e_min_value_avci_karagoz.append(e)
                        return print("Relative roughness value: ",e,",","Reynolds value: ",r,",","minimum relative error: ",self.min_relative_error_value(),"%")

    def value_r_e_max_relative_error(self):
        x0 = 1/sqrt(fseed)
        for e in relative_roughness_list_avci_karagoz:
            for r in reynolds_list_avci_karagoz:
                iteration = 0
                max_iteration = 7

                if (r <= 2000 and r != 0):
                    f = 64/r
                    laminar_avci_karagoz.append(f)

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
                    turbulent_avci_karagoz.append(f)

                if (r != 0):
                    f1 = self.avci_karagoz_equation(r,e)
                    relative_error_avci_karagoz = (abs(f - f1) / f) * 100

                    if (relative_error_avci_karagoz == self.max_relative_error_value()):
                        r_max_value_avci_karagoz.append(r)
                        e_max_value_avci_karagoz.append(e)
                        return print("Relative roughness value: ",e,",","Reynolds value: ",r,",","maximum relative error: ",self.max_relative_error_value(),"%")

    def lenght_values(self):
        for r in reynolds_list_avci_karagoz:
            if (r <= 2000):
                return laminar_avci_karagoz
            elif (r >= 4000 and r <= 10**8):
                return turbulent_avci_karagoz



