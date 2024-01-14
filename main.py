from Functions.Cheng import *
from Functions.Churchill import *
from Functions.Damacillo_plascencia import *
from Functions.Mbpls import *
from Functions.Modified_avci_karagoz import *
from Functions.Swamee import *
from Functions.Russian import *
from Functions.Brkic_and_Praks import *
from Functions.Avci_karagoz import *


# Welcome !!! #
# Algorithm made by Camilo Castro and Juan Rincon

def main():
    init()


# Initialization of classes (Functions to evaluate)


def init():
    try:

        print()
        print("Functions to evaluate")
        print()
        print("1: Churchill function")
        print("2: Swamee function")
        print("3: Russian function")
        print("4: Brkic and Praks function")
        print("5: Damacillio and Placencia function")
        print("6: Avci Karagoz function")
        print("7: Cheng function")
        print("8: Modified Avci Karagoz function")
        print("9: Mbpls function")
        print()

        typed_value = input("Please enter a value from 1 to 9: ")
        typed_number_value = int(typed_value)

        if typed_number_value > 9 or typed_number_value == 0 or typed_number_value < 0:
            print(f"Invalid number (there's no function associated to number: {typed_number_value})")
        else:
            typed_re_lower = input("Enter the value for Reynolds's lower limit: ")

            typed_re_lower_value = float(typed_re_lower)

            typed_re_upper = input("Enter the value for Reynolds's upper limit: ")

            typed_re_upper_value = float(typed_re_upper)

            typed_m_value = input("Enter the value for Sobol's potency : ")

            typed_m = int(typed_m_value)

            typed_r_lower_value = input("Enter the value for relative roughness lower limit: ")

            typed_rl = float(typed_r_lower_value)

            typed_r_upper_value = input("Enter the value for relative roughness upper limit: ")

            typed_ru = float(typed_r_upper_value)

            typed_points = input("What kind of range you prefer for the iterations?, type same or variable:  ")
            points = iterations(typed_points)

            typed_optimization = input("Do you want to optimize the algorithm? type yes or no: ")
            optim = optimization(typed_optimization)

            if typed_number_value == 1:

                churchill = Churchill(typed_re_lower_value, typed_re_upper_value, typed_m, typed_rl, typed_ru, points, optim)

                churchill.CRP()
                churchill.laminar_relative_error()
                churchill.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print(f"Results {optimization_message(optim)}")
                print()
                print("laminar churchill flow list: ", laminar_churchill)
                print()
                print("turbulent churchill flow list :", turbulent_churchill)
                print()
                print("churchill list: ", churchill_list)
                print()
                print("relative error list (%): ", relative_error_churchill)
                print()
                print(vector_message_regime_churchill[0], "list length: ", len(churchill.lenght_values()))
                print("number of indexes of the churchill's list: ", len(churchill_list))
                print("number of relative error indexes: ", len(relative_error_churchill))
                print()
                print(vector_message_regime_churchill[0], "minimum error: ", churchill.min_relative_error_value(), "%")
                print(vector_message_regime_churchill[0], "average error: ", churchill.avg_relative_error(), "%")
                print(vector_message_regime_churchill[0], "maximum error: ", churchill.max_relative_error_value(), "%")




            elif typed_number_value == 2:

                swamee = Swamee(typed_re_lower_value, typed_re_upper_value, typed_m, typed_rl, typed_ru, points, optim)

                swamee.CRP()
                swamee.laminar_relative_error()
                swamee.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print(f"Results {optimization_message(optim)}")
                print()
                print("laminar Swamee flow list: ", laminar_swamee)
                print()
                print("turbulent Swamee flow list :", turbulent_swamee)
                print()
                print("Swamee list: ", swamee_list)
                print()
                print("relative error list (%): ", relative_error_swamee)
                print()
                print(vector_message_regime_swamee[0], "list length: ", len(swamee.lenght_values()))
                print("number of indexes of the Swamee's list: ", len(swamee_list))
                print("number of relative error indexes: ", len(relative_error_swamee))
                print()
                print(vector_message_regime_swamee[0], "minimum error: ", swamee.min_relative_error_value(), "%")
                print(vector_message_regime_swamee[0], "average error: ", swamee.avg_relative_error(), "%")
                print(vector_message_regime_swamee[0], "maximum error: ", swamee.max_relative_error_value(), "%")

            elif typed_number_value == 3:

                russian = Russian(typed_re_lower_value, typed_re_upper_value, typed_m, typed_rl, typed_ru, points, optim)

                russian.CRP()
                russian.laminar_relative_error()
                russian.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print(f"Results {optimization_message(optim)}")
                print()
                print("laminar russian flow list: ", laminar_russian)
                print()
                print("turbulent russian flow list :", turbulent_russian)
                print()
                print("russian list: ", russian_list)
                print()
                print("relative error list (%): ", relative_error_russian)
                print()
                print(vector_message_regime_russian[0], "list length: ", len(russian.lenght_values()))
                print("number of indexes of the russian's list: ", len(russian_list))
                print("number of relative error indexes: ", len(relative_error_russian))
                print()
                print(vector_message_regime_russian[0], "minimum error: ", russian.min_relative_error_value(), "%")
                print(vector_message_regime_russian[0], "average error: ", russian.avg_relative_error(), "%")
                print(vector_message_regime_russian[0], "maximum error: ", russian.max_relative_error_value(), "%")

            elif typed_number_value == 4:

                brkic_p = Brkic_praks(typed_re_lower_value, typed_re_upper_value, typed_m, typed_rl, typed_ru, points, optim)

                brkic_p.CRP()
                brkic_p.laminar_relative_error()
                brkic_p.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print(f"Results {optimization_message(optim)}")
                print()
                print("laminar Brkic and Praks flow list: ", laminar_brkic_p)
                print()
                print("turbulent Brkic and Praks flow list :", turbulent_brkic_p)
                print()
                print("Brkic and Praks list: ", brkic_praks_list)
                print()
                print("relative error list (%): ", relative_error_brkic_p)
                print()
                print(vector_message_regime_brkic_p[0], "list length: ", len(brkic_p.lenght_values()))
                print("number of indexes of the Brkic and Praks's list: ", len(brkic_praks_list))
                print("number of relative error indexes: ", len(relative_error_brkic_p))
                print()
                print(vector_message_regime_brkic_p[0], "minimum error: ", brkic_p.min_relative_error_value(), "%")
                print(vector_message_regime_brkic_p[0], "average error: ", brkic_p.avg_relative_error(), "%")
                print(vector_message_regime_brkic_p[0], "maximum error: ", brkic_p.max_relative_error_value(), "%")

            elif typed_number_value == 5:

                damacillo = Damacillo_plasc(typed_re_lower_value, typed_re_upper_value, typed_m, typed_rl, typed_ru, points, optim)

                damacillo.CRP()
                damacillo.laminar_relative_error()
                damacillo.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print(f"Results {optimization_message(optim)}")
                print()
                print("laminar Damacillo and Plascencia flow list: ", laminar_damacillo_plascencia)
                print()
                print("turbulent Damacillo and Plascencia flow list :", turbulent_damacillo_plascencia)
                print()
                print("Damacillo and Plascencia list: ", damacillo_plascencia_list)
                print()
                print("relative error list (%): ", relative_error_damacillo_plascencia)
                print()
                print(vector_message_regime_damacillo_plascencia[0], "list length: ", len(damacillo.lenght_values()))
                print("number of indexes of the Damacillo and Plascencia's list: ", len(damacillo_plascencia_list))
                print("number of relative error indexes: ", len(relative_error_damacillo_plascencia))
                print()
                print(vector_message_regime_damacillo_plascencia[0], "minimum error: ", damacillo.min_relative_error_value(), "%")
                print(vector_message_regime_damacillo_plascencia[0], "average error: ", damacillo.avg_relative_error(), "%")
                print(vector_message_regime_damacillo_plascencia[0], "maximum error: ", damacillo.max_relative_error_value(), "%")

            elif typed_number_value == 6:

                avci = Avci_karagoz(typed_re_lower_value, typed_re_upper_value, typed_m, typed_rl, typed_ru, points, optim)

                avci.CRP()
                avci.laminar_relative_error()
                avci.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print(f"Results {optimization_message(optim)}")
                print()
                print("laminar Avci and Karagoz flow list: ", laminar_avci_karagoz)
                print()
                print("turbulent Avci and Karagoz flow list :", turbulent_avci_karagoz)
                print()
                print("Avci and Karagoz list: ", avci_karagoz_list)
                print()
                print("relative error list (%): ", relative_error_avci_karagoz)
                print()
                print(vector_message_regime_avci_karagoz[0], "list length: ", len(avci.lenght_values()))
                print("number of indexes of the Avci and Karagoz's list: ", len(avci_karagoz_list))
                print("number of relative error indexes: ", len(relative_error_avci_karagoz))
                print()
                print(vector_message_regime_avci_karagoz[0], "minimum error: ", avci.min_relative_error_value(), "%")
                print(vector_message_regime_avci_karagoz[0], "average error: ", avci.avg_relative_error(), "%")
                print(vector_message_regime_avci_karagoz[0], "maximum error: ", avci.max_relative_error_value(), "%")

            elif typed_number_value == 7:

                cheng = Cheng(typed_re_lower_value, typed_re_upper_value, typed_m, typed_rl, typed_ru, points, optim)

                cheng.CRP()
                cheng.laminar_relative_error()
                cheng.turbulent_relative_error()

                print()
                print("------------------------------------------")
                print(f"Results {optimization_message(optim)}")
                print()
                print("laminar Cheng flow list: ", laminar_cheng)
                print()
                print("turbulent Cheng flow list :", turbulent_cheng)
                print()
                print("Cheng list: ", cheng_list)
                print()
                print("relative error list (%): ", relative_error_cheng)
                print()
                print(vector_message_regime_cheng[0], "list length: ", len(cheng.lenght_values()))
                print("number of indexes of the Cheng's list: ", len(cheng_list))
                print("number of relative error indexes: ", len(relative_error_cheng))
                print()
                print(vector_message_regime_cheng[0], "minimum error: ", cheng.min_relative_error_value(), "%")
                print(vector_message_regime_cheng[0], "average error: ", cheng.avg_relative_error(), "%")
                print(vector_message_regime_cheng[0], "maximum error: ", cheng.max_relative_error_value(), "%")

            elif typed_number_value == 8:

                modified = Modified_avci_karagoz(typed_re_lower_value, typed_re_upper_value, typed_m, typed_rl, typed_ru, points, optim)

                modified.CRP()
                modified.laminar_relative_error()
                modified.turbulent_relative_error()

                print()
                print("------------------------------------------")
                print(f"Results {optimization_message(optim)}")
                print()
                print("laminar modified Avci and Karagoz flow list: ", laminar_modified_avci_karagoz)
                print()
                print("turbulent modified Avci and Karagoz  flow list :", turbulent_modified_avci_karagoz)
                print()
                print("modified Avci and Karagoz  list: ", modified_avci_karagoz_list)
                print()
                print("relative error list (%): ", relative_error_modified_avci_karagoz)
                print()
                print(vector_message_regime_modified_avci_karagoz[0], "list length: ", len(modified.lenght_values()))
                print("number of indexes of the modified Avci and Karagoz 's list: ", len(modified_avci_karagoz_list))
                print("number of relative error indexes: ", len(relative_error_modified_avci_karagoz))
                print()
                print(vector_message_regime_modified_avci_karagoz[0], "minimum error: ", modified.min_relative_error_value(), "%")
                print(vector_message_regime_modified_avci_karagoz[0], "average error: ", modified.avg_relative_error(), "%")
                print(vector_message_regime_modified_avci_karagoz[0], "maximum error: ", modified.max_relative_error_value(), "%")

            elif typed_number_value == 9:

                mbp = Mbpls(typed_re_lower_value, typed_re_upper_value, typed_m, typed_rl, typed_ru, points, optim)

                mbp.CRP()
                mbp.laminar_relative_error()
                mbp.turbulent_relative_error()

                print()
                print("------------------------------------------")
                print(f"Results {optimization_message(optim)}")
                print()
                print("laminar Mbpls flow list: ", laminar_mbpls)
                print()
                print("turbulent Mbpls flow list :", turbulent_mbpls)
                print()
                print("Mbpls list: ", modified_mbpls_list)
                print()
                print("relative error list (%): ", relative_error_mbpls)
                print()
                print(vector_message_regime_mbpls[0], "list length: ", len(mbp.lenght_values()))
                print("number of indexes of the Mbpls's list: ", len(modified_mbpls_list))
                print("number of relative error indexes: ", len(relative_error_mbpls))
                print()
                print(vector_message_regime_mbpls[0], "minimum error: ", mbp.min_relative_error_value(), "%")
                print(vector_message_regime_mbpls[0], "average error: ", mbp.avg_relative_error(), "%")
                print(vector_message_regime_mbpls[0], "maximum error: ", mbp.max_relative_error_value(), "%")




    except ValueError:
        print("Invalid value, please try again")
        init()


def iterations(x):
    if x == "same":
        return False
    else:
        return True


def optimization(x):
    if x == "yes":
        val = input("Select a kind of optimization: type 1 for random-cd or type 2 for lloyd: ")

        if val == "1":
            return "random-cd"
        else:
            return "lloyd"

    elif x == "no":
        return None


def optimization_message(x):
    if x == "random-cd":
        return "optimized with random-cd"
    elif x == "lloyd":
        return "optimized with lloyd"
    elif x is None:
        return "unoptimized"


main()
