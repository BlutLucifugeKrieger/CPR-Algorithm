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
# Algorithm made by CRP

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

            typed_re_lower_value = int(typed_re_lower)

            typed_re_upper = input("Enter the value for Reynolds's upper limit: ")

            typed_re_upper_value = int(typed_re_upper)

            typed_m_value = input("Enter the value for Sobol's potency : ")

            typed_m = int(typed_m_value)


            typed_r_lower_value = input("Enter the value for relative roughness lower limit: ")

            typed_rl = int(typed_r_lower_value)

            typed_r_upper_value = input("Enter the value for relative roughness upper limit: ")

            typed_ru = float(typed_r_upper_value)


            if typed_number_value == 1:



                churchill = Churchill()

                churchill.sobol(typed_re_lower_value,typed_re_upper_value,typed_m,typed_rl,typed_ru)
                churchill.CRP()
                churchill.laminar_relative_error()
                churchill.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print("Results")
                print()
                print("laminar flow list: ", laminar_churchill)
                print()
                print("turbulent flow list :", turbulent_churchill)
                print()
                print("Churchill list: ", churchill_list)
                print()
                print("relative error list (%): ", relative_error_churchill)
                print()
                print("number of indexes of the churchill's list: ", len(churchill_list))
                print("number of churchill's indexes and", vector_message_regime_churchill[0], ":",
                      len(relative_error_churchill))
                print()
                print(vector_message_regime_churchill[0], "list length: ", len(churchill.lenght_values()))
                print()
                print(vector_message_regime_churchill[0], "minimum relative error: ", churchill.min_relative_error_value()," %")
                print(vector_message_regime_churchill[0], "average relative error: ", churchill.avg_relative_error()," %")
                print(vector_message_regime_churchill[0], "maximum relative error: ", churchill.max_relative_error_value()," %")
                print()
                churchill.value_r_e_min_relative_error()
                churchill.value_r_e_max_relative_error()


            elif typed_number_value == 2:

                swamee = Swamee()

                swamee.sobol(typed_re_lower_value,typed_re_upper_value,typed_m,typed_rl,typed_ru)
                swamee.CRP()
                swamee.laminar_swamee_relative_error_swamee()
                swamee.turbulent_swamee_relative_error_swamee()
                print()
                print("------------------------------------------")
                print("Results")
                print()
                print("laminar flow list: ", laminar_swamee)
                print()
                print("turbulent flow list :", turbulent_swamee)
                print()
                print("Swamee list: ", swamee_list)
                print()
                print("relative error list (%): ", relative_error_swamee)
                print()
                print("number of indexes of the swamee's list: ", len(swamee_list))
                print("number of swamee's indexes and", vector_message_regime_swamee[0], ":",
                      len(relative_error_swamee))
                print()
                print(vector_message_regime_swamee[0], "list length: ", len(swamee.lenght_values()))
                print()
                print(vector_message_regime_swamee[0], "minimum relative error: ", swamee.min_relative_error_swamee_value()," %")
                print(vector_message_regime_swamee[0], "average relative error: ", swamee.avg_relative_error_swamee()," %")
                print(vector_message_regime_swamee[0], "maximum relative error: ", swamee.max_relative_error_swamee_value()," %")
                print()
                swamee.value_r_e_min_relative_error_swamee()
                swamee.value_r_e_max_relative_error_swamee()

            elif typed_number_value == 3:


                russian = Russian()

                russian.sobol(typed_re_lower_value,typed_re_upper_value,typed_m,typed_rl,typed_ru)
                russian.CRP()
                russian.laminar_relative_error_russian()
                russian.turbulent_relative_error_russian()
                print()
                print("------------------------------------------")
                print("Results")
                print()
                print("laminar flow list: ", laminar_russian)
                print()
                print("turbulent flow list :", turbulent_russian)
                print()
                print("Russian list: ", russian_list)
                print()
                print("relative error list (%): ", relative_error_russian)
                print()
                print("number of indexes of the russian's list: ", len(russian_list))
                print("number of russian's indexes and", vector_message_regime_russian[0], ":", len(relative_error_russian))
                print()
                print(vector_message_regime_russian[0], "list length: ", len(russian.lenght_values()))
                print()
                print(vector_message_regime_russian[0], "minimum relative error: ", russian.min_relative_error_russian_value()," %")
                print(vector_message_regime_russian[0], "average relative error: ", russian.avg_relative_error_russian()," %")
                print(vector_message_regime_russian[0], "maximum relative error: ", russian.max_relative_error_russian_value()," %")
                print()
                russian.value_r_e_min_relative_error_russian()
                russian.value_r_e_max_relative_error_russian()

            elif typed_number_value == 4:


                brkic_p = Brkic_praks()

                brkic_p.sobol(typed_re_lower_value,typed_re_upper_value,typed_m,typed_rl,typed_ru)
                brkic_p.CRP()
                brkic_p.laminar_brkic_p_relative_error_brkic_p()
                brkic_p.turbulent_brkic_p_relative_error_brkic_p()
                print()
                print("------------------------------------------")
                print("Results")
                print()
                print("laminar flow list: ", laminar_brkic_p)
                print()
                print("turbulent flow list :", turbulent_brkic_p)
                print()
                print("Brkic and Praks list: ", brkic_praks_list)
                print()
                print("relative error list (%): ", relative_error_brkic_p)
                print()
                print("number of indexes of the Brkic and Praks's list: ", len(brkic_praks_list))
                print("number of Brkic and Praks's indexes and", vector_message_regime_brkic_p[0], ":",
                      len(relative_error_brkic_p))
                print()
                print(vector_message_regime_brkic_p[0], "list length: ", len(brkic_p.lenght_values()))
                print()
                print(vector_message_regime_brkic_p[0], "minimum relative error: ", brkic_p.min_relative_error_brkic_p_value()," %")
                print(vector_message_regime_brkic_p[0], "average relative error: ", brkic_p.avg_relative_error_brkic_p()," %")
                print(vector_message_regime_brkic_p[0], "maximum relative error: ", brkic_p.max_relative_error_brkic_p_value()," %")
                print()
                brkic_p.value_r_e_min_relative_error_brkic_p()
                brkic_p.value_r_e_max_relative_error_brkic_p()

            elif typed_number_value == 5:


                damacillio = Damacillo_plasc()
                damacillio.sobol(typed_re_lower_value,typed_re_upper_value,typed_m,typed_rl,typed_ru)

                damacillio.CRP()
                damacillio.laminar_relative_error()
                damacillio.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print("Results")
                print()
                print("laminar flow list: ", laminar_damacillo_plascencia)
                print()
                print("turbulent flow list :", turbulent_damacillo_plascencia)
                print()
                print("damacillo plascencia list: ", damacillo_plascencia_list)
                print()
                print("relative error list (%): ", relative_error_damacillo_plascencia)
                print()
                print(vector_message_regime_damacillo_plascencia[0], "list length: ", len(damacillio.lenght_values()))
                print("number of indexes of the damacillo plascencia list: ", len(damacillo_plascencia_list))
                print("number of relative error indexes: ", len(relative_error_damacillo_plascencia))
                print()
                print(vector_message_regime_damacillo_plascencia[0], "minimum relative error: ",damacillio.min_relative_error_value(), "%")
                print(vector_message_regime_damacillo_plascencia[0], "average relative error: ", damacillio.avg_relative_error(),"%")
                print(vector_message_regime_damacillo_plascencia[0], "maximum  relative error: ",
                      damacillio.max_relative_error_value(), "%")
                print()
                damacillio.value_r_e_min_relative_error()
                damacillio.value_r_e_max_relative_error()

            elif typed_number_value == 6:

                avci = Avci_karagoz()
                avci.sobol(typed_re_lower_value,typed_re_upper_value,typed_m,typed_rl,typed_ru)

                avci.CRP()
                avci.laminar_relative_error()
                avci.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print("Results")
                print()
                print("laminar flow list: ",laminar_avci_karagoz)
                print()
                print("turbulent flow list :",turbulent_avci_karagoz)
                print()
                print("avci karagoz list: ", avci_karagoz_list)
                print()
                print("relative error list (%): ",relative_error_avci_karagoz)
                print()
                print(vector_message_regime_avci_karagoz[0],"list length: ",len(avci.lenght_values()))
                print("number of indexes of the avci karagoz list: ",len(avci_karagoz_list))
                print("number of relative error indexes: ",len(relative_error_avci_karagoz))
                print()
                print(vector_message_regime_avci_karagoz[0],"minimum relative error: ",avci.min_relative_error_value(),"%")
                print(vector_message_regime_avci_karagoz[0],"average relative error: ",avci.avg_relative_error(),"%")
                print(vector_message_regime_avci_karagoz[0],"maximum relative error: ",avci.max_relative_error_value(),"%")
                print()
                avci.value_r_e_min_relative_error()
                avci.value_r_e_max_relative_error()

            elif typed_number_value == 7:

                cheng = Cheng()
                cheng.sobol(typed_re_lower_value,typed_re_upper_value,typed_m,typed_rl,typed_ru)

                cheng.CRP()
                cheng.laminar_relative_error()
                cheng.turbulent_relative_error()

                print()
                print("------------------------------------------")
                print("Results")
                print()
                print("laminar flow list: ",laminar_cheng)
                print()
                print("turbulent flow list :",turbulent_cheng)
                print()
                print("cheng list: ", cheng_list)
                print()
                print("relative error list (%): ",relative_error_cheng)
                print()
                print(vector_message_regime_cheng[0],"list length: ",len(cheng.lenght_values()))
                print("number of indexes of the cheng list: ",len(cheng_list))
                print("number of relative error indexes: ",len(relative_error_cheng))
                print()
                print(vector_message_regime_cheng[0],"minimum relative error: ",cheng.min_relative_error_value(),"%")
                print(vector_message_regime_cheng[0],"average relative error: ",cheng.avg_relative_error(),"%")
                print(vector_message_regime_cheng[0],"maximum relative error: ",cheng.max_relative_error_value(),"%")
                print()
                cheng.value_r_e_min_relative_error()
                cheng.value_r_e_max_relative_error()

            elif typed_number_value == 8:

                modified = Modified_avci_karagoz()
                modified.sobol(typed_re_lower_value,typed_re_upper_value,typed_m,typed_rl,typed_ru)

                modified.CRP()
                modified.laminar_relative_error()
                modified.turbulent_relative_error()
                print()
                print("------------------------------------------")
                print("Results")
                print()
                print("laminar flow list: ",laminar_modified_avci_karagoz)
                print()
                print("turbulent flow list :",turbulent_modified_avci_karagoz)
                print()
                print("modified avci karagoz list: ", modified_avci_karagoz_list)
                print()
                print("relative error list (%): ",relative_error_modified_avci_karagoz)
                print()
                print(vector_message_regime_modified_avci_karagoz[0],"list length: ",len(modified.lenght_values()))
                print("number of indexes of the modified avci karagoz list: ",len(modified_avci_karagoz_list))
                print("number of relative error indexes: ",len(relative_error_modified_avci_karagoz))
                print()
                print(vector_message_regime_modified_avci_karagoz[0],"minimum relative error: ",modified.min_relative_error_value(),"%")
                print(vector_message_regime_modified_avci_karagoz[0],"average relative error: ",modified.avg_relative_error(),"%")
                print(vector_message_regime_modified_avci_karagoz[0],"maximum relative error: ",modified.max_relative_error_value(),"%")
                print()
                modified.value_r_e_min_relative_error()
                modified.value_r_e_max_relative_error()

            elif typed_number_value == 9:

                mbp = Mbpls()
                mbp.sobol(typed_re_lower_value,typed_re_upper_value,typed_m,typed_rl,typed_ru)

                mbp.CRP()
                mbp.laminar_relative_error()
                mbp.turbulent_relative_error()

                print()
                print("------------------------------------------")
                print("Results")
                print()
                print("laminar flow list: ",laminar_mbpls)
                print()
                print("turbulent flow list :",turbulent_mbpls)
                print()
                print("mbpls list: ", modified_mbpls)
                print()
                print("relative error list (%): ",relative_error_mbpls)
                print()
                print(vector_message_regime_mbpls[0],"list length: ",len(mbp.lenght_values()))
                print("number of indexes of the mbpls list: ",len(modified_mbpls))
                print("number of relative error indexes: ",len(relative_error_mbpls))
                print()
                print(vector_message_regime_mbpls[0],"minimum relative error: ",mbp.min_relative_error_value(),"%")
                print(vector_message_regime_mbpls[0],"average relative error: ",mbp.avg_relative_error(),"%")
                print(vector_message_regime_mbpls[0],"maximum relative error: ",mbp.max_relative_error_value(),"%")
                print()
                mbp.value_r_e_min_relative_error()
                mbp.value_r_e_max_relative_error()


    except ValueError:
        print("Invalid value, please try again")
        init()


main()
