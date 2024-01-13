# CRP Algorithm
[![Python Version](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-311/) [![SciPy Version](https://img.shields.io/badge/scipy-1.11.4-brightgreen.svg)](https://pypi.org/project/scipy/1.11.4/) [![NumPy Version](https://img.shields.io/badge/numpy-1.26.0-blue.svg)](https://pypi.org/project/numpy/1.26.0/)
[![Matplotlib Version](https://img.shields.io/badge/matplotlib-3.8.2-blue.svg)](https://pypi.org/project/matplotlib/3.8.2/)


This repository contains an algorithm that is capable to calculate the relative error of a laminar flow and turbulent flow depends on a the roughness and reynolds values


## What should i know before running the program?
First of all, check your current python version to avoid issues while the program is running,
if your are wondering what kind of libraries you need to have to run this program, we listed all of this just right here!:

 * Python version -> 3.11 or 3.12
 * SciPy version -> 1.11.4
 * Numpy version -> 1.26.0
 * Matplotlib -> 3.8.2

## How to use?

Well, to execute this program you need to clone this repository and after that, you should to open this project using a Python IDE or you can execute this project using the command prompt

### First case -> using an IDE (Pycharm) 


![image_ssd](https://github.com/BlutLucifugeKrieger/CRP-Algorithm/assets/130005378/50832fb3-45b7-434a-8da9-1a64d7398c9f)


### Second case -> using a command prompt

![image2_x](https://github.com/BlutLucifugeKrieger/CRP-Algorithm/assets/130005378/914576ec-9521-49ec-834e-ea6502f059bc)



## Program execution

Using the Pycharm IDE or the command prompt we got this:

![image](https://github.com/BlutLucifugeKrieger/CRP-Algorithm/assets/130005378/8ad45689-7934-440e-958c-add91063486b)

After that, you can choose a function to evaluate using the keyboard:

![image](https://github.com/BlutLucifugeKrieger/CRP-Algorithm/assets/130005378/860dcc66-c830-432a-820d-9ec1c38a8cb0)

Now, you can enter the values to define the range of Reynolds, roughness, and the exponent for the quantity of samples of sobol.

### Example of using the Churchill function

**--> main.py executed**


  **-->** Functions to evaluate
  
  **-->** 1: Churchill function
  
  **-->** 2: Swamee function
  
  **-->** 3: Russian function
  
  **-->** 4: Brkic and Praks function
  
  **-->** 5: Damacillio and Placencia function
  
  **-->** 6: Avci Karagoz function
  
  **-->** 7: Cheng function
  
  **-->** 8: Modified Avci Karagoz function
  
  **-->** 9: Mbpls function
  
  **-->** Please enter a value from 1 to 9: 1
  
  **-->** Enter the value for Reynolds's lower limit: 0
  
  **-->** Enter the value for Reynolds's upper limit: 2000
  
  **-->** Enter the value for Sobol's potency : 3
  
  **-->** Enter the value for relative roughness lower limit: 0
  
  **-->** Enter the value for relative roughness upper limit: 0.05
  
  ------------------------------------------
  **-->** Results
  
  **-->** laminar flow list:  [0.064, 0.042666666666666665, 0.128, 0.08533333333333333, 0.036571428571428574, 0.0512, 0.256, 0.064, 0.042666666666666665, 0.128, 0.08533333333333333, 0.036571428571428574, 0.0512, 0.256, 0.064, 0.042666666666666665, 0.128, 0.08533333333333333, 0.036571428571428574, 0.0512, 0.256, 0.064, 0.042666666666666665, 0.128, 0.08533333333333333, 0.036571428571428574, 0.0512, 0.256, 0.064, 0.042666666666666665, 0.128, 0.08533333333333333, 0.036571428571428574, 0.0512, 0.256, 0.064, 0.042666666666666665, 0.128, 0.08533333333333333, 0.036571428571428574, 0.0512, 0.256, 0.064, 0.042666666666666665, 0.128, 0.08533333333333333, 0.036571428571428574, 0.0512, 0.256]
  
  **-->** turbulent flow list : []
  
  **-->** Churchill list:  [0.06400000000000129, 0.042666668520312735, 0.12800000000000003, 0.08533333333333336, 0.03657183701380633, 0.05120000000313824, 0.25600000000000006, 0.06400000000000129, 0.04266666852030939, 0.12800000000000003, 0.08533333333333336, 0.036571837001153525, 0.05120000000313824, 0.25600000000000006, 0.06400000000000129, 0.04266666852031389, 0.12800000000000003, 0.08533333333333336, 0.03657183701786132, 0.05120000000313824, 0.25600000000000006, 0.06400000000000129, 0.042666668520311535, 0.12800000000000003, 0.08533333333333336, 0.03657183700938751, 0.05120000000313824, 0.25600000000000006, 0.06400000000000129, 0.04266666852031417, 0.12800000000000003, 0.08533333333333336, 0.036571837018818984, 0.05120000000313824, 0.25600000000000006, 0.06400000000000129, 0.04266666852030526, 0.12800000000000003, 0.08533333333333336, 0.036571836984444314, 0.05120000000313824, 0.25600000000000006, 0.06400000000000129, 0.04266666852031345, 0.12800000000000003, 0.08533333333333336, 0.03657183701633656, 0.05120000000313824, 0.25600000000000006]
  
  **-->** relative error list (%):  [2.0166160408230382e-12, 4.3444829764954784e-06, 2.168404344971009e-14, 3.252606517456513e-14, 0.001116834626667916, 6.1293743518209265e-09, 2.168404344971009e-14, 2.0166160408230382e-12, 4.344475137713771e-06, 2.168404344971009e-14, 3.252606517456513e-14, 0.0011168000291613073, 6.1293743518209265e-09, 2.168404344971009e-14, 2.0166160408230382e-12, 4.344485676158888e-06, 2.168404344971009e-14, 3.252606517456513e-14, 0.0011168457145379045, 6.1293743518209265e-09, 2.168404344971009e-14, 2.0166160408230382e-12, 4.344480162990841e-06, 2.168404344971009e-14, 3.252606517456513e-14, 0.0011168225439584087, 6.1293743518209265e-09, 2.168404344971009e-14, 2.0166160408230382e-12, 4.344486342943224e-06, 2.168404344971009e-14, 3.252606517456513e-14, 0.0011168483331517806, 6.1293743518209265e-09, 2.168404344971009e-14, 2.0166160408230382e-12, 4.344465461209382e-06, 2.168404344971009e-14, 3.252606517456513e-14, 0.0011167543399141083, 6.1293743518209265e-09, 2.168404344971009e-14, 2.0166160408230382e-12, 4.344484651587835e-06, 2.168404344971009e-14, 3.252606517456513e-14, 0.0011168415452736868, 6.1293743518209265e-09, 2.168404344971009e-14]
  
  **-->** number of indexes of the churchill's list:  49
  **-->** number of churchill's indexes and laminar : 49
  
  **-->** laminar list length:  49
  
  **-->** laminar minimum relative error:  2.168404344971009e-14  %
  **-->** laminar average relative error:  0.0001601673757824948  %
  **-->** laminar maximum relative error:  0.0011168483331517806  %
  
  **-->** Relative roughness value:  0.025 , Reynolds value:  500.0 , minimum relative error:  2.168404344971009e-14 %
  **-->** Relative roughness value:  0.043750000000000004 , Reynolds value:  1750.0 , maximum relative error:  0.0011168483331517806 %
  
  Process finished with exit code 0










