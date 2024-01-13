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

![image_li1](https://github.com/BlutLucifugeKrieger/CRP-Algorithm/assets/130005378/b058a313-6fd5-4d70-904e-8257b241e5d7)


After that, you can choose a function to evaluate using the keyboard:

![image_li2](https://github.com/BlutLucifugeKrieger/CRP-Algorithm/assets/130005378/11ecbec5-0344-44b8-85c0-3ed9f07db8ab)


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




## Analysis

**As you can see, the result reflects a laminar flow behavior because of the range selected for the Reynolds values.
Otherwise, if we select a 4000 or more as as a lower limit of Reynolds and 10000 or more as a upper limit, the result will show us a turbulet behavior.**

![image](https://github.com/BlutLucifugeKrieger/CRP-Algorithm/assets/130005378/f37351a9-c9ba-49b3-a937-3fc529df0db7)


**-> main.py executed**

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
**-->** Enter the value for Reynolds's lower limit: 4000
**-->** Enter the value for Reynolds's upper limit: 10000
**-->** Enter the value for Sobol's potency : 3
**-->** Enter the value for relative roughness lower limit: 0
**-->** Enter the value for relative roughness upper limit: 0.05

------------------------------------------
**-->** Results

**-->** laminar flow list:  []

**-->** turbulent flow list : [0.05759092262863681, 0.05684275233413723, 0.058713777379475404, 0.05808977339718978, 0.056554897679344664, 0.05718284955980753, 0.05951737118942298, 0.0473123724124838, 0.046317287036489074, 0.04877311574964676, 0.047965861874627175, 0.04592921590433478, 0.04677196726829221, 0.04979757523310998, 0.06653881957100381, 0.0659160066009069, 0.0674826312771987, 0.06695681706460098, 0.06567773751883775, 0.0661984982819381, 0.06816440042038048, 0.05267767778374759, 0.05183053663883548, 0.053938374682274796, 0.05323927335902978, 0.051502945924172046, 0.052216374447379595, 0.05483353133550177, 0.07070910779898001, 0.07012883325651365, 0.07159108436511338, 0.07109933867075391, 0.06990721878388641, 0.07039185446527378, 0.07223005357963169, 0.04125727864805727, 0.040013745696316404, 0.04303758342314546, 0.04205976575334951, 0.03952082328880941, 0.04058540743217594, 0.04425966783080628, 0.0621841405319751, 0.06150727054049174, 0.06320575402677074, 0.06263718647953287, 0.06124771514479607, 0.061814560038628126, 0.06394085275087075]

**-->** Churchill list:  [0.05900128439869199, 0.058084292790041625, 0.06036305199779492, 0.05960840638508816, 0.05772894139088041, 0.058502202080809994, 0.06131339696086385, 0.048505927589905615, 0.047385840212644546, 0.05014662213546404, 0.04924020780926434, 0.04694788076816056, 0.04789805201280897, 0.0512939710429984, 0.06805113971307533, 0.06723373192071717, 0.069268265839006, 0.06859419992307979, 0.06691809183503435, 0.06760575841704701, 0.07007890919806403, 0.05400486949375304, 0.05300754888531328, 0.05547881910475935, 0.05466279361206875, 0.05261979639139028, 0.05346263812822697, 0.05651291974724113, 0.07225682418049642, 0.07147246519444654, 0.07342273827629772, 0.07277819311334355, 0.07116988953706362, 0.07182932276700377, 0.07415835178000532, 0.04219578966999687, 0.040859996275441936, 0.044117077708924105, 0.043060378381364076, 0.04033152088034523, 0.041473545175571135, 0.0454422386824249, 0.06365251999380181, 0.06279250984151356, 0.06493265290518185, 0.06422316494688089, 0.06245994207807107, 0.06318413596273686, 0.06581203464516983]

**-->** relative error list (%):  [2.4489306746301063, 2.1841666789923946, 2.8090078545960666, 2.614286300472722, 2.0759364081822715, 2.3072521414354417, 3.0176496971359072, 2.5227125941097004, 2.307028853640885, 2.816113682110262, 2.656776892632595, 2.217901707591825, 2.407606115982551, 3.0049571748896664, 2.2728388507970396, 1.9990976209898876, 2.6460654067747367, 2.4454311454187496, 1.8885460478001943, 2.125818820111916, 2.8086637099079255, 2.5194575118778832, 2.2708857033054373, 2.8559340758014655, 2.6738160820475017, 2.168517639481389, 2.3867296303830567, 3.0627033693378927, 2.1888501067167083, 1.915947942579093, 2.5584944374399097, 2.3612799696549374, 1.8062093945986786, 2.042094660880298, 2.6696618717688367, 2.274776845912477, 2.114899678595789, 2.5082595255524867, 2.3790256795115017, 2.0513175689975416, 2.1883179191422504, 2.6718927402241053, 2.3613407683453733, 2.089572971988282, 2.7321861830485967, 2.532007831907093, 1.9792198458492096, 2.2156202733674495, 2.9264262420610225]

**-->** number of indexes of the churchill's list:  49
**-->** number of churchill's indexes and turbulent : 49

**-->** turbulent list length:  49

**-->** turbulent minimum relative error:  1.8062093945986786  %
**-->** turbulent average relative error:  2.4098415275220235  %
**-->** turbulent maximum relative error:  3.0627033693378927  %

**-->** Relative roughness value:  0.043750000000000004 , Reynolds value:  9250.0 , minimum relative error:  1.8062093945986786 %
**-->** Relative roughness value:  0.018750000000000003 , Reynolds value:  4750.0 , maximum relative error:  3.0627033693378927 %

Process finished with exit code 0




