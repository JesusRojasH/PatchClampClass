# PatchClampClass
This repository has a class that can help with the patch clamp analisys technique using the pyabf and efel libraries avaliable in Python.

## Patch Clamp with Pyabf and Efel
### Introduction
This repository has the porpouse of bring a class that can help with the Patch Clamp technique used in neural synaptic comunication analysis, as well as in others kind of analysis that include electrical phenomena in neural activity. This class was created based on the Python libraries Pyabf and Efel (libraries focused in the previous mentioned analisys), the first part is dedicated to action potencials.

## Structure of PatchClamp()
The \_\_init\_\_ special function in python builds up the initial arguments in a funtion, those argument in the class are the follows:

__PatchClamp(DPath, T0, TF, S_Path):___

- _DPath_: Is the local path to the data that will be use in the analysis.

- _T0_ : Initial time of the data analisys.

- _TF_ : Final time of the data analisys.

- _S_Path_ : Local path where the analisys results will be saved.

### Class Modules

These are the modules of the class that help with potential analysis, most of them need two or three extra arguments to permorm correctly the analisis, every extra argument will be explained in the correspondent module. 

__SFX(self[^1], \*\*kwargs[^2] )  (Spike Features):__ 

__Extra arguments \*\*kwargs:__

- _SwpNumber_: Is the number of the sweep highlighted on the plot.
- _Save_ : This is an option to save the analysis plot in a .pgn format and the results table in a .csv format. 


[^1]: The self argument in a module allows the use of global variables of the class in that module. 
[^2]: \*\*kwargs is a special argument that allows to give a variable number of arguments in the module, you can modify the value of the variables if you use the variable key (name) correctly, \*\*kwargs default values is None.  
