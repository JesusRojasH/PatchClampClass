# PatchClampClass
This repository has a class that can help with the patch clamp analisys technique using the pyabf and efel libraries avaliable in Python.

## Patch Clamp with Pyabf and Efel
### Introduction
This repository has the porpouse of bring a class that can help with the Patch Clamp technique used in neural synaptic comunication analysis, as well as in others kind of analysis that include electrical phenomena in neural activity. This class was created based on the Python libraries Pyabf and Efel (libraries focused in the previous mentioned analisys), the first part is dedicated to action potencials. 

***Warnings:***

- Some errors can arise during the install process, frecuently IPFX has the `metadata-generation-failed` error, this can be fixed using the command `pip install --no-deps ipfx` to install the library.

- This class has been tested in a venv with *Python 3.9.13* version, *efel 5.0.5*, *IPFX 1.0.8* and *pyabf 2.3.8*.

## Structure of PatchClamp()
The \_\_init\_\_ special function in python builds up the initial arguments in a funtion, those argument in the class are the follows:

__PatchClamp(DPath, T0, TF, S_Path):___

- _DPath_: Is the local path to the data that will be use in the analysis.

- _T0_ : Initial time of the data analisys.

- _TF_ : Final time of the data analisys.

- _S_Path_ : Local path where the analisys results will be saved.

### Class Modules

These are the modules of the class that help with potential analysis, most of them need two or three extra arguments to perform correctly the analisis, every extra argument will be explained in the correspondent module. 

__SFX(self[^1], \*\*kwargs[^2] )  (Spike Features from one or several traces):__ 

__Extra arguments \*\*kwargs:__

- _SwpNumber_: Is the number of the sweep highlighted on the plot.
- _Save_ : This is an option to save the analysis plot in a .pgn format and the results table in a .csv format. ('Y', default = __None__)

__STFX(self, **kwargs) (Spike Train Features from one or several trace):__ 

__Extra arguments \*\*kwargs:__

- _Ntraces_: This argument determine if the analysis will be perform for one  or several traces. ('One', 'Sev')
- _SwpNumber_: Is the number of the sweep highlighted on the plot.
- _Save_ : This is an option to save the analysis plot in a .pgn format and the results table in a .csv format. ('Y', default = __None__)

__STF_Spreads(self, CSVPath) (Spike Train Features form data in spreadsheets):__

- _CSVPath_: This is the local path to the csv file that containst the data for the analysis.

__SFX_Efel(self, SN) (Spike Features from one sweep):__

- _SN_ : Sweep number to highlight.

__STFX_Efel(self, SN) (Spike Train Features from several sweep):__  

- _SN_ : Sweep number to highlight.

__STFX_NA(self, TH, Save_T, F_Name) (Spike Train Features from data in Numpy Arrays):__

- _TH_ : Trace to highlight.
- _Save_T_ : Is an option to save the table analysis table in a csv document. ('Y', default = __None__)
- _F_Name_ : Name to save the analysis table.

__Phases(self, Data, Swp, SR, HW_ms, Save) (Phase-Plane Plot):__

- _Data_ : Data from which phase space will be obtained.
- _Swp_ : Select the trace in your recording.
- _SR_ : Is the rate of sample taken by unit of time (ms).
- _HW_ms_: Half window por milisecond. 
- _Save_: Option to save the phase plane in png format. 

  
[^1]: The self argument in a module allows the use of global variables of the class in that module. 
[^2]: \*\*kwargs is a special argument that allows to give a variable number of arguments in the module, you can modify the value of the variables if you use the variable key (name) correctly, \*\*kwargs default values is None.  
