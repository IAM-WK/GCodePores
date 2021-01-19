# GCodePores

A tool to process GCode from various 3D printers and correlate information about pores obtained from µCT images.
Further development is continued under https://sourceforge.net/p/ctfam/wiki/GCodePoreCorrelation/.
Various related tools will be availabe under https://sourceforge.net/p/ctfam/

- - - -
## Supported GCode Flavours
* Cura
* Orlas Suite
* Arburg Freeformer

- - - -

## Getting Started
Prerequisites:
* .mhd file : for physical image dimensions
* list of analysed pores from 
* GCode file
* [Paraview](https://www.paraview.org/)

<p align="center">
  <img src="https://github.com/IAM-WK/GCodePores/blob/master/doc/GCodePores_ProgramFlow.png" width="350" title="Flowchart">
</p>

- - - -

### Preprocessing µCT Images
Rationale:
1. Segment image in pores and material
2. label pores with [MorphoLibJ](https://imagej.net/MorphoLibJ)
3. analyse pores with [MorphoLibJ](https://imagej.net/MorphoLibJ)
4. save pore list to a file

### Analysing Data with GCodePores
---
Step 1: registering the GCode to your µCT Scan:

Execute GCodePores with parameters like: 
GCodePores.exe -i image.mhd -g sampleobject.gcode --interpolationoff --vtkfileonly -a 0 -u "0 0 0" -k 1 -p MorphoLibJoutput.csv \
(-k and -p are needed for the program to run but in this step irrelevant)

Open the image.mhd in Paraview. The program will generate a coordinatefile.vtk which can be opened alongside in Paraview. Now register the coordinatefile with Paraview (Transforming->Translation and Transforming->Rotation) until the path matches onto the image.\
Please only use rotation around z-axis! If other rotations are needed for registration, register the image until at least rotations in x- and y-axis aren't needed, save the image in this orientation and do the MorphoLibJ analysis on this registered Image!

Execute GCodePores again with the obtained Translation under -u and z-axis rotation under -a, for example: \
GCodePores.exe -i image.mhd -g sampleobject.gcode --interpolationoff --vtkfileonly -a 1.2423 -u "3.2 4.2 1.2" -k 1 -p MorphoLibJoutput.csv \
Check again in Paraview if the registration result with the new coordinatefile.vtk is good.

---

Step 2: starting an analysis \
Execute GCodePores with parameters like: \
GCodePores.exe -i image.mhd -g sampleobject.gcode -d 0.1 -a 1.2423 -u "3.2 4.2 1.2" -k 3 -p MorphoLibJoutput.csv -j 7  -o AnalysisResults.csv\
Check GCodePores.exe --help for additional options and descriptions.
The Analysis results can be found in AnalysisResults.csv

### Troubleshooting

Help is available through GCodePores.exe -h or --help


### Dependencies 

Boost: program_options


## Authors

tbd
