## wat-lii-analysis

This is modular program to model and analyze time-resolved laser-induced
incandescence (TiRe-LII) signals.


### Details

This program is built to simulate signals from various materials,
including soot, silicon, germanium, iron, silver, and molybdenum.
Signals are generated predominantly using absorption, conduction,
and evaporation submodels, with capabilities to do other cooling
modes.

### Classes

##### @HTModel

This class is designed to generate temperature 
decay curves. 

##### @SModel

This class is design to simulate incandescence 
from a temperature trace, incorporating blackbody
radiation and the optical properties. Methods exist
to both forward and backward (i.e. pyrometry) 
calculations. 

----------------------------------------------------------------------

#### License

This software is licensed under an MIT license (see the corresponding file
for details).


#### Contact information

The primary author of the code is Timothy A. Sipkens, who can be
emailed at [tsipkens@uwaterloo.ca](mailto:tsipkens@uwaterloo.ca).
