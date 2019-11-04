
## MATLAB tools for LII analysis (wat-lii)

This is modular program to model and analyze time-resolved laser-induced
incandescence (TiRe-LII) signals.


## Details

This program is built to simulate signals from various materials,
including soot, silicon, germanium, iron, silver, and molybdenum.
Signals are generated predominantly using absorption, conduction,
and evaporation submodels, with capabilities to do other cooling
modes.

## Classes

#### @HTModel

This class is designed to generate temperature
decay curves. This is done by solving, at the very least, an ordinary
differential equation for temperature. Mass and annealed fraction can also
be solved simultaneously.

Various aspects of the `HTModel` class are controlled using an `opts` structure,
which is a property of the object.

#### @SModel

This class is design to simulate incandescence
from a temperature trace, incorporating blackbody
radiation and the optical properties. Methods exist
to both forward and backward (i.e. pyrometry)
calculations.

#### @Prop

This class contains the thermosphysical, optical, and other model parameters
to be used in evaluating both the spectroscopic and heat transfer models. These
parameters are also the ones that can be selected as dependent variables
when creating instances of the `HTModel` and `SModel` classes.

Note that this class is subclass of the `dynamicprops` class, such
that new properties can be added using the inherited `addprops` method
described in the MATLAB literature.

Note that this class is also currently a subclass of the `handle` class, such
that instances act like as "pass-by-reference". This means that changes to the
instance within other functions act on the original instance. Accordingly,
a `copy` method is included with the `Prop` class to allow for independent
copies to be made. Simple assignments (e.g. `prop2 = prop;`) will be
insufficient to create an independent copy, as changes to `prop2` will
affect the original `prop`.

----------------------------------------------------------------------

#### License

This software is licensed under an MIT license (see the corresponding file
for details).


#### Contact information

The primary author of the code is Timothy A. Sipkens, who can be
emailed at [tsipkens@uwaterloo.ca](mailto:tsipkens@uwaterloo.ca).
