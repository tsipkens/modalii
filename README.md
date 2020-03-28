
## MATLAB tools for LII analysis (wat-lii)

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

This is modular program to model and analyze time-resolved laser-induced incandescence (TiRe-LII) signals, developed at the University of Waterloo.

This program is built to simulate signals from various materials, including soot, silicon, germanium, iron, silver, and molybdenum. Signals are generated predominantly using absorption, conduction, and evaporation submodels, with capabilities to do other cooling
modes. The program contain the following components.

## 1. Upper directory and other main\*.m scripts

Throughout the program, `main*.m` scripts are used to create instances of the classes and perform the analysis of LII signals.

## 2. Classes

#### 2.1 @HTModel

This class is designed to generate temperature decay curves. This is done by solving, at the very least, an ordinary differential equation for temperature. Mass and annealed fraction can also be solved simultaneously.

One can create an instance of the HTModel class by calling the 
construction methods as follows:

```Matlab
htmodel = HTModel(prop,x_fields,t,opts);
```

Here, `prop` is an instant of the `Prop` class, which contains all of
the physical parameters required to define the heat transfer model.
The parameter `x_fields` then contains a cell of strings, where each entry
is a property of the given `prop` variable. The `t` input if a 
vector of time for which the heat transfer model will be evaluated. 
Finally, various aspects of the `HTModel` class are controlled using an `opts` structure,
which is a property of the object.

The two key methods for evaluating the heat transfer model
are the `evaluate` and `de_solve` methods. The `de_solve` method
solves the ODEs without altering the default physical properties. 
In contrast, the `evaluate` method solves the ODEs for a vector of
property values given by the `x` input to the method. The latter
is particularly useful in optimization scenerio. 

#### 2.2 @SModel

This class is design to simulate incandescence
from a temperature trace, incorporating blackbody
radiation and the optical properties. Methods exist
to both forward and backward (i.e. pyrometry)
calculations.

#### 2.3 @Prop

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

#### @Signal

This is depreciated class that is to be replaced by structured arrays. It packages together a series of signals and related information, such as the time and wavelengths. The class is still included for legacy purposes. The class's previous methods were largely moved to the `+data` package, described below. 

## 3. Packages

#### 3.1 +data

The data package is available to filter or otherwise process data. 

----------------------------------------------------------------------

#### License

This software is licensed under an MIT license (see the corresponding file for details).


#### Contact information and acknowledgements

The primary author of the code is Timothy A. Sipkens, who can be emailed at  [tsipkens@uwaterloo.ca](mailto:tsipkens@uwaterloo.ca). The code was developed at the University of Waterloo.  Kyle Daun contributed significantly to the ideas summarized in this code. 
