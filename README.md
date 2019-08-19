# LIME-juicer
Python wrapper on C LIME radiative transfer code. 

LIME (https://github.com/lime-rt/lime) is a 3D radiative transfer code originally written by Christian Brinch, updated and maintained by Ian Stewart. 

Provided here is a python wrapper on the C code which can be used to produce LIME executables for the desired input model parameters. The wrapper can be run as a standalone for quick RT simulations with a single set of parameters, or can be imported into external python scripts for more complicated RT calculations e.g. over a grid of parameters.

Currently, the wrapper is written specifically to compute radiative transfer through protoplanetary disk models, but in the future, will be extended to work more generally and for a broader set of model geometries.

Stay tuned for more detailed installation and use instructions [Under Construction].
