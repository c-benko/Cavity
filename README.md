# Cavity
Simulations of various important properties of femtosecond enhancement cavities.

# File descriptions
`mirrors.py` is a class that contains various mirror properties like dispersion and reflectivities.

`cavity_geo.py` runs a simulation using the ABCD matrices to determine geometric properties of the beam inside a cavity.

`cavity_pulse.py` runs a simulation of the circulating intracavity pulse as a function of cavity mirrors, input pulse duration and spectral bandwidth and outputs the circulating intracavity spectrum and pulse. 

`\Data` contains dispersion data of the mirrors used in the experiments. 

`plotter.py` is just a display function.

`hole.py` runs a simulation of the optical loss induced by a circular aperture in one of the mirrors. 
