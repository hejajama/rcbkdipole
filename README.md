# rcbkdipole
This code provides rcBK evolved dipole amplitude

Based on 
T. Lappi, H. Mäntysaari
Phys.Rev. D88 (2013) 114020, arXiv:1309.6963 

The "mvgamma" fit parameters are from Albacete, Armesto, Milhano, Quiroga-Arias, Salgado, Eur. Phys. J C71, 1705, arXiv:1012.4408

Initial condition for the LO BK evolution is fitted to the HERA data.
The resulting dipole amplitudes can be found in data tables in data directory

This code interpolates in dipole size and rapidity/Bjorken-x

## Building 
Requires
* Cmake
* GSL

How to compile:
```
 mkdir build
 cd build
 cmake ..
 make
```

This generates a library build/lib/libamplitude.a that you can link in your own program.

## How to use 
See the example code in file ``src/dipole_amplitude.cpp``

## Data tables for protons
* data/proton/mv.dat MV parametrization
* data/proton/mve.dat MVe parametrization
* data/proton/mvgamma.dat MVgamma parametrization

## Data tables for Pb
* data/Pb/mve/glauber_mve_X 

Here X refers to the impact parameter in GeV^(-1), all impact parameters are evolved independently.

## Questions and comments
This should be robust and maintained code, if you find bugs, please send an email to heikki.mantysaari@jyu.fi
