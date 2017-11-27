[![Build Status](https://travis-ci.org/JuliaInv/FWI.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/FWI.jl) [![Coverage Status](https://coveralls.io/repos/github/JuliaInv/FWI.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/FWI.jl?branch=master)

# FWI.jl
A Julia package for solving the Full Waveform Inversion on a regular rectangular mesh.
For forward modelling and sensitivities it uses either a direct solver or a block Shifted Laplacian 
Multigrid preconditioner with BiCGSTAB.

More detailed information about this package can be found in the papers:

(1) Lars Ruthotto, Eran Treister and Eldad Haber, jInv--a flexible Julia package for PDE parameter estimation, SIAM Journal on Scientific Computing, 39 (5), S702-S722, 2017.

(2) Eran Treister and Eldad Haber, Full waveform inversion guided by travel time tomography, SIAM Journal on Scientific Computing, 39 (5), S587-S609, 2017. 

Specifically, paper (2) describes most of the algorithms in the package in details. For the travel time tomography, see EikonalInv.jl.  The package also has an option for data generation in time domain (in parallel) and conversion to frequency domain.



# Requirements

This package is intended to use with julia versions 0.5.x.

This package is an add-on for jInv, which needs to be installed.

This package uses: EikonalInv.jl (for seismic utils), Multigrid.jl, ForwardHelmholtz.jl, MAT.jl. 
# Installation

```
Pkg.clone("https://github.com/simonster/MAT.jl.git","MAT");
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.clone("https://github.com/JuliaInv/FactoredEikonalFastMarching.jl","FactoredEikonalFastMarching")
Pkg.clone("https://github.com/JuliaInv/EikonalInv.jl","EikonalInv")
Pkg.clone("https://github.com/JuliaInv/Multigrid.jl","Multigrid")
Pkg.clone("https://github.com/JuliaInv/ForwardHelmholtz.jl","ForwardHelmholtz")

Pkg.test("FWI")
```

# Examples
Under "examples/" you can find the 2D experiments presented in the paper (2) above.

