[![Build Status](https://travis-ci.org/JuliaInv/FWI.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/FWI.jl) [![Coverage Status](https://coveralls.io/repos/github/JuliaInv/FWI.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/FWI.jl?branch=master)

# FWI.jl
A Julia package for solving the Full Waveform Inversion on a regular rectangular mesh.
For forward modelling and sensitivities it uses either a direct solver or a block Shifted Laplacian 
Multigrid preconditioner with BiCGSTAB.

# Requirements

This package is intended to use with julia versions 0.4.x.

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

Under "examples/" you can find the 2D experiment.