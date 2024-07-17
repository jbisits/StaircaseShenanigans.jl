# StaircaseShenanigans.jl

<!--
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jbisits.github.io/StaircaseShenanigans.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jbisits.github.io/StaircaseShenanigans.jl/dev/)
-->
[![Build Status](https://github.com/jbisits/StaircaseShenanigans.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jbisits/StaircaseShenanigans.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package for setting up Direct Numerical Simulations (DNS) of thermohaline staircases using [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl).
This is currently work in progress and breaking changes are likley.

To add the package (assuming Julia is already installed):

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/jbisits/StaircaseShenanigans.jl")
```

To then use the package

```julia
julia> using StaircaseShenanigans
```

The package provides functions and methods to setup DNS experiments with salt and temperature staircases as initial conditions.