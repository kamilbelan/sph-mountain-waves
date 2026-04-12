using DrWatson
@quickactivate "SPH"

# load all packages used by formulations
using Printf
using SmoothedParticles
using DataFrames
using Plots
using JLD2
using HDF5
using Glob
using CSV
using Parameters
using LinearAlgebra
using TOML

# Import formulation modules to trigger JIT compilation
formulations = [
	"EvolWBalancedWCSPH",
	"EvolWBalancedHopkins",
	"EvolWBalancedThetaHopkins"
]

for formulation in formulations
	model_path = srcdir("formulations", "evolutionary", "$(formulation).jl")
	if isfile(model_path)
		include(model_path)
		println("Precompiled: $formulation")
	else
		println("Warning: $model_path not found")
	end
end

println("-------sysimage precompilation complete-------")
