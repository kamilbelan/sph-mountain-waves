"""
 Static atmosphere above a mountain with the Witch of Agnesi profile:

 h(x) = (h_m a²) / (x² + a²),

 all thermodynamic processes are adiabatic
"""

module StaticTEMPLATE

export run_sim

using DrWatson
using Dates
@quickactivate "SPH"

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

const FLUID = 0.0
const WALL = 1.0
const MOUNTAIN = 2.0


# ==============
# INCLUDE CORE SCRIPTS
# ==============

include(srcdir("core", "stationary_domain.jl"))
include(srcdir("core", "diagnostics.jl"))
include(srcdir("core", "time_loop.jl"))
include(srcdir("io", "data_storage.jl"))

# ==============
# INCLUDE UTILS SCRIPTS
# ==============

#include(srcdir("utils", "new_packing.jl"))

# ==============
# PARTICLE CONSTRUCTOR
# ==============

mutable struct Particle <: AbstractParticle
	h::Float64        # smoothing length
	Dh::Float64       # rate of smoothing length
	x::RealVector     # position
	m::Float64        # mass
	v::RealVector     # velocity
	Dv::RealVector    # acceleration
	ρ_bg::Float64     # background density
	ρ′::Float64       # density perturbation
	ρ::Float64        # total density
	Dρ::Float64       # density rate
	n::Float64        # number density
	Omega::Float64    # ∇h correction factor
	P_bg::Float64     # background pressure
	P′::Float64       # pressure perturbation
	P::Float64        # total pressure
	c::Float64        # local speed of sound
	θ_bg::Float64     # bakcground potential temperature
	θ′::Float64       # potential temperature perturbation
	θ::Float64        # total potential temperature
	T_bg::Float64     # background temperature
	T′::Float64       # temperature perturbation
	T::Float64        # total temperature
	type::Float64     # particle type

	function Particle(x::RealVector, v::RealVector, type::Float64, global_params::Dict, sim_params::Dict)
		# unpack all parameters
		@unpack g, R_mass, cp, cv, γ, R_gas, T_bg, ρ0, N = global_params
		@unpack dom_height, dom_length, a, z_t = global_params
		@unpack h_m, z_β = sim_params
		@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
		@unpack η, dr, dt_rel, t_end = sim_params

		# Derived values within the constructor
		h0 = η * dr
		obj = new(
			h0,             # h 
			0.0,            # Dh
			x,              # x 	 	
			0.0,            # m
			v,              # v
			VEC0,           # Dv
			0.0,            # ρ_bg
			0.0,            # ρ′
			0.0,            # ρ
			0.0,            # Dρ
			0.0,            # n
			0.0,            # Omega
			0.0,            # P_bg
			0.0,            # P′
			0.0,            # P
			0.0,            # c
			0.0,            # θ_bg 
			0.0,            # θ′ 
			0.0,            # θ 
			0.0,            # T_bg
			0.0,            # T′
			0.0,            # T
			type,           # type
		)

		# initialization
		obj.T_bg = T_bg
		obj.ρ_bg = 0.0 # needs SPH sum!
		obj.P_bg = background_pressure(obj.x[2],ρ0, T_bg, g, R_mass)
		obj.θ_bg = background_pot_temperature(obj.x[2],ρ0, T_bg, g, R_mass, R_gas)

		obj.ρ′ = 0.0
		obj.P′ = 0.0
		obj.T′ = 0.0
		obj.θ′ = 0.0

		obj.T = obj.T′ + T_bg
		obj.ρ = obj.ρ′ + obj.ρ_bg
		obj.P = obj.P′ + obj.P_bg
		obj.θ = obj.θ′ + obj.θ_bg

		obj.m = ρ0 * dr^2
		obj.c = sqrt(γ * obj.P / obj.ρ)
		return obj
	end
end

# ==============
# Background values (e.g density, pressure, potential temperature, entropy)
# ==============

function background_density(y::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)
	return ρ0 * exp(-y * g / (R_mass * T_bg))
end

function background_pressure(y::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)
	ρ_bg = background_density(y, ρ0, T_bg, g, R_mass)
	return R_mass * T_bg * ρ_bg
end

function background_pot_temperature(y::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64, R_gas::Float64)
	P_bg = background_pressure(y, ρ0, T_bg, g, R_mass)
	return T_bg * (((T_bg * R_gas * ρ0) / P_bg))^(2 / 7)
end

# ==============
# Pressure computation (e.g followed by sound speed computation)
# ==============

# ==============
# Thermodynamics (e.g determining temperature and potential temperature)
# ==============

# ==============
# Smoothing-length evolution (e.g. computing the SPH sum, setting the adaptive h, reseting the rate...)
# ==============

@inbounds function update_smoothing!(p::Particle, η::Float64, h0::Float64)
	omega = max(p.Omega, 0.01)
	# Newton iteration to find the suitable h_new
	h_new = p.h - (p.h - η / sqrt(p.n)) / omega 
	p.h = clamp(h_new, 0.8 * p.h, 1.2 * p.h)
	p.h = min(p.h, 3 * h0) # precaution for massive expansion
end

# ==============
# Additional forces: Rayleigh damping & gravity/buyoancy
# ==============

function damping_structure(z::Float64, v::RealVector, z_t::Float64, z_β::Float64, γ_r::Float64)
	z_bottom = z_t - z_β
	if z >= z_bottom
		ζ = (z - z_bottom) / z_β 

		profile = (sin(π / 2 * ζ))^2 

		# friction force
		return -γ_r * profile * v 
	else
		return VEC0
	end
end

function buyoancy_force(p::Particle, g::Float64)
	return -g * VECY 

end


# ==============
# Momentum balance
# ==============



# ==============
# Mass balance 
# ==============



# ==============
# Move & accelerate
# ==============

function move!(p::Particle, dt::Float64)
	if p.type == FLUID
		p.x += dt * p.v
	end
end

function accelerate!(p::Particle, dt::Float64, g::Float64, z_t::Float64, z_β::Float64, γ_r::Float64)
	if p.type == FLUID
		p.v += 0.5 * dt * (p.Dv + buyoancy_force(p, g) + damping_structure(p.x[2], p.v, z_t, z_β, γ_r)) # this is a vector sum
	end
end

function reset_acceleration!(p::Particle)
	p.Dv = VEC0
end

# ==============
# Verlet step
# ==============

function verlet_step!(sys, global_params, sim_params)
	# unpack all parameters
	@unpack g, R_mass, cp, cv, γ, R_gas, T_bg, ρ0, N = global_params
	@unpack dom_height, dom_length, a, z_t = global_params
	@unpack h_m, z_β = sim_params
	@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel = sim_params

	# compute derived parameters
	h0 = η * dr
	c = sqrt(65e3 * (γ) / ρ0)
	dt = dt_rel * h0 / c
	γ_r = γ_r_rel * N
	
end

# ==============
# Main Entry Point
# ==============

function run_sim(global_params::Dict, sim_params::Dict)
	# ==============
	# Parameters initialization
	# ==============
	
	# unpack all parameters
	@unpack g, R_mass, cp, cv, γ, R_gas, T_bg, ρ0, N = global_params
	@unpack dom_height, dom_length, a, z_t = global_params
	@unpack h_m, z_β = sim_params
	@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel = sim_params

	# compute derived parameters
	h0 = η * dr
	bc_width = 6*dr
	c = sqrt(65e3 * (γ) / ρ0)
	dt = dt_rel * h0 / c
	dt_frame = t_end / 100
	γ_r = γ_r_rel * N

	# create the particle system
	sys = make_system(Particle, global_params, sim_params)

	# ==============
	# Initialization of the physical fields
	# ==============
	
	# initialization of the pressure

	# compute temperature and potential temperature

	# compute acceleration (balance of momentum)

	# ==============
	# Output handling
	# ==============
	
	# initialize simulation output directory with metadata
	run_dir = initialize_run_directory(sim_params)

	# initialize data accumulation for diagnostics
	average_velocities, maximum_velocities, energies = initialize_diagnostics_arrays()

	# ==============
	# Time loop
	# ==============
	
	# use the chosen verlet_step! to advance in time
	step_function!(sys) = verlet_step!(sys, global_params, sim_params)

	# iterate in time
	time_loop(global_params, run_dir, step_function!, t_end, dt, dt_frame, sys, average_velocities, maximum_velocities, energies)

	# create a xdmf file to use with the generated h5 files
	generate_sph_xdmf(run_dir)

	# write CSVs once at the end, create diagnostic plots
	finalize_diagnostics(run_dir, average_velocities, maximum_velocities, energies)

	println("\nEnd of the road 🏵.")
	println("Simulation $(basename(run_dir)) complete.")
	println("All artifacts (data, plots, metadata) securely saved to:\n  $run_dir")

	return run_dir
end
end # module

