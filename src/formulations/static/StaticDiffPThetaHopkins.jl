"""
Static atmosphere above a mountain with the Witch of Agnesi profile:

h(x) = (h_m a²) / (x² + a²),

all thermodynamic processes are adiabatic
"""

module StaticDiffPThetaHopkins

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

include(srcdir("core", "domain.jl"))
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
	x::RealVector     # position
	m::Float64        # mass
	v::RealVector     # velocity
	Dv::RealVector    # acceleration
	ρ::Float64        # total density
	n::Float64        # number density
	Omega::Float64    # ∇h correction factor
	P_bg::Float64     # background pressure
	P::Float64        # total pressure
	θ_bg::Float64     # bakcground potential temperature
	θ′::Float64       # potential temperature perturbation
	θ::Float64        # total potential temperature
	T_bg::Float64     # background temperature
	T::Float64        # total temperature
	type::Float64     # particle type

	function Particle(x::RealVector, v::RealVector, type::Float64, global_params::Dict, sim_params::Dict)
		# unpack all parameters
		@unpack g, R_mass, cp, cv, γ, R_gas, T_bg, ρ0, N = global_params
		@unpack dom_height, dom_length, h_m, a, z_t, z_β = global_params
		@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
		@unpack η, dr, dt_rel, t_end = sim_params

		# Derived values within the constructor
		h0 = η * dr
		obj = new(
			h0,             # h 
			x,              # x 	 	
			0.0,            # m
			v,              # v
			VEC0,           # Dv
			0.0,            # ρ
			0.0,            # n
			0.0,            # Omega
			0.0,            # P_bg
			0.0,            # P
			0.0,            # θ_bg 
			0.0,            # θ′ 
			0.0,            # θ 
			0.0,            # T_bg
			0.0,            # T
			type,           # type
		)

		# initialization

		obj.ρ = 0.0 # needs SPH sum!
		obj.m = ρ0 * dr^2

		obj.T_bg = T_bg
		obj.T = T_bg

		obj.P_bg = 0.0 # needs SPH sum!
		obj.P = 0.0 # needs SPH sum!

		obj.θ_bg = background_pot_temperature(obj.x[2],ρ0, T_bg, g, R_mass, R_gas)
		obj.θ′ = 0.0
		obj.θ = background_pot_temperature(obj.x[2],ρ0, T_bg, g, R_mass, R_gas)

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

@inbounds function reset_pressure!(p::Particle, θ_0::Float64)
	ker_self = wendland2(p.h, 0.0)
	p.P = (p.m / θ_0) * p.θ * ker_self
	p.P_bg = (p.m / θ_0) * p.θ_bg * ker_self
end

@inbounds function compute_pressure!(p::Particle, q::Particle, r::Float64, θ_0::Float64)
	ker = wendland2(p.h, r)
	p.P += (q.m / θ_0) * p.θ * ker
	p.P_bg += (q.m / θ_0) * p.θ_bg * ker
end

@inbounds function finalize_pressure!(p::Particle, γ::Float64)
	p.P = p.P^γ
	p.P_bg = p.P_bg^γ
end

# ==============
# Thermodynamics 
# ==============

@inbounds function compute_pot_temperature!(p::Particle, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64, R_gas::Float64)
	p.θ_bg = background_pot_temperature(p.x[2], ρ0, T_bg, g, R_mass, R_gas)
	p.θ′ = p.θ - p.θ_bg
end

@inbounds function find_temperature!(p::Particle, R_mass::Float64)
	p.T = p.P / (R_mass * p.ρ)
end

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

# NO BUYOANCY EXPLICITELY!

# ==============
# Momentum balance
# ==============

@inbounds function balance_of_momentum!(p::Particle, q::Particle, r::Float64, α::Float64, β::Float64, ε::Float64, rho_floor::Float64, P_floor::Float64, θ_0::Float64, γ::Float64)
	x_pq = p.x - q.x
	v_pq = p.v - q.v
	dot_product = SmoothedParticles.dot(x_pq, v_pq)

	prefac = q.m * p.θ * q.θ / θ_0^2
	expfac = 1.0 - 2.0 / γ
	ker_i = rDwendland2(p.h, r)
	ker_j = rDwendland2(q.h, r)
	f_i = 1.0 / p.Omega
	f_j = 1.0 / q.Omega
	pP = max(P_floor, p.P)
	qP = max(P_floor, q.P)

	# acceleration due the gradient of total pressure
	a_tot = -prefac * (f_i * pP^expfac * ker_i + f_j * qP^expfac * ker_j) * x_pq

	prefac_bg = q.m * p.θ_bg * q.θ_bg  / θ_0^2
	pP_bg = max(P_floor, p.P_bg)
	qP_bg = max(P_floor, q.P_bg)

	# acceleration due to the gradient of background pressure
	a_bg = -prefac_bg * (f_i * pP_bg^expfac * ker_i + f_j * qP_bg^expfac * ker_j) * x_pq

	# total acceleration
	p.Dv += a_tot - a_bg

	# artificial viscous force
	if dot_product < 0.0
		h_ij = 0.5 * (p.h + q.h)
		ker_ij = rDwendland2(h_ij, r)
		prho = max(p.ρ, rho_floor)
		qrho = max(q.ρ, rho_floor)
		c_i = sqrt(γ * p.P / prho)
		c_j = sqrt(γ * q.P / qrho)
		c_ij = 0.5 * (c_i + c_j)
		ρ_ij = 0.5 * (prho + qrho)
		μ_ij = (h_ij * dot_product) / (r * r + ε * h_ij * h_ij)
		π_ij = (-α * c_ij * μ_ij + β * μ_ij * μ_ij) / ρ_ij

		# artificial viscous force
		p.Dv += -q.m * π_ij * ker_ij * x_pq
	end
end


# ==============
# Mass balance 
# ==============

@inbounds function reset_density!(p::Particle)
	ker_self = wendland2(p.h, 0.0)
	p.ρ = p.m * ker_self
end

@inbounds function compute_density!(p::Particle, q::Particle, r::Float64)
	p.ρ += q.m * wendland2(p.h, r)
end

# ==============
# Number density
# ==============

@inbounds function reset_number_density!(p::Particle)
	p.n = wendland2(p.h, 0.0)
	p.Omega = 0.0
end

@inbounds function compute_number_density!(p::Particle, q::Particle, r::Float64)
	p.n += wendland2(p.h, r)
	p.Omega += -0.5 * r^2 * rDwendland2(p.h, r) # derived through chain rule
end

@inbounds function finalize_number_density!(p::Particle)
	p.Omega = 1.0 + (p.Omega / p.n)
end


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
	@unpack dom_height, dom_length, h_m, a, z_t, z_β = global_params
	@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel = sim_params

	# compute derived parameters
	h0 = η * dr
	c = sqrt(65e3 * (γ) / ρ0)
	dt = dt_rel * h0 / c
	γ_r = γ_r_rel * N
	P_0 = ρ0 * R_mass * T_bg
	θ_0 = (P_0^((γ - 1) / γ)) / R_mass

	# half-step acceleration & drift
	apply!(sys, p -> accelerate!(p, dt, g, z_t, z_β, γ_r))
	apply!(sys, p -> move!(p, dt))
	create_cell_list!(sys)

	# set the number density and adjust the smoothing length accordingly
	max_iter = 3
	for iter = 1:max_iter
		apply!(sys, p -> reset_number_density!(p))
		apply!(sys, (p, q, r) -> compute_number_density!(p, q, r))
		apply!(sys, p -> finalize_number_density!(p))

		# we in fact do not need the last values, since they will be taken into the balance of momentum
		if iter < max_iter
			apply!(sys, p -> update_smoothing!(p, η, h0))
			create_cell_list!(sys)
		end
	end

	# compute density and smoothing length
	apply!(sys, p -> reset_density!(p))
	apply!(sys, (p, q, r) -> compute_density!(p, q, r))

	# compute pressure
	apply!(sys, p -> reset_pressure!(p, θ_0))
	apply!(sys, p -> compute_pot_temperature!(p, ρ0, T_bg, g, R_mass, R_gas))
	apply!(sys, (p, q, r) -> compute_pressure!(p, q, r, θ_0))
	apply!(sys, p -> finalize_pressure!(p, γ))

	# compute temperature 
	apply!(sys, p -> find_temperature!(p, R_mass))
	# compute the forces
	apply!(sys, p -> reset_acceleration!(p))
	apply!(sys, (p, q, r) -> balance_of_momentum!(p, q, r, α, β, ϵ, rho_floor, P_floor, θ_0, γ))
	apply!(sys, p -> accelerate!(p, dt, g, z_t, z_β, γ_r))
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
	@unpack dom_height, dom_length, h_m, a, z_t, z_β = global_params
	@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel = sim_params

	# compute derived parameters
	h0 = η * dr
	c = sqrt(65e3 * (γ) / ρ0)
	dt = dt_rel * h0 / c
	dt_frame = t_end / 100
	γ_r = γ_r_rel * N
	P_0 = ρ0 * R_mass * T_bg
	θ_0 = (P_0^((γ - 1) / γ)) / R_mass

	# create the particle system
	sys = make_system(Particle, global_params, sim_params)

	# ==============
	# Initialization of the physical fields
	# ==============

	# initialize the number density and the smoothing length
	max_iter = 3
	for iter in 1:max_iter
		apply!(sys, p -> reset_number_density!(p))
		apply!(sys, (p, q, r) -> compute_number_density!(p, q, r))
		apply!(sys, p -> finalize_number_density!(p))

		if iter < max_iter
			apply!(sys, p -> update_smoothing!(p, η, h0))
			create_cell_list!(sys)
		end
	end

	# initialize the density
	apply!(sys, p -> reset_density!(p))
	apply!(sys, (p, q, r) -> compute_density!(p, q, r))

	# initialization of the pressure
	apply!(sys, p -> reset_pressure!(p, θ_0))
	apply!(sys, p -> compute_pot_temperature!(p, ρ0, T_bg, g, R_mass, R_gas))
	apply!(sys, (p, q, r) -> compute_pressure!(p, q, r, θ_0))
	apply!(sys, p -> finalize_pressure!(p, γ))

	# compute temperature 
	apply!(sys, p -> find_temperature!(p, R_mass))

	# compute acceleration (balance of momentum)
	apply!(sys, p -> reset_acceleration!(p))
	apply!(sys, (p, q, r) -> balance_of_momentum!(p, q, r, α, β, ϵ, rho_floor, P_floor, θ_0, γ))
	apply!(sys, p -> accelerate!(p, dt, g, z_t, z_β, γ_r))

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

