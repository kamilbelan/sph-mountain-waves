"""
 An isentropic flow around a mountain with the Witch of Agnesi profile:

 h(x) = (hₘ a²) / (x² + a²),
"""

module EvolWBalancedWCSPH

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


# ==============
# INCLUDE CORE SCRIPTS
# ==============

include(srcdir("core", "evolutionary_domain.jl"))
include(srcdir("core", "diagnostics.jl"))
include(srcdir("core", "time_loop.jl"))
include(srcdir("io", "data_storage.jl"))
include(srcdir("core", "bc", "ebc_shared.jl"))

# ==============
# INCLUDE FORMULATION-SPECIFIC SCRIPTS
# ==============

include(srcdir("core", "bc", "ebc_wcsph.jl"))

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
	ρ::Float64        # total density
	ρ_bg::Float64     # background density
	δρ::Float64       # density perturbation
	Dρ::Float64       # density rate
	P::Float64        # total pressure
	P_bg::Float64     # background pressure
	δP::Float64       # pressure perturbation
	c::Float64        # local speed of sound
	θ::Float64        # total potential temperature
	θ_bg::Float64     # bakcground potential temperature
	δθ::Float64       # potential temperature perturbation
	T_bg::Float64     # background temperature
	T::Float64        # total temperature
	type::Float64     # particle type
	id::Int64              # the index of the particle in sys.particles
	spawn_y::Float64       # the altitude it was generated at
	grad_δρ::RealVector     # ∇δρ
	grad_u::RealVector     # ∇u (x-velocity gradient)
	grad_w::RealVector     # ∇w (y-velocity gradient)
	best_match_id::Int64   # rd of the fluid particle 'e'
	min_dist_sq::Float64   # argmin distance tracker

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
			0.0,            # ρ
			0.0,            # ρ_bg
			0.0,            # δρ
			0.0,            # Dρ
			0.0,            # P
			0.0,            # P_bg
			0.0,            # δP
			0.0,            # c
			0.0,            # θ 
			0.0,            # θ_bg 
			0.0,            # δθ
			0.0,            # T_bg
			0.0,            # T
			type,            # type
	                0,              # the index of the particle in sys.particles
			0.0,            #spawn_y::Float64
			VEC0,           #grad_ρ::RealVector
			VEC0,           #grad_u::RealVector
			VEC0,           #grad_w::RealVector
			0,              #best_match_id::Int64
			0.0,            #min_dist_sq::Float64
		)

		# initialization
		
		obj.ρ_bg = background_density(obj.x[2], ρ0, T_bg, g, R_mass)
		obj.δρ = 0.0
		obj.ρ = obj.δρ + obj.ρ_bg
		obj.m = ρ0 * dr^2

		obj.T_bg = T_bg
		obj.T = T_bg

		obj.P_bg = background_pressure(obj.x[2],ρ0, T_bg, g, R_mass)
		obj.δP = 0.0
		obj.P = obj.δP + obj.P_bg
		obj.c = sqrt(γ * obj.P / obj.ρ)

		obj.θ_bg = background_pot_temperature(obj.x[2],ρ0, T_bg, g, R_mass, R_gas)
		obj.δθ = 0.0
		obj.θ = obj.δθ + obj.θ_bg

		obj.spawn_y = obj.x[2]

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

function set_bulk_velocity!(p::Particle, v_initial::Float64)
	if p.type == FLUID
		p.v = v_initial * VECX
	end
end

# ==============
# Pressure computation (e.g followed by sound speed computation)
# ==============

@inbounds function compute_pressure!(p::Particle, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64, P_floor)
	p.P_bg = background_pressure(p.x[2], ρ0, T_bg, g, R_mass)
	p.δP = p.c^2 * p.δρ
	pP = p.P_bg + p.δP
	p.P = max(pP, P_floor)
end


# ==============
# Thermodynamics (e.g determining temperature and potential temperature)
# ==============

@inbounds function find_temperature!(p::Particle, R_mass::Float64)
	p.T = p.P / (R_mass * p.ρ)
end

@inbounds function find_pot_temp!(p::Particle, ρ0::Float64, T_bg::Float64, g::Float64, R_gas::Float64, R_mass::Float64)
	p.θ = p.T * (((T_bg * R_gas * ρ0) / p.P))^(2 / 7)
	p.θ_bg = background_pot_temperature(p.x[2], ρ0, T_bg, g, R_mass, R_gas)
	p.δθ = p.θ - p.θ_bg
end


# ==============
# Smoothing-length evolution (e.g.  setting the adaptive h)
# ==============

@inbounds function balance_of_smoothing!(p::Particle)
	if p.type == FLUID || p.type == INFLOW_INCOMING
		p.Dh = -0.5 * (p.h / p.ρ) * p.Dρ
	end
end

@inbounds function compute_smoothing!(p::Particle, dt::Float64, h_top::Float64)
	if p.type == FLUID || p.type == INFLOW_INCOMING
		p.h += p.Dh * dt
		p.h = min(p.h, h_top)
	end
end

@inbounds function reset_smoothing_rate!(p::Particle)
	p.Dh = 0.0
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
	return -g * VECY * p.δρ / p.ρ 

end

# ==============
# Momentum balance
# ==============
@inbounds function balance_of_momentum!(p::Particle, q::Particle, r::Float64, α::Float64, β::Float64, ϵ::Float64, rho_floor::Float64, γ::Float64)
	if p.type == FLUID
		x_pq = p.x - q.x
		v_pq = p.v - q.v
		dot_product = SmoothedParticles.dot(x_pq, v_pq)

		h_ij = 0.5 * (p.h + q.h)
		ker = rDwendland2(h_ij, r)

		prho = max(p.ρ, rho_floor)
		qrho = max(q.ρ, rho_floor)

		# pairwise conservative force
		p.Dv += -q.m * (p.δP / prho^2 + q.δP / qrho^2) * ker * x_pq

		# artificial viscous force
		if dot_product < 0.0
			c_i = sqrt(γ * p.P / prho)
			c_j = sqrt(γ * q.P / qrho)
			c_ij = 0.5 * (c_i + c_j)
			ρ_ij = 0.5 * (prho + qrho)
			μ_ij = (h_ij * dot_product) / (r * r + ϵ * h_ij * h_ij)
			π_ij = (-α * c_ij * μ_ij + β * μ_ij * μ_ij) / ρ_ij

			# artificial viscous force
			p.Dv += -q.m * π_ij * ker * x_pq
		end
	end

end

# ==============
# Mass balance 
# ==============

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID || p.type == INFLOW_INCOMING
		ker = (q.m / q.ρ) * rDwendland2(p.h, r)
		p.Dρ += p.ρ * ker * SmoothedParticles.dot(p.x - q.x, p.v - q.v)
	end
end

@inbounds function compute_density!(p::Particle, dt::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)
	if p.type == FLUID || p.type == INFLOW_INCOMING
		p.ρ += p.Dρ * dt
		p.ρ_bg = background_density(p.x[2], ρ0, T_bg, g, R_mass)
		p.δρ = p.ρ - p.ρ_bg
	end
end

@inbounds function reset_density_rate!(p::Particle)
	p.Dρ = 0.0
end

@inbounds function update_density_perturbation!(p::Particle, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)
	if (p.type == INFLOW_GHOST) || (p.type == OUTFLOW_GHOST) || (p.type == MOUNTAIN)
		# Background density is strictly a function of the ghost's fixed y-coordinate
		p.ρ_bg = background_density(p.x[2], ρ0, T_bg, g, R_mass)
		p.δρ = p.ρ - p.ρ_bg
	end
end


# ==============
# Move & accelerate
# ==============

function move!(p::Particle, dt::Float64)
	if p.type == FLUID || p.type == INFLOW_INCOMING
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
function verlet_step!(sys::ParticleSystem, global_params::Dict, sim_params::Dict )
	# unpack all parameters
	@unpack g, R_mass, cp, cv, γ, R_gas, T_bg, ρ0, N = global_params
	@unpack dom_height, dom_length, a, z_t = global_params
	@unpack h_m, z_β = sim_params
	@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel, v_initial = sim_params

	# compute derived parameters
	h0 = η * dr
	K = g / (R_mass * T_bg)
	c = sqrt(65e3 * (γ) / ρ0)
	dt = dt_rel * h0 / c
	γ_r = γ_r_rel * N
        h_top = η * dr * exp(K * dom_height / 2.0) * (4 / 3)^(1 / 4)
	

	# half-step acceleration & drift
	apply!(sys, p -> accelerate!(p, dt, g, z_t, z_β, γ_r))
	apply!(sys, p -> move!(p, dt))

	# lifecycle management: teleport particles at the outflow to the inflow and set correct values for them
	manage_particle_lifecycle!(sys, dr, K, dom_length)
	apply!(sys, p -> set_inflow_values!(p, v_initial, ρ0, T_bg, g, R_mass))

	# create the cell list after particles move and teleport
	create_cell_list!(sys)
	
	# extrapolate δρ, v to the ghost and mountain particles at the boundaries
	apply!(sys, p -> reset_ebc_gradients!(p))
	apply!(sys, (p, q, r) -> compute_ebc_gradients!(p, q, r))
	apply!(sys, p -> reset_ghost_ebc_search!(p))
	apply!(sys, (p, q, r) -> search_ghost_extrapolator!(p, q, dr, K))
	apply!(sys, p -> reset_mountain_ebc_search!(p))
	apply!(sys, (p, q, r) -> search_mountain_extrapolator!(p, q, r))
	apply_extrapolation!(sys, ρ0, T_bg, g, R_mass)
	apply_mountain_freeslip!(sys, h_m, a)

	# set correct density values to the ghost and mountain particles
	apply!(sys, p -> update_density_perturbation!(p, ρ0, T_bg, g, R_mass))

	# reset rates 
	apply!(sys, p -> reset_acceleration!(p))
	apply!(sys, p -> reset_density_rate!(p))
	apply!(sys, p -> reset_smoothing_rate!(p))

	# compute density (balance of mass) and smoothing length
	apply!(sys, (p, q, r) -> balance_of_mass!(p, q, r))
	apply!(sys, p -> balance_of_smoothing!(p))
	apply!(sys, p -> compute_density!(p, dt, ρ0, T_bg, g, R_mass))
	apply!(sys, p -> compute_smoothing!(p, dt, h_top))
	create_cell_list!(sys)

	# compute pressure
	apply!(sys, p -> compute_pressure!(p, ρ0, T_bg, g, R_mass, P_floor))

	# compute temperature and potential temperature
	apply!(sys, p -> find_temperature!(p, R_mass))
	apply!(sys, p -> find_pot_temp!(p, ρ0, T_bg, g, R_gas, R_mass))

	# compute acceleration (balance of momentum)
	apply!(sys, (p, q, r) -> balance_of_momentum!(p, q, r, α, β, ϵ, rho_floor, γ ))
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
	@unpack dom_height, dom_length, a, z_t = global_params
	@unpack h_m, z_β = sim_params
	@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel, v_initial = sim_params

	# compute derived parameters
	h0 = η * dr
	c = sqrt(65e3 * (γ) / ρ0)
	dt = dt_rel * h0 / c
	dt_frame = t_end / 100
	K = g / (R_mass * T_bg)
	γ_r = γ_r_rel * N

	# create the particle system
	sys = make_system(Particle, global_params, sim_params)

	# ==============
	# Initialization of the physical fields
	# ==============
	
	# initialization of the pressure
	apply!(sys, p -> compute_pressure!(p, ρ0, T_bg, g, R_mass, P_floor))

	# compute temperature and potential temperature
	apply!(sys, p -> find_temperature!(p, R_mass))
	apply!(sys, p -> find_pot_temp!(p, ρ0, T_bg, g, R_gas, R_mass))
	
	# extrapolate δρ, v to the ghost and mountain particles at the boundaries
	apply!(sys, p -> reset_ebc_gradients!(p))
	apply!(sys, (p, q, r) -> compute_ebc_gradients!(p, q, r))
	apply!(sys, p -> reset_ghost_ebc_search!(p))
	apply!(sys, (p, q, r) -> search_ghost_extrapolator!(p, q, dr, K))
	apply!(sys, p -> reset_mountain_ebc_search!(p))
	apply!(sys, (p, q, r) -> search_mountain_extrapolator!(p, q, r))
	apply_extrapolation!(sys, ρ0, T_bg, g, R_mass)
	apply_mountain_freeslip!(sys, h_m, a)

	# compute acceleration (balance of momentum)
	apply!(sys, p -> reset_acceleration!(p))
	apply!(sys, (p, q, r) -> balance_of_momentum!(p, q, r, α, β, ϵ, rho_floor, γ ))

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

