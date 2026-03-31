"""
Isentropic atmospheric flow around a mountain with the Witch of Agnesi profile:

 h(x) = (h_m a²) / (x² + a²).

"""

module EvolWBalancedHopkins

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

include(srcdir("core", "bc", "ebc_hopkins.jl"))

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
	θ_bg::Float64     # background potential temperature
	θ′::Float64       # potential temperature perturbation
	θ::Float64        # total potential temperature
	T_bg::Float64     # background temperature
	T::Float64        # total temperature
        A::Float64        # entropy-like variable
        A_bg::Float64     # entropy-like variable
	type::Float64     # particle type
	id::Int64              # the index of the particle in sys.particles
	spawn_y::Float64       # the altitude it was generated at
	grad_A::RealVector     # ∇A
	grad_u::RealVector     # ∇u (x-velocity gradient)
	grad_w::RealVector     # ∇w (y-velocity gradient)
	best_match_id::Int64   # rd of the fluid particle 'e'
	min_dist_sq::Float64   # argmin distance tracker

	function Particle(x::RealVector, v::RealVector, type::Float64, global_params::Dict, sim_params::Dict)
		# unpack all parameters
		@unpack g, R_mass, cp, cv, γ, R_gas, T_bg, ρ0, N = global_params
		@unpack dom_height, dom_length, a, z_t, z_β = global_params
		@unpack h_m = sim_params
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
			0.0,            # A
                        0.0,            # A_bg
			type,           # type
	                0,              # the index of the particle in sys.particles
			0.0,            #spawn_y::Float64
			VEC0,           #grad_A::RealVector
			VEC0,           #grad_u::RealVector
			VEC0,           #grad_w::RealVector
			0,              #best_match_id::Int64
			0.0,            #min_dist_sq::Float64
		)

		# initialization

		obj.ρ = 0.0 # needs SPH sum!
		obj.m = ρ0 * dr^2

		obj.T_bg = T_bg
		obj.T = T_bg

		obj.P_bg = 0.0 # needs SPH sum!
		obj.P = 0.0 # needs SPH sum!

		obj.A_bg = background_entropy(obj.x[2],ρ0, T_bg, g, R_mass, γ)
		obj.A = background_entropy(obj.x[2],ρ0, T_bg, g, R_mass, γ)

		obj.θ_bg = background_pot_temperature(obj.x[2],ρ0, T_bg, g, R_mass, R_gas)
		obj.θ′ = 0.0
		obj.θ = obj.θ′ + obj.θ_bg

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

function background_entropy(y::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64, γ::Float64)
	P_bg = background_pressure(y, ρ0, T_bg, g, R_mass)
	ρ_bg = background_density(y, ρ0, T_bg, g, R_mass)
	return P_bg / ρ_bg^γ
end


# ==============
# Pressure computation (e.g followed by sound speed computation)
# ==============

@inbounds function reset_pressure!(p::Particle, γ::Float64)
	ker_self = wendland2(p.h, 0.0)
        p.P = p.m * p.A^(1 / γ) * ker_self
        p.P_bg = p.m * p.A_bg^(1 / γ) * ker_self
end

@inbounds function compute_pressure!(p::Particle, q::Particle, r::Float64, γ::Float64)
        ker = wendland2(p.h, r)
        p.P += q.m * q.A^(1 / γ) * ker
        p.P_bg += q.m * q.A_bg^(1 / γ) * ker
end

@inbounds function finalize_pressure!(p::Particle, γ::Float64)
        p.P = p.P^γ
        p.P_bg = p.P_bg^γ
end

# ==============
# Thermodynamics (e.g determining temperature and potential temperature)
# ==============

@inbounds function compute_entropy!(p::Particle, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64, γ::Float64)
	p.A_bg = background_entropy(p.x[2], ρ0, T_bg, g, R_mass, γ)
end

@inbounds function find_temperature!(p::Particle, R_mass::Float64)
        p.T = p.P / (R_mass * p.ρ)
end

@inbounds function find_pot_temp!(p::Particle, ρ0::Float64, T_bg::Float64, g::Float64, R_gas::Float64, R_mass::Float64)
        p.θ = p.T * (((T_bg * R_gas * ρ0) / p.P))^(2 / 7)
	p.θ_bg = background_pot_temperature(p.x[2], ρ0, T_bg, g, R_mass, R_gas)
	p.θ′ = p.θ - p.θ_bg
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

@inbounds function balance_of_momentum!(p::Particle, q::Particle, r::Float64, α::Float64, β::Float64, ε::Float64, rho_floor::Float64, P_floor::Float64, γ::Float64)
	x_pq = p.x - q.x
	v_pq = p.v - q.v
	dot_product = SmoothedParticles.dot(x_pq, v_pq)

	prefac = q.m * (p.A * q.A)^(1 / γ)
	expfac = 1.0 - 2.0 / γ
	ker_i = rDwendland2(p.h, r)
	ker_j = rDwendland2(q.h, r)
	f_i = 1.0 / p.Omega
	f_j = 1.0 / q.Omega
	pP = max(P_floor, p.P)
	qP = max(P_floor, q.P)

	# acceleration due the gradient of total pressure
	a_tot = -prefac * (f_i * pP^expfac * ker_i + f_j * qP^expfac * ker_j) * x_pq

	prefac_bg = q.m * (p.A_bg * q.A_bg)^(1 / γ)
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
	if p.type == FLUID || p.type == INFLOW_INCOMING
		p.x += dt * p.v
	end
end

function accelerate!(p::Particle, dt::Float64, z_t::Float64, z_β::Float64, γ_r::Float64)
	if p.type == FLUID
		p.v += 0.5 * dt * (p.Dv + damping_structure(p.x[2], p.v, z_t, z_β, γ_r)) # this is a vector sum
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
	@unpack dom_height, dom_length, a, z_t, z_β = global_params
	@unpack h_m = sim_params
	@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel, v_initial = sim_params

	# compute derived parameter
	h0 = η * dr
	c = sqrt(65e3 * (γ) / ρ0)
	dt = dt_rel * h0 / c
	γ_r = γ_r_rel * N
	K = g / (R_mass * T_bg)
	
	# half-step acceleration & drift
	apply!(sys, p -> accelerate!(p, dt, z_t, z_β, γ_r))
	apply!(sys, p -> move!(p, dt))
	
	# lifecycle management: teleport particles at the outflow to the inflow and set correct values for them
	manage_particle_lifecycle!(sys, dr, K, dom_length)
	apply!(sys, p -> set_inflow_values!(p, v_initial, ρ0, T_bg, g, R_mass, γ))

	# create the cell list after particles move and teleport
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
	
	# extrapolate A, v to the ghost particles at the boundaries
	apply!(sys, p -> reset_ebc_gradients!(p))
	apply!(sys, (p, q, r) -> compute_ebc_gradients!(p, q, r))
	apply!(sys, p -> reset_ghost_ebc_search!(p))
	apply!(sys, (p, q, r) -> search_ghost_extrapolator!(p, q, dr, K))
	apply!(sys, p -> reset_mountain_ebc_search!(p))
	apply!(sys, (p, q, r) -> search_mountain_extrapolator!(p, q, r))
	apply_extrapolation!(sys)
	apply_mountain_freeslip!(sys, h_m, a)

	# compute pressure
	apply!(sys, p -> reset_pressure!(p, γ))
	apply!(sys, p -> compute_entropy!(p, ρ0, T_bg, g, R_mass, γ))
	apply!(sys, (p, q, r) -> compute_pressure!(p, q, r, γ))
	apply!(sys, p -> finalize_pressure!(p, γ))

	# compute temperature and potential temperature
	apply!(sys, p -> find_temperature!(p, R_mass))
	apply!(sys, p -> find_pot_temp!(p, ρ0, T_bg, g, R_gas, R_mass))

	# compute the forces
	apply!(sys, p -> reset_acceleration!(p))
	apply!(sys, (p, q, r) -> balance_of_momentum!(p, q, r, α, β, ϵ, rho_floor, P_floor, γ))
	apply!(sys, p -> accelerate!(p, dt, z_t, z_β, γ_r))
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
	@unpack dom_height, dom_length, a, z_t, z_β = global_params
	@unpack h_m = sim_params
	@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel = sim_params

	# compute derived parameters
	h0 = η * dr
	c = sqrt(65e3 * (γ) / ρ0)
	dt = dt_rel * h0 / c
	dt_frame = t_end / 100

	# create the particle system
	sys = make_system(Particle, global_params, sim_params)

	# ==============
	# Initialization of the physical fields
	# ==============
	
	# set the initial velocity to the whole bulk
	#apply!(sys, p -> set_bulk_velocity!(p, v_initial))
	
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
	apply!(sys, p -> reset_pressure!(p, γ))
	apply!(sys, p -> compute_entropy!(p, ρ0, T_bg, g, R_mass, γ))
	apply!(sys, (p, q, r) -> compute_pressure!(p, q, r, γ))
	apply!(sys, p -> finalize_pressure!(p, γ))

	# compute temperature and potential temperature
	apply!(sys, p -> find_temperature!(p, R_mass))
	apply!(sys, p -> find_pot_temp!(p, ρ0, T_bg, g, R_gas, R_mass))

	# compute acceleration (balance of momentum)
	apply!(sys, p -> reset_acceleration!(p))
	apply!(sys, (p, q, r) -> balance_of_momentum!(p, q, r, α, β, ϵ, rho_floor, P_floor, γ ))

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
