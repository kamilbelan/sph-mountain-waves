"""
 Static atmosphere above a mountain with the Witch of Agnesi profile:

 h(x) = (h_m a²) / (x² + a²),

 all thermodynamic processes are adiabatic
"""

module StaticStandardWCSPH

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
#include(srcdir("utils", "new_packing.jl"))
include(srcdir("utils", "make_xdmf.jl"))

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
		@unpack dom_height, dom_length, h_m, a, z_t, z_β = global_params
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
		obj.ρ_bg = background_density(obj.x[2], ρ0, T_bg, g, R_mass)
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

		obj.m = obj.ρ * dr^2
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

@inbounds function compute_pressure!(p::Particle, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64, P_floor)
	p.P_bg = background_pressure(p.x[2], ρ0, T_bg, g, R_mass)
	p.P′ = p.c^2 * p.ρ′
	pP = p.P_bg + p.P′
	p.P = max(pP, P_floor)
end


# ==============
# Thermodynamics (e.g determning temperature and potential temperature)
# ==============

@inbounds function find_temperature!(p::Particle, R_mass::Float64)
	p.T = p.P / (R_mass * p.ρ)
	p.T′ = p.T - p.T_bg
end

@inbounds function find_pot_temp!(p::Particle, ρ0::Float64, T_bg::Float64, g::Float64, R_gas::Float64, R_mass::Float64)
	p.θ = p.T * (((T_bg * R_gas * ρ0) / p.P))^(2 / 7)
	p.θ_bg = background_pot_temperature(p.x[2], ρ0, T_bg, g, R_mass, R_gas)
	p.θ′ = p.θ - p.θ_bg
end


# ==============
# Smoothing-length evolution (e.g.  setting the adaptive h)
# ==============

@inbounds function balance_of_smoothing!(p::Particle)
	p.Dh = -0.5 * (p.h / p.ρ) * p.Dρ
end

@inbounds function compute_smoothing!(p::Particle, dt::Float64)
	p.h += p.Dh * dt

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
	return -g * VECY * p.ρ′ / p.ρ # the (density) of gravity is - g * VECY

end


# ==============
# Momentum balance
# ==============
@inbounds function balance_of_momentum!(p::Particle, q::Particle, r::Float64, α::Float64, β::Float64, ϵ::Float64, rho_floor::Float64, γ::Float64)
	x_pq = p.x - q.x
	v_pq = p.v - q.v
	dot_product = SmoothedParticles.dot(x_pq, v_pq)

	h_ij = 0.5 * (p.h + q.h)
	ker = rDwendland2(h_ij, r)

	prho = max(p.ρ, rho_floor)
	qrho = max(q.ρ, rho_floor)

	# pairwise conservative force
	p.Dv += -q.m * (p.P′ / prho^2 + q.P′ / qrho^2) * ker * x_pq

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

# ==============
# Mass balance 
# ==============

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	ker = (q.m / q.ρ) * rDwendland2(p.h, r)
	p.Dρ += p.ρ * ker * SmoothedParticles.dot(p.x - q.x, p.v - q.v)
end

@inbounds function compute_density!(p::Particle, dt::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)
	p.ρ += p.Dρ * dt
	p.ρ_bg = background_density(p.x[2], ρ0, T_bg, g, R_mass)
	p.ρ′ = p.ρ - p.ρ_bg
end

@inbounds function reset_density_rate!(p::Particle)
	p.Dρ = 0.0
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
# Diagnostics
# ==============

function avg_velocity(sys::ParticleSystem)::Float64
	v = 0.0
	for p in sys.particles
		v += SmoothedParticles.norm(p.v)
	end
	v = v / length(sys.particles)
	return v
end

function max_velocity(sys::ParticleSystem)::Float64
	v = maximum(SmoothedParticles.norm(p.v) for p in sys.particles)
	return v
end

function energy(sys::ParticleSystem, g::Float64)::Float64
	E = 0.0
	for p in sys.particles
		E += 0.5 * p.ρ * SmoothedParticles.dot(p.v, p.v) + p.ρ * g * p.x[2]
	end
	return E
end


# ==============
# Main Entry Point
# ==============

function run_sim(global_params::Dict, sim_params::Dict)
	# unpack all parameters
	@unpack g, R_mass, cp, cv, γ, R_gas, T_bg, ρ0, N = global_params
	@unpack dom_height, dom_length, h_m, a, z_t, z_β = global_params
	@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel = sim_params

	# compute derived parameters
	h0 = η * dr
	bc_width = 6*dr
	c = sqrt(65e3 * (γ) / ρ0)
	dt = dt_rel * h0 / c
	dt_frame = t_end / 100
	γ_r = γ_r_rel * N

	# system construction 
	function make_system()
		grid = Grid(dr, :hexagonal)
		domain = Rectangle(-dom_length / 2.0, 0.0, dom_length / 2.0, dom_height)
		fence = BoundaryLayer(domain, grid, bc_width)
		witch_profile(x) = (h_m * a^2) / (x^2 + a^2)
		mountain = Specification(domain, x -> (x[2] <= witch_profile(x[1])))

		sys = ParticleSystem(Particle, domain + fence, h0)
		# passing params to the Particle constructor
		generate_particles!(sys, grid, domain - mountain, x -> Particle(x, VEC0, FLUID, global_params, sim_params))
		generate_particles!(sys, grid, fence, x -> Particle(x, VEC0, WALL, global_params, sim_params))
		generate_particles!(sys, grid, mountain, x -> Particle(x, VEC0, FLUID, global_params, sim_params))

		create_cell_list!(sys)
		return sys
	end

	# ==============
	# Verlet step
	# ==============
	function verlet_step!(sys)
		# half-step acceleration & drift
		apply!(sys, p -> accelerate!(p, dt, g, z_t, z_β, γ_r))
		apply!(sys, p -> move!(p, dt))
		create_cell_list!(sys)

		# reset rates 
		apply!(sys, p -> reset_acceleration!(p))
		apply!(sys, p -> reset_density_rate!(p))
		apply!(sys, p -> reset_smoothing_rate!(p))

		# compute density (balance of mass) and smoothing length
		apply!(sys, (p, q, r) -> balance_of_mass!(p, q, r))
		apply!(sys, p -> balance_of_smoothing!(p))
		apply!(sys, p -> compute_density!(p, dt, ρ0, T_bg, g, R_mass))
		apply!(sys, p -> compute_smoothing!(p, dt))
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

	# execution loop
	sys = make_system()

	# ==============
	# Initialization
	# ==============
	
	# initialization of the pressure
	apply!(sys, p -> compute_pressure!(p, ρ0, T_bg, g, R_mass, P_floor))

	# compute temperature and potential temperature
	apply!(sys, p -> find_temperature!(p, R_mass))
	apply!(sys, p -> find_pot_temp!(p, ρ0, T_bg, g, R_gas, R_mass))

	# compute acceleration (balance of momentum)
	apply!(sys, p -> reset_acceleration!(p))
	apply!(sys, (p, q, r) -> balance_of_momentum!(p, q, r, α, β, ϵ, rho_floor, γ ))

	# ==============
	# Output handling
	# ==============
	
	# choose the parameters for the folder name
	name_keys = [:dr, :dt_rel, :t_end, :gamma_r_rel]
	name_params = filter(p -> p.first in name_keys, sim_params)
	short_name = savename(name_params)

	# extract model name safely 
	model_name_str = string(sim_params[:model])

	# generate a hex hash for uniqueness
	full_hash = string(hash(sim_params), base=16)[1:6]
	run_name = "$(short_name)_$(full_hash)"

	# build the directory tree: data/sims/modelname/YYYY-MM-DD/run_name
	date_str = Dates.format(now(), "yyyy-mm-dd")
	run_dir = datadir("sims", model_name_str, date_str, run_name)
	mkpath(run_dir)

	println("Output directory: $run_dir")

	# save metadata in the run directory
	metadata_dict = Dict(String(k) => v for (k, v) in sim_params)
	metadata_dict["module"] = model_name_str 
	@tagsave(joinpath(run_dir, "metadata.jld2"), metadata_dict)

	# initialize data accumulation for diagnostics
	average_velocities = DataFrame(t=Float64[], u=Float64[])
	maximum_velocities = DataFrame(t=Float64[], u=Float64[])
	ene = DataFrame(t=Float64[], E=Float64[])

	# ==============
	# Time loop
	# ==============
	nsteps = Int(round(t_end / dt))
	frame_counter = 0

	println("---------------------------")

	for k = 1:nsteps
		t = k * dt
		verlet_step!(sys)

		save_interval = max(1, Int(round(dt_frame / dt)))
		if k % save_interval == 0
			@show t
			println("num. of particles = ", length(sys.particles))

			u_avg = avg_velocity(sys)
			@show u_avg
			push!(average_velocities, (t, u_avg))

			u_max = max_velocity(sys)
			@show u_max
			push!(maximum_velocities, (t, u_max))

			E = energy(sys, g)
			@show E
			push!(ene, (t, E))

			frame_file = joinpath(run_dir, "frame_$(lpad(frame_counter, 4, '0')).h5")

			# pre-allocate flat arrays
			N_parts = length(sys.particles)
			pos_matrix = zeros(Float64, 3, N_parts)
			vel_matrix = zeros(Float64, 3, N_parts)
			densities = zeros(Float64, N_parts)
			densities_pert = zeros(Float64, N_parts)
			pressures = zeros(Float64, N_parts)
			pressures_pert = zeros(Float64, N_parts)
			temperatures = zeros(Float64, N_parts)
			pot_temperatures = zeros(Float64, N_parts)
			types = zeros(Float64, N_parts)

			@inbounds for (i, p) in enumerate(sys.particles)
				pos_matrix[1, i] = p.x[1]
				pos_matrix[2, i] = p.x[2]
				pos_matrix[3, i] = 0.0  

				vel_matrix[1, i] = p.v[1]
				vel_matrix[2, i] = p.v[2]
				vel_matrix[3, i] = 0.0  

				densities[i] = p.ρ
				densities_pert[i] = p.ρ′
				pressures[i] = p.P
				pressures_pert[i] = p.P′
				temperatures[i] = p.T
				pot_temperatures[i] = p.θ
				types[i] = p.type
			end

			#  write securely to HDF5
			h5open(frame_file, "w") do file
				# save metadata as HDF5 attributes 
				attributes(file)["time"] = t
				attributes(file)["frame_counter"] = frame_counter
				attributes(file)["n_particles"] = N_parts

				# save the arrays as primary datasets
				file["positions"] = pos_matrix
				file["velocities"] = vel_matrix
				file["densities"] = densities
				file["densities_pert"] = densities_pert
				file["pressures"] = pressures
				file["pressures_pert"] = pressures_pert
				file["temperatures"] = temperatures
				file["pot_temperatures"] = pot_temperatures
				file["types"] = types
			end
			frame_counter += 1
		end
	end

	# create a xdmf file to use with the generated h5 files
	generate_sph_xdmf(run_dir)

	# write CSVs once at the end
	CSV.write(joinpath(run_dir, "average_velocities.csv"), average_velocities)
	CSV.write(joinpath(run_dir, "maximum_velocities.csv"), maximum_velocities)
	CSV.write(joinpath(run_dir, "energies.csv"), ene)

	println("Generating diagnostic plots...")

	# fallback to headless plotting if not in an interactive REPL (eg for cluster use)
	if !isinteractive()
		ENV["GKSwstype"] = "nul"
	end

	p1 = plot(average_velocities.t, average_velocities.u; xlabel="t (s)", ylabel="avg. velocity (m/s)", lc=:blue)
	p2 = plot(maximum_velocities.t, maximum_velocities.u; xlabel="t (s)", ylabel="max. velocity (m/s)", lc=:purple)
	p3 = plot(ene.t, ene.E; xlabel="t (s)", ylabel="Total energy J", lc=:orange)

	plt = plot(p1, p2, p3; layout=(3, 1), size=(800, 900))

	# display the plot if local, otherwise skip displaying on HPC
	if isinteractive()
		display(plt)
	end

	# save the plot into the run directory
	savefig(plt, joinpath(run_dir, "diagnostics.pdf"))

	println("\nEnd of the road 🏵.")
	println("Simulation $run_name complete.")
	println("All artifacts (data, plots, metadata) securely saved to:\n  $run_dir")

	return run_dir

end
end # module

