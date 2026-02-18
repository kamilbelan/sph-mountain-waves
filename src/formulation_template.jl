"""
 Static atmosphere above a mountain with the Witch of Agnesi profile:

 h(x) = (h_m a²) / (x² + a²),

 all thermodynamic processes are adiabatic
"""

module StaticTEMPLATE

export run_sim

using DrWatson
@quickactivate "SPH"

using Printf
using SmoothedParticles
using DataFrames
using Plots
using JLD2
using CSV
using Parameters 
using LinearAlgebra

const FLUID = 0.0
const WALL = 1.0
const MOUNTAIN = 2.0
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
	ρ_bg::Float64     # background density
	ρ′::Float64       # density perturbation
	ρ::Float64        # total density
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
			x,              # x 	 	
			0.0,            # m
			v,              # v
			VEC0,           # Dv
			0.0,            # ρ_bg
			0.0,            # ρ′
			0.0,            # ρ
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


# ==============
# Pressure computation (e.g followed by sound speed computation)
# ==============


# ==============
# Thermodynamics (e.g determning temperature and potential temperature)
# ==============

# ==============
# Smoothing-length evolution (e.g. calculating the smth. rate, setting the adaptive h)
# ==============

# ==============
# Additional forces: Rayleigh damping & gravity/buyoancy
# ==============

function damping_structure(z::Float64, z_t::Float64, z_β::Float64, γ_r::Float64)
	if z >= (z_t - z_β)
		return -γ_r * (sin(π / 2 * (1 - (z_t - z_β) / z_β)))^2 * VECY
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
		p.v += 0.5 * dt * (p.Dv + buyoancy_force(p, g) + damping_structure(p.x[2], z_t, z_β, γ_r)) # this is a vector sum
	end
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

		# compute density and smoothing length
		apply!(sys, p -> reset_density!(p))
		apply!(sys, compute_density!)
		apply!(sys, p -> finalize_density!(p, ρ0, T_bg, g, R_mass))
		apply!(sys, p -> update_smoothing!(p, η, rho_floor))
		create_cell_list!(sys)

		# pressure–entropy: build P̄ from A
		apply!(sys, p -> compute_pressure!(p, ρ0, T_bg, g, R_mass, P_floor))
		apply!(sys, p -> compute_sound_speed!(p, rho_floor, γ))

		# thermodynamics from P̄ and ρ
		apply!(sys, p -> find_temperature!(p, R_mass))
		apply!(sys, p -> find_pot_temp!(p, ρ0, T_bg, g, R_gas, R_mass))

		# forces
		apply!(sys, (p, q, r) -> balance_of_momentum!(p, q, r, α, β, ϵ, rho_floor, γ ))
		apply!(sys, p -> accelerate!(p, dt, g, z_t, z_β, γ_r))
	end

	# execution loop
	sys = make_system()

	run_name = savename(sim_params) # DrWatson magic 
	run_dir = datadir("sims", run_name)
	mkpath(run_dir)

	println("Output directory: $run_dir")

	# a great help from DrWatson with metadata processing
	metadata_dict = Dict(String(k) => v for (k, v) in sim_params)
	metadata_dict["module"] = "StaticWCSPH" # add module metadata manually
	@tagsave(joinpath(run_dir, "metadata.jld2"), metadata_dict)

	# UNCOMMENT TO SAVE VTK FRAMES FOR PARAVIEW!
	#outpath = joinpath(RESULTS_DIR, folder_name)
	#out = new_pvd_file(outpath)
	#save_frame!(out, sys, export_vars...)

	average_velocities = DataFrame(t=Float64[], u=Float64[])
	maximum_velocities = DataFrame(t=Float64[], u=Float64[])
	ene = DataFrame(t=Float64[], E=Float64[])

	nsteps = Int(round(t_end / dt))
	frame_counter = 0

	@show T_bg
	@show ρ0
	@show c
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

			# the output is saved into a jld2 file, which is more suitable for analysis (in julia) than vtk
			frame_file = joinpath(run_dir, "frame_$(lpad(frame_counter,4,'0')).jld2")
			jldsave(frame_file;
	   t, frame_counter,
	   n_particles=length(sys.particles),
	   positions=[p.x for p in sys.particles],
	   velocities=[p.v for p in sys.particles],
	   densities=[p.ρ for p in sys.particles],
	   densities_pert=[p.ρ′ for p in sys.particles],
	   pressures=[p.P for p in sys.particles],
	   pressures_pert=[p.P′ for p in sys.particles],
	   temperatures=[p.T for p in sys.particles],
	   pot_temperatures=[p.θ for p in sys.particles],
	   types=[p.type for p in sys.particles]
	   )

			# UNCOMMENT TO SAVE VTK FRAMES FOR PARAVIEW!
			#save_frame!(out, sys, export_vars...)
			CSV.write(joinpath(run_dir, "average_velocities.csv"), average_velocities)
			CSV.write(joinpath(run_dir, "maximum_velocities.csv"), maximum_velocities)
			CSV.write(joinpath(run_dir, "energies.csv"), ene)
		end
	end

	#save_pvd_file(out)

	p1 = plot(
		average_velocities.t, average_velocities.u;
		xlabel="t (s)",
		ylabel="avg. velocity (m/s)",
		lc=:blue,
	)
	p2 = plot(
		maximum_velocities.t, maximum_velocities.u;
		xlabel="t (s)",
		ylabel="max. velocity (m/s)",
		lc=:purple,
	)
	p3 = plot(
		ene.t, ene.E;
		xlabel="t (s)",
		ylabel="Total energy J",
		lc=:orange,)

	#savefig(joinpath(outpath,"velocities.pdf"))

	plot(p1, p2, p3; layout=(3, 1), size=(800, 900))

	# save the plots for immediate diagnostics
	savefig(plotsdir("diagnostics_$(run_name).pdf"))

	println("\n End of the road 🏵. \n Simulation $run_name complete. ")
	println("\n RESULTS: $run_dir")
	println("\n PLOTS: ", plotsdir("diagnostics_$(run_name).png"))
end
end # module

