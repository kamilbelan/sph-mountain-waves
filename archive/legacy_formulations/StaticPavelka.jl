"""
 Static atmosphere above a mountain with the Witch of Agnesi profile:

h(x)=(h_ma²)/(x²+a²),

all thermodynamic proccesses are adiabatic
"""

module StaticPavelka

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

# ============================================================
# PARTICLE CONSTRUCTOR
# ============================================================

mutable struct Particle <: AbstractParticle
        h::Float64 #smoothing length
        Dh::Float64 #rate of the smoothing length
        x::RealVector  #position
        m::Float64     #mass
        v::RealVector  #velocity
        Dv::RealVector #acceleration
        ρ::Float64    #density
        Dρ::Float64    #rate of density
        P::Float64     #pressure
        θ::Float64     #potential temperature
        S::Float64     #entropy
        s::Float64     #entropy density
        T::Float64     #temperature 
        type::Float64  #particle type

        function Particle(x::RealVector, u::RealVector, type::Float64, global_params::Dict, sim_params::Dict)
		# unpack all parameters
		@unpack g, R_mass, cp, cv, γ, R_gas, T_bg, ρ0, N = global_params
		@unpack dom_height, dom_length, h_m, a, z_t, z_β = global_params
		@unpack rho_floor, P_floor, ϵ, α, β  = sim_params
		@unpack η, dr, dt_rel, t_end = sim_params

		# Derived values within the constructor
		h0 = η * dr
                obj = new(
			h0,
			0.0,
			x,
			0.0,
			u,
			VEC0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			type,
		)
                obj.T = T_bg # prescribe initial temperature
                obj.ρ   = ρ0 * exp(-obj.x[2] * g / (R_mass * obj.T))
                obj.m = obj.ρ * dr^2
                obj.P = R_mass * obj.T * obj.ρ
                obj.θ = obj.T * (((T_bg * R_gas * ρ0) / obj.P))^(2 / 7) # 2/7 is exactly R_gas/cp
                obj.S = obj.m * cv * log((cv * obj.T * (γ - 1)) / (γ * obj.ρ^(γ - 1))) # calculate initial entropy from the initial temperature
                return obj
        end
end


# ==============
# Pressure computation
# ==============

@inbounds function set_pressure!(p::Particle, R_mass::Float64, P_floor::Float64)
        if p.type == FLUID
                pP = R_mass * p.ρ * p.T
		p.P = max(pP, P_floor)
        end
end


# ==============
# Thermodynamics 
# ==============

@inbounds function find_s!(p::Particle)
        if p.type == FLUID
                p.s = p.S * p.ρ / p.m
        end
end

@inbounds function set_temperature!(p::Particle, γ::Float64, cv::Float64)
        if p.type == FLUID
                p.T = (p.ρ^(γ - 1.0)) * exp(p.s / (p.ρ * cv)) / (cv * (γ - 1.0))
        end
end

@inbounds function find_pot_temp!(p::Particle, ρ0::Float64, T_bg::Float64, R_gas::Float64)
        if p.type == FLUID
                p.θ = p.T * (((T_bg * R_gas * ρ0) / p.P)^(2))^(1 / 7) #2/7 is exactly R_gas/cp
        end
end

@inbounds function entropy_production!(p::Particle, q::Particle, r::Float64, dt::Float64, μ::Float64)
        if p.type == FLUID && q.type == FLUID
                ker = (q.m / q.ρ) * rDwendland2(0.5 * (p.h + q.h), r)
                x_pq = p.x - q.x
                u_pq = p.v - q.v
                p.S += -4.0 * p.m * q.m * p.ρ * ker * μ / (p.T * p.ρ * q.ρ) * SmoothedParticles.dot(u_pq, x_pq)^2 / (r * r + 0.01 * p.h * q.h) * dt #viscous
        end
end

# ==============
# Smoothing-length & density evolution
# ==============

@inbounds function balance_of_smoothing!(p::Particle)
        p.Dh += -0.5 * (p.h / p.ρ) * p.Dρ
end

@inbounds function update_smoothing!(p::Particle, dt::Float64)
        if p.type == FLUID
                p.h += dt * p.Dh
        end
        p.Dh = 0.0
end

@inbounds function update_density!(p::Particle, dt::Float64, rho_floor::Float64)
        if p.type == FLUID
                p.ρ += dt * p.Dρ
		p.ρ = max(p.ρ, rho_floor)
        end
        p.Dρ = 0.0
end

# ==============
# Momentum balance
# ==============

@inbounds function balance_of_momentum!(p::Particle, q::Particle, r::Float64, μ::Float64)
        ker = (q.m / q.ρ) * rDwendland2(0.5 * (p.h + q.h), r)
        x_pq = p.x - q.x
        p.Dv += -p.ρ * ker * (p.P / p.ρ^2 + q.P / q.ρ^2) * x_pq
        p.Dv += p.ρ * 8.0 * ker * μ::Float64 / (p.ρ * q.ρ) * SmoothedParticles.dot(p.v - q.v, x_pq) / (r * r + 0.0025 * (p.h + q.h)^2) * x_pq
end


# ==============
# Mass balance
# ==============

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64, ν::Float64)
        ker = (q.m / q.ρ) * rDwendland2(0.5 * (p.h + q.h), r)
        p.Dρ += p.ρ * ker * (SmoothedParticles.dot(p.x - q.x, p.v - q.v))
        if p.type == FLUID && q.type == FLUID
                p.Dρ += 2 * ν / p.ρ * (p.ρ - q.ρ)
        end
end


# ==============
# Rayleigh damping
# ==============

function damping_structure(z::Float64, z_t::Float64, z_β::Float64, γ_r::Float64)
	if z >= (z_t - z_β)
		return -γ_r * (sin(π / 2 * (1 - (z_t - z_β) / z_β)))^2 * VECY
	else
		return VEC0
	end
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
		p.v += 0.5 * dt * (p.Dv -g*VECY + damping_structure(p.x[2], z_t, z_β, γ_r)) # this is a vector sum
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


function verlet_step!(sys::ParticleSystem)
        apply!(sys, accelerate!)
        apply!(sys, move!)
        create_cell_list!(sys)

        apply!(sys, balance_of_mass!)
        apply!(sys, balance_of_smoothing!)
        apply!(sys, update_smoothing!)
        apply!(sys, update_density!)
        create_cell_list!(sys)
        apply!(sys, find_s!)
        apply!(sys, set_temperature!)
        apply!(sys, set_pressure!)
        apply!(sys, entropy_production!)
        apply!(sys, balance_of_momentum!)
        apply!(sys, accelerate!)
end

"""
Computes the average velocity of the particles
"""

function avg_velocity(sys::ParticleSystem)::Float64
        sum = 0.0
        for p in sys.particles
                sum += norm(p.u)
        end
        avg = sum / length(sys.particles)
        return avg
end

"""
Finds the largest velocity of the particles
"""

function max_velocity(sys::ParticleSystem)::Float64
        velocities = [norm(p.u) for p in sys.particles]
        max_vel = maximum(velocities)
        return max_vel
end

function energy(sys::ParticleSystem)::Float64
        E = 0.0
        for p in sys.particles
		E += 0.5 * p.ρ * SmoothedParticles.dot(p.u,p.u) + p.ρ * g * p.x[2]
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

	# verlet step 
	function verlet_step!(sys)
		# half-step acceleration & drift
		apply!(sys, p -> accelerate!(p, dt, g, z_t, z_β, γ_r))
		apply!(sys, p -> move!(p, dt))
		create_cell_list!(sys)

		# compute density and smoothing length
		apply!(sys, (p, q, r) -> balance_of_mass!(p, q, r, ν))
		apply!(sys, balance_of_smoothing!)
		apply!(sys, p -> update_smoothing!(p, dt))
		apply!(sys, p -> update_density!(p, dt, rho_floor))
		create_cell_list!(sys)

		# thermodynamics: dissipate, find tmep
		apply!(sys, find_s!)
		apply!(sys, p -> set_temperature!(p, γ, cv))
		apply!(sys, p -> set_pressure!(p, R_mass, P_floor))
		apply!(sys, (p,q,r) -> entropy_production!(p, q, r, dt, μ))
		apply!(sys, (p,q,r) -> balance_of_momentum!(p, q, r, μ))
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
	metadata_dict["module"] = "StaticPavelka" # add model metadata manually
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

	println("\n End of the road 🏵. Simulation $run_name complete. ")
	println("\n RESULTS: $run_dir")
	println("\n PLOTS: ", plotsdir("diagnostics_$(run_name).png"))
end
end # module
