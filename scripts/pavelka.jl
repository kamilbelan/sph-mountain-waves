"""
 Static atmosphere above a mountain with the Witch of Agnesi profile:

h(x)=(hₘa²)/(x²+a²),

all thermodynamic proccesses are adiabatic
"""

module AdiabaticStaticWitch
export main

using DrWatson
@quickactivate "SPH"
using Printf
using SmoothedParticles
using DataFrames
using Plots
using JLD2
using CSV


# UNCOMMENT IF ONE WISHES TO SAVE FRAMES FOR PARAVIEW
#const folder_name = "pavelka_total_witch"
#const export_vars = (:u, :ρ, :P, :θ, :T, :type)

#=
Declare constants
=#

#geometry parameters
const dom_height = 26e3   #height of the domain 
const dom_length = 400e3  #length of the domain
const dr = dom_height / 75  #average particle distance (decrease to make finer simulation)
const h0 = 1.8 * dr          #smoothing length    
const bc_width = 6 * dr     #boundary width
const hₘ = 0#100            #parameters for the Witch of Agnesi profile; mountain height
const a = 0#10e3            #parameters for the Witch of Agnesi profile; mountain width

#physical parameters
const ρ0 = 1.393 #referential fluid density
const mu = 1.0# 15.98e-6 #dynamic viscosity
const c = sqrt(65e3 * (7 / 5) / ρ0) #speed of sound
const nu = 0.1 * h0 * c

#meteorological parameters
const N = sqrt(0.0196)     #Brunt-Vaisala frequency
const g = 9.81             #gravity
const R_mass = 287.05      #specific molar gas constant
const γᵣ = 10 * N            #damping coefficient
const zᵦ = 12e3             #bottom part of the damping layer
const zₜ = dom_height      #top part of the damping layer

#thermodynamical parameters
const R_gas = 8.314        #universal molar gas constant
const cp = 7 * R_mass / 2      #specific molar heat capacity at a constant pressure
const cv = cp - R_mass #specific molar heat capacity at a constant volume
const γ = cp / cv #poisson constant
const T0 = 250 #initial temperature

#temporal parameters
const dt = 0.01 * h0 / c   #time step
const t_end = 20 #end of simulation
const dt_frame = t_end / 100 #how often data is saved

#particle types
const FLUID = 0.0
const WALL = 1.0
const MOUNTAIN = 2.0


include(scriptsdir("utils", "new_packing.jl"))
"""
Declare the struct Particle <: AbstractParticle
"""
mutable struct Particle <: AbstractParticle
        h::Float64 #smoothing length
        Dh::Float64 #rate of the smoothing length
        x::RealVector  #position
        m::Float64     #mass
        u::RealVector  #velocity
        Du::RealVector #acceleration
        ρ::Float64    #density
        Dρ::Float64    #rate of density
        P::Float64     #pressure
        θ::Float64     #potential temperature
        S::Float64     #entropy
        s::Float64     #entropy density
        T::Float64     #temperature 
        type::Float64  #particle type

        function Particle(x::RealVector, u::RealVector, type::Float64)
                obj = new(h0, 0.0, x, 0.0, u, VEC0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, type)
                obj.T = T0 #prescribe initial temperature
                obj.ρ   = ρ0 * exp(-obj.x[2] * g / (R_mass * obj.T))
                obj.m = obj.ρ * dr^2
                obj.P = R_mass * obj.T * obj.ρ
                obj.θ = obj.T * (((T0 * R_gas * ρ0) / obj.P))^(2 / 7) #2/7 is exactly R_gas/cp
                obj.S = obj.m * cv * log((cv * obj.T * (γ - 1)) / (γ * obj.ρ^(γ - 1))) #calculate initial entropy from the initial temperature
                return obj
        end
end

"""
Define geometry and make particles
"""

function make_system()
        grid = Grid(dr, :hexagonal)
        domain = Rectangle(-dom_length / 2.0, 0.0, dom_length / 2.0, dom_height)
        fence = BoundaryLayer(domain, grid, bc_width)

        witch_profile(x) = (hₘ * a^2) / (x^2 + a^2)
        mountain = Specification(domain, x -> (x[2] <= witch_profile(x[1])))

        sys = ParticleSystem(Particle, domain + fence, h0)
        generate_particles!(sys, grid, domain - mountain, x -> Particle(x, VEC0, FLUID))
        generate_particles!(sys, grid, fence, x -> Particle(x, VEC0, WALL))
        generate_particles!(sys, grid, mountain, x -> Particle(x, VEC0, FLUID))

        create_cell_list!(sys)
        create_cell_list!(sys)
        apply!(sys, balance_of_mass!)
        apply!(sys, balance_of_smoothing!)
        apply!(sys, find_s!)
        apply!(sys, set_temperature!)
        apply!(sys, set_pressure!)
        apply!(sys, balance_of_momentum!)
        return sys
end

"""
Define particle interactions 	
"""

@inbounds function balance_of_momentum!(p::Particle, q::Particle, r::Float64)
        ker = (q.m / q.ρ) * rDwendland2(0.5 * (p.h + q.h), r)
        x_pq = p.x - q.x
        p.Du += -p.ρ * ker * (p.P / p.ρ^2 + q.P / q.ρ^2) * x_pq
        p.Du += p.ρ * 8.0 * ker * mu / (p.ρ * q.ρ) * dot(p.u - q.u, x_pq) / (r * r + 0.0025 * (p.h + q.h)^2) * x_pq
end

"""
Calculate density, pressure, temperature, potential temperature, entropy, entropy density
"""


@inbounds function find_s!(p::Particle)
        if p.type == FLUID
                p.s = p.S * p.ρ / p.m
        end
end

@inbounds function set_temperature!(p::Particle)
        if p.type == FLUID
                p.T = (p.ρ^(γ - 1.0)) * exp(p.s / (p.ρ * cv)) / (cv * (γ - 1.0))
        end
end

@inbounds function set_pressure!(p::Particle)
        if p.type == FLUID
                p.P = R_mass * p.ρ * p.T
        end
end

@inbounds function find_pot_temp!(p::Particle)
        if p.type == FLUID
                p.θ = p.T * (((T0 * R_gas * ρ0) / p.P)^(2))^(1 / 7) #2/7 is exactly R_gas/cp
        end
end

@inbounds function entropy_production!(p::Particle, q::Particle, r::Float64)
        if p.type == FLUID && q.type == FLUID
                ker = (q.m / q.ρ) * rDwendland2(0.5 * (p.h + q.h), r)
                x_pq = p.x - q.x
                u_pq = p.u - q.u
                p.S += -4.0 * p.m * q.m * p.ρ * ker * mu / (p.T * p.ρ * q.ρ) * dot(u_pq, x_pq)^2 / (r * r + 0.01 * p.h * q.h) * dt #viscous
        end
end

@inbounds function balance_of_smoothing!(p::Particle)
        p.Dh += -0.5 * (p.h / p.ρ) * p.Dρ
end

@inbounds function update_smoothing!(p::Particle)
        if p.type == FLUID
                p.h += dt * p.Dh
        end
        p.Dh = 0.0
end

@inbounds function update_density!(p::Particle)
        if p.type == FLUID
                p.ρ += dt * p.Dρ
        end
        p.Dρ = 0.0
end

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
        ker = (q.m / q.ρ) * rDwendland2(0.5 * (p.h + q.h), r)
        p.Dρ += p.ρ * ker * (dot(p.x - q.x, p.u - q.u))
        if p.type == FLUID && q.type == FLUID
                p.Dρ += 2 * nu / p.ρ * (p.ρ - q.ρ)
        end
end



"""
Rayleigh damping
"""

function damping_structure(z, zₜ, zᵦ, γᵣ)
        #=
                if z >= (zₜ - zᵦ)
                        return γᵣ * (sin(π / 2 * (1 - (zₜ - zᵦ) / zᵦ)))^2
                else
                        return 0
                end
        	=#
        return 0.0
end

"""
Move and accelerate
"""

function move!(p::Particle)
        if p.type == FLUID
                p.x += dt * p.u
                #p.ρ = 0.0     #we want to reset the density only for the fluid, so we call it here
        end
end

function accelerate!(p::Particle)
        if p.type == FLUID
                p.u += 0.5 * dt * (p.Du - g * VECY - damping_structure(p.x[2], zₜ, zᵦ, γᵣ) * VECY)
        end
        p.Du = VEC0
end

"""
Modifed Verlet scheme
"""

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
		E += 0.5 * p.ρ * dot(p.v,p.v) + p.ρ * g * p.x[2]
        end
        return E
end



"""
Put everything into a time loop
"""

function main()
        # save the parameters as a dictionary for future analysis
        params = @dict(
                dom_height, dom_length, dr, bc_width,
                hₘ, a,  h0,
                c, 
                γᵣ, zᵦ, zₜ,
                dt, t_end, dt_frame
        )
        sys = make_system()

        run_name = savename(params)
        run_dir = datadir("sims", run_name)
        mkpath(run_dir)

        println("Output directory: $run_dir")

        # a great help from DrWatson with metadata processing
        @tagsave(
                joinpath(run_dir, "metadata.jld2"),
                merge(params, Dict(
                        :script => basename(@__FILE__),
                        :module => "PerturbedStaticWitch"
                ))
        )


        # UNCOMMENT TO SAVE VTK FRAMES FOR PARAVIEW!
        #outpath = joinpath(RESULTS_DIR, folder_name)
        #out = new_pvd_file(outpath)
        #save_frame!(out, sys, export_vars...)
        average_velocities = DataFrame(t=Float64[], u=Float64[])
        maximum_velocities = DataFrame(t=Float64[], u=Float64[])
        ene = DataFrame(t=Float64[], E=Float64[])

        nsteps = Int(round(t_end / dt))
        frame_counter = 0

        @show T0
        @show ρ0
        @show mu
        @show c
        println("---------------------------")

        #a modified Verlet scheme
        for k = 1:nsteps
                t = k * dt
                verlet_step!(sys)
                #save data at selected 
                if (k % Int64(round(dt_frame / dt)) == 0)
                        @show t
                        println("num. of particles = ", length(sys.particles))
                        u_avg = avg_velocity(sys)
                        @show u_avg
                        push!(average_velocities, (t, u_avg))

                        u_max = max_velocity(sys)
                        @show u_max
                        push!(maximum_velocities, (t, u_max))

                        E = energy(sys)
                        @show E
                        push!(ene, (t, E))

                        E = energy(sys)
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
                                pressures=[p.P for p in sys.particles],
                                temperatures=[p.T for p in sys.particles],
                                pot_temperatures=[p.θ for p in sys.particles],
                                types=[p.type for p in sys.particles]
                        )


                        # UNCOMMENT TO SAVE VTK FRAMES FOR PARAVIEW!
                        #save_frame!(out, sys, export_vars...)

                end
        end

        CSV.write(joinpath(run_dir, "average_velocities.csv"), average_velocities)
        CSV.write(joinpath(run_dir, "maximum_velocities.csv"), maximum_velocities)
        CSV.write(joinpath(run_dir, "energies.csv"), ene)


	#UNCOMMENT IF ONE WISHES TO SAVE VTK FRAMES
        #save_pvd_file(out)
        p1 = plot(average_velocities.t, average_velocities.u;
                xlabel="t (s)",
                ylabel="avg. velocity (m/s)",
                lc=:blue,
        )
        p2 = plot(maximum_velocities.t, maximum_velocities.u;
                xlabel="t (s)",
                ylabel="max. velocity (m/s)",
                lc=:orange,
        )
        p3 = plot(
                ene.t, ene.E;
                xlabel="t (s)",
                ylabel="Total energy J",
                lc=:orange,)

        plot(p1, p2, p3; layout=(3, 1), size=(800, 900))

        # save the plots for immediate diagnostics
        savefig(plotsdir("diagnostics_$(run_name).pdf"))

        println("\n End of the road 🏵 ")
        println("\n sesults: $run_dir")
        println("\n plots: ", plotsdir("diagnostics_$(run_name).png"))
end


end # module

