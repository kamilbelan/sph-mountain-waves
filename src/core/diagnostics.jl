using Plots
using CSV
using DataFrames
using SmoothedParticles

# ==============
# Diagnostics helper functions
# ==============

"""
    avg_velocity(sys::ParticleSystem)

Calculate and return the average velocity magnitude of all particles in the system.
"""

function avg_velocity(sys::ParticleSystem)::Float64
	v = 0.0
	count = 0.0
	for p in sys.particles
		if p.type == FLUID
			v += SmoothedParticles.norm(p.v)
			count += 1.0
		end
	end
	v = v / count
	return v
end

"""
    max_velocity(sys::ParticleSystem)

Calculate and return the maximal velocity magnitude of all particles in the system.
"""

function max_velocity(sys::ParticleSystem)::Float64
	v = maximum(SmoothedParticles.norm(p.v) for p in sys.particles if p.type == FLUID)
	return v
end

"""
    energy(sys::ParticleSystem, g::Float64)

Compute the total energy (kinetic + potential) of the fluid system.
"""
function energy(sys::ParticleSystem, g::Float64)::Float64
	E = 0.0
	for p in sys.particles
		if p.type == FLUID
			E += 0.5 * p.ρ * SmoothedParticles.dot(p.v, p.v) + p.ρ * g * p.x[2]
		end
	end
	return E
end

# ==============
# Diagnostics helpers
# ==============


"""
    initialize_diagnostics_arrays()

Create and return three empty DataFrames to track average velocity, 
maximum velocity, and total energy over time.
"""

function initialize_diagnostics_arrays()
	# initialize data accumulation for diagnostics
	average_velocities = DataFrame(t=Float64[], u=Float64[])
	maximum_velocities = DataFrame(t=Float64[], u=Float64[])
	energies = DataFrame(t=Float64[], E=Float64[])
	return average_velocities, maximum_velocities, energies
end

"""
    finalize_diagnostics(run_dir, average_velocities, maximum_velocities, energies)

Write the diagnostic DataFrames to CSV files in  `run_dir and generate 
a  plot of the results.
"""

function finalize_diagnostics(run_dir::String, average_velocities::DataFrame,  maximum_velocities::DataFrame, energies::DataFrame)
	# write CSVs once at the end
	CSV.write(joinpath(run_dir, "average_velocities.csv"), average_velocities)
	CSV.write(joinpath(run_dir, "maximum_velocities.csv"), maximum_velocities)
	CSV.write(joinpath(run_dir, "energies.csv"), energies)

	println("Generating diagnostic plots...")

	# fallback to headless plotting if not in an interactive REPL (eg for cluster use)
	if !isinteractive()
		ENV["GKSwstype"] = "nul"
	end

	p1 = plot(average_velocities.t, average_velocities.u; xlabel="t (s)", ylabel="avg. velocity (m/s)", lc=:blue)
	p2 = plot(maximum_velocities.t, maximum_velocities.u; xlabel="t (s)", ylabel="max. velocity (m/s)", lc=:purple)
	p3 = plot(energies.t, energies.E; xlabel="t (s)", ylabel="Total energy J", lc=:orange)

	plt = plot(p1, p2, p3; layout=(3, 1), size=(800, 900))

	# display the plot if local, otherwise skip displaying on HPC
	if isinteractive()
		display(plt)
	end

	# save the plot into the run directory
	savefig(plt, joinpath(run_dir, "diagnostics.pdf"))
end


