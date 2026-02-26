using SmoothedParticles
using DataFrames
using Parameters

"""
    time_loop(global_params::Dict, 
		   run_dir::String,
		   step_function!::Function,
		   t_end::Float64, 
		   dt::Float64,
		   dt_frame::Float64,
		   sys::ParticleSystem,
		   average_velocities::DataFrame,
		   maximum_velocities::DataFrame,
		   energies::DataFrame,
		   )

Advances the simulation by one time-step.
"""
function time_loop(global_params::Dict, 
		   run_dir::String,
		   step_function!::Function,
		   t_end::Float64, 
		   dt::Float64,
		   dt_frame::Float64,
		   sys::ParticleSystem,
		   average_velocities::DataFrame,
		   maximum_velocities::DataFrame,
		   energies::DataFrame,
		   )

	@unpack g = global_params
	nsteps = Int(round(t_end / dt))
	frame_counter = 0

	println("\n" * "="^60)

	# save the initial frame
	t = 0.0
	write_frame!(run_dir, sys, frame_counter, t)

	for k = 1:nsteps
		t = k * dt
		step_function!(sys)

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
			push!(energies, (t, E))

			# write the data at the given time to a h5 file
			write_frame!(run_dir, sys, frame_counter, t)
			frame_counter += 1
		end # if
	end # for
end # function
