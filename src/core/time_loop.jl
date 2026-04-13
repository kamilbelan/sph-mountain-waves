using SmoothedParticles
using DataFrames
using Parameters

"""
    time_loop(global_params, run_dir, step_function!, t_end, dt, dt_frame, sys,
              average_velocities, maximum_velocities, energies;
              k_start=0, frame_counter=0)

Advances the simulation from step `k_start` to `nsteps`.
When resuming from a checkpoint, pass `k_start` and `frame_counter`
so the loop continues where it left off.
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
		   energies::DataFrame;
		   k_start::Int=0,
		   frame_counter::Int=0,
		   )

	@unpack g = global_params
	nsteps = Int(round(t_end / dt))

	println("\n" * "="^60)

	save_interval = max(1, Int(round(dt_frame / dt)))

	# if simulation is already complete (e.g. chained job with nothing left to do), exit early
	if k_start >= nsteps
		println("Simulation already complete (step $k_start / $nsteps). Nothing to do.")
		return
	end

	# only save the initial frame on a fresh run
	if k_start == 0
		t = 0.0
		write_frame!(run_dir, sys, frame_counter, t)
		frame_counter += 1
	else
		println("Resuming from step $k_start / $nsteps (t=$(k_start * dt))")
	end

	# introduce a try-catch block to obtain the xdmf file even when the sim crashes
	try
		for k = (k_start + 1):nsteps
			t = k * dt
			step_function!(sys)

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
				flush(stdout) # force the output

				# write the data at the given time to a h5 file
				write_frame!(run_dir, sys, frame_counter, t)
				frame_counter += 1

				# save checkpoint (overwrites previous)
				save_checkpoint(run_dir, sys, k, t, frame_counter,
						average_velocities, maximum_velocities, energies)
			end # if
		end # for
	catch e
		@warn "Time loop terminated TOO early" exception=(e, catch_backtrace())
	end
end # function
