using JLD2
using CSV
using DataFrames
using SmoothedParticles

"""
    save_checkpoint(run_dir, sys, k, t, frame_counter, average_velocities, maximum_velocities, energies)

Save a full simulation checkpoint to `run_dir/checkpoint.jld2`.
Overwrites any existing checkpoint (we only keep the latest).
"""
function save_checkpoint(run_dir::String,
                         sys::ParticleSystem,
                         k::Int,
                         t::Float64,
                         frame_counter::Int,
                         average_velocities::DataFrame,
                         maximum_velocities::DataFrame,
                         energies::DataFrame)

    checkpoint_path = joinpath(run_dir, "checkpoint.jld2")
    tmp_path = checkpoint_path * ".tmp"

    # write to a temp file first, then atomically rename —
    # protects against SIGTERM arriving mid-write
    jldsave(tmp_path;
        particles = sys.particles,
        k = k,
        t = t,
        frame_counter = frame_counter,
        average_velocities = average_velocities,
        maximum_velocities = maximum_velocities,
        energies = energies,
    )
    mv(tmp_path, checkpoint_path, force=true)
    println("   checkpoint saved at t=$t (step $k)")
end

"""
    load_checkpoint(run_dir) -> (particles, k, t, frame_counter, avg_vel, max_vel, energies)

Load a checkpoint from `run_dir/checkpoint.jld2`.
Returns `nothing` if no checkpoint exists.
"""
function load_checkpoint(run_dir::String)
    checkpoint_path = joinpath(run_dir, "checkpoint.jld2")
    if !isfile(checkpoint_path)
        return nothing
    end

    data = load(checkpoint_path)
    return (
        particles = data["particles"],
        k = data["k"],
        t = data["t"],
        frame_counter = data["frame_counter"],
        average_velocities = data["average_velocities"],
        maximum_velocities = data["maximum_velocities"],
        energies = data["energies"],
    )
end

"""
    restore_particles!(sys, saved_particles)

Replace the particles in `sys` with those loaded from a checkpoint,
then rebuild the cell list.
"""
function restore_particles!(sys::ParticleSystem, saved_particles)
    empty!(sys.particles)
    append!(sys.particles, saved_particles)
    create_cell_list!(sys)
    println("   restored $(length(sys.particles)) particles from checkpoint")
end
