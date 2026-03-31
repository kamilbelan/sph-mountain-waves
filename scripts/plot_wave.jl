using CairoMakie, Glob, HDF5, JLD2

function plot_wave(run_dir::String)
	# create a directory plots next to the simulation output
	plots_dir = joinpath(run_dir, "plots")
	mkpath(plots_dir)

	# pick the last frame in a robust way
	last_frame = last(sort(glob("frame_*.h5", run_dir)))

	# load the data from the last frame
	pos, vel, types = h5open(last_frame, "r") do fid
		(read(fid["positions"]),
			read(fid["velocities"]),
			read(fid["types"]))
	end

	# filter to fluid only
	mask = types .== 0.0
	x  = pos[1, mask] ./ 1e3
	y  = pos[2, mask] ./ 1e3
	v_y = vel[2, mask]

	#  plot
	fig = Figure(size=(900, 450))
	ax  = Axis(fig[1,1], xlabel="Distance (km)", ylabel="Height (km)")

	sc = scatter!(ax, x, y,
		color=v_y, colormap=:RdBu, colorrange=(-0.5, 0.5),
		markersize=3)
	Colorbar(fig[1,2], sc, label="w (m s⁻¹)")

	# save as a pdf
	filepath = joinpath(plots_dir, "mountain_wave.pdf")
	save(filepath, fig)

	println("Plot saved to $filepath")
end

# allow running from command line
if !isempty(ARGS)
	plot_wave(ARGS[1])
end
