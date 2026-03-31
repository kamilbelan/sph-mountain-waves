using CairoMakie, Glob, HDF5, NaturalNeighbours 

function plot_wave(run_dir::String)
	# create a directory plots next to the simulation output
	plots_dir = joinpath(run_dir, "plots")
	mkpath(plots_dir)

	# pick the last frame in a robust way
	last_frame = last(sort(glob("frame_*.h5", run_dir)))

	# load the data from the last frame

	pos, vel, theta, types = h5open(last_frame, "r") do fid
		(read(fid["positions"]),
			read(fid["velocities"]),
			read(fid["pot_temperatures"]),
			read(fid["types"]))
	end

	# filter to fluid only
	mask = types .== 0.0
	x = pos[1, mask] ./ 1e3
	y = pos[2, mask] ./ 1e3
	v_y = vel[2, mask]
	θ = theta[mask]

	# interpolate for a line plot
	points = Matrix(hcat(x,y)') # create a 2 x N points matrix
	itp_v_y = interpolate(points, v_y)
	itp_θ = interpolate(points, θ)

	# create a regular grid
	xi = LinRange(minimum(x), maximum(x), 800)
	yi = LinRange(0, 26, 400)
	v_y_grid = [itp_v_y(xj, yj) for yj in yi, xj in xi]
	θ_grid = [itp_θ(xj, yj) for yj in yi, xj in xi]

	#  plot as in Doyle et. al
	fig = Figure(size=(900, 450))
	ax  = Axis(fig[1,1], xlabel="Distance (km)",
	    ylabel="Height (km)")

	cf = contourf!(ax, xi, yi, v_y_grid',
		levels=-0.5:0.05:0.5,
		colormap=:RdBu)
	Colorbar(fig[1,2], cf, label="w (m s⁻¹)")

	contour!(ax, xi, yi, θ_grid',
	  levels=250:10:400,
	  color=:black, linewidth=0.5)

	# save as a pdf
	filepath = joinpath(plots_dir, "mountain_wave.pdf")
	save(filepath, fig)

	println("Plot saved to $filepath")
end
