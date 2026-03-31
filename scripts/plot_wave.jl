using CairoMakie, Glob, HDF5

function bin_to_grid(x, y, vals, xi, yi)
	grid   = fill(NaN32, length(yi), length(xi))
	counts = zeros(Int, length(yi), length(xi))
	dx = xi[2] - xi[1]
	dy = yi[2] - yi[1]
	for k in eachindex(x)
		i = round(Int, (x[k] - xi[1]) / dx) + 1
		j = round(Int, (y[k] - yi[1]) / dy) + 1
		if 1 <= i <= length(xi) && 1 <= j <= length(yi)
			grid[j, i]   = isnan(grid[j, i]) ? vals[k] : grid[j, i] + vals[k]
			counts[j, i] += 1
		end
	end
	for idx in eachindex(grid)
		counts[idx] > 0 && (grid[idx] /= counts[idx])
	end
	return grid
end

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

	# bin particles onto a regular grid
	xi = LinRange(minimum(x), maximum(x), 400)
	yi = LinRange(0, 26, 200)
	v_y_grid = bin_to_grid(x, y, v_y, xi, yi)
	θ_grid   = bin_to_grid(x, y, θ,   xi, yi)

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
