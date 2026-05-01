using DrWatson
@quickactivate "SPH"
using HDF5
using JLD2
using CairoMakie
using Printf

# Constants matching your domain setup [cite: 52]
const FLUID = 0.0

# ── Argument Parsing ──────────────────────────────────────────────────────────
# Usage: julia plot_tensile_crash.jl <run_dir> <frame_file.h5> <formulation>
if length(ARGS) < 3
	error("Usage: julia plot_tensile_crash.jl <run_dir> <frame_file.h5> <formulation>")
end
run_dir     = ARGS[1]
crash_frame = ARGS[2]
formulation = ARGS[3] # e.g., SPH or PTH

# Extract frame number from filename (e.g., "frame_0150.h5" -> "0150")
frame_num = replace(basename(crash_frame), r"frame_(\d+)\.h5" => s"\1")

# ── 1. Fast Velocity Timeline Scanner ─────────────────────────────────────────
println("Scanning frames for velocity explosion...")
all_files = sort(filter(f -> occursin(r"^frame_\d{4}\.h5$", basename(f)), readdir(run_dir; join=true)))

stride = 10 
scanned_files = all_files[1:stride:end]

if scanned_files[end] != all_files[end]
	push!(scanned_files, all_files[end])
end

times = Float64[]
l_inf = Float64[]

for fpath in scanned_files
	h5open(fpath, "r") do fid
		t_raw = read(HDF5.attributes(fid)["time"])
		t     = t_raw isa AbstractArray ? first(t_raw) : t_raw

		vel   = read(fid["velocities"])
		types = read(fid["types"])
		mask  = types .== FLUID

		vnorm = hypot.(vel[1, mask], vel[2, mask])

		push!(times, t)
		push!(l_inf, isempty(vnorm) ? 0.0 : maximum(vnorm))
	end
end

# ── 2. Load the Specific Crash Frame (Physical Space) ─────────────────────────
filepath = joinpath(run_dir, crash_frame)
t_crash, x_p, z_p = h5open(filepath, "r") do fid
	t_raw = read(HDF5.attributes(fid)["time"])
	t     = t_raw isa AbstractArray ? first(t_raw) : t_raw

	pos   = read(fid["positions"])
	types = read(fid["types"])
	mask  = types .== FLUID

	t, pos[1, mask], pos[2, mask]
end

# ── 3. Publication-Quality Figure ─────────────────────────────────────────────
fig = Figure(size=(1200, 500), fontsize=20)

# -- Left Panel: Particle Clumping --
ax1 = Axis(fig[1, 1], 
	   xlabel="distance x [km]", ylabel="height z [km]") # Cleaned labels

# Increased markersize as requested to make clumping obvious
scatter!(ax1, x_p ./ 1e3, z_p ./ 1e3, color=(:black, 0.7), markersize=5)

xlims!(ax1, -20.0, 20.0) 
ylims!(ax1, 0.0, 10.0)

# -- Right Panel: Log-Scale Velocity Explosion --
ax2 = Axis(fig[1, 2], 
	   xlabel="time t [s]", ylabel=L"|v|_{\infty}\;[\mathrm{m/s}]", 
	   xticks=0:200:1200,
	   yscale=log10)

lines!(ax2, times, l_inf, color=:red, linewidth=2.5, label=L"L^\infty \text{ norm}")
vlines!(ax2, [t_crash], color=:black, linestyle=:dash, label="snapshot")

axislegend(ax2; position=:lt, framevisible=false)

# ── Save Logic ────────────────────────────────────────────────────────────────
out_dir = plotsdir("evolutionary")
mkpath(out_dir)
out_filename = "crash_$(formulation)_$(frame_num).png"
out_path = joinpath(out_dir, out_filename)
save(out_path, fig; px_per_unit = 3)
println("Analysis complete. Plot saved to: $out_path")
