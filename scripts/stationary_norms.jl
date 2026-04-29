# stationary_norms.jl — Figure 1: Velocity-norm time series
#
# Compares the L∞ and L² (RMS) norms of |v| over time for the Standard and
# well-balanced (WB) SPH hydrostatic-equilibrium tests.  Both norms are
# computed directly from the HDF5 frames so the result is independent of the
# diagnostic CSVs.  Only FLUID particles (type == 0.0) are included.
#
# Usage:
#   julia --project=. scripts/stationary_norms.jl <std_run_dir> <wb_run_dir> [out_prefix]
#
# Output: figures/stationary/<out_prefix>.{pdf,png}

using HDF5
using CairoMakie
using LaTeXStrings
using Printf

# ── constants ─────────────────────────────────────────────────────────────────
const FLUID = 0.0
const FLOOR = 1e-16   # log-scale floor for near-zero WB norms (see note above)

# Paul Tol colorblind-safe palette
const COL_STD = "#882255"   # wine   → Standard scheme
const COL_WB  = "#332288"   # indigo → WB scheme

# ── data loading ──────────────────────────────────────────────────────────────

function load_velocity_norms(run_dir::String)
    all_files = readdir(run_dir; join=true)
    frames = sort(filter(f -> occursin(r"^frame_\d{4}\.h5$", basename(f)), all_files))
    isempty(frames) && error("No frame_*.h5 files found in: $run_dir")
    println("  Found $(length(frames)) frames in $(basename(run_dir))")

    times = Float64[]
    l_inf = Float64[]
    l2    = Float64[]

    for fpath in frames
        h5open(fpath, "r") do fid
            t_raw = read(HDF5.attributes(fid)["time"])
            t     = t_raw isa AbstractArray ? first(t_raw) : t_raw # <--- Array bug fix
            vel   = read(fid["velocities"])    # 3 × N matrix
            types = read(fid["types"])         # N vector
            mask  = types .== FLUID
            vx    = vel[1, mask]
            vz    = vel[2, mask]
            vnorm = hypot.(vx, vz)
            N     = length(vnorm)

            linf_val = N > 0 ? maximum(vnorm)                       : 0.0
            l2_val   = N > 0 ? sqrt(sum(v^2 for v in vnorm) / N)    : 0.0

            push!(times, t)
            push!(l_inf, linf_val)
            push!(l2,    l2_val)
        end
    end

    return times, l_inf, l2
end

# ── argument parsing ──────────────────────────────────────────────────────────
length(ARGS) >= 2 ||
    error("Usage: julia --project=. scripts/stationary_norms.jl <std_dir> <wb_dir> [out_prefix]")

std_dir = ARGS[1]
wb_dir  = ARGS[2]

# Dynamic filename extraction
if length(ARGS) >= 3
    out_prefix = ARGS[3]
else
    variant = basename(normpath(std_dir))
    out_prefix = "stationary_norms_$(variant)"
end

println("Loading Standard run: $std_dir")
t_std, linf_std, l2_std = load_velocity_norms(std_dir)
println(@sprintf("  L∞ in [%.4g, %.4g] m/s", minimum(linf_std), maximum(linf_std)))

println("Loading WB run: $wb_dir")
t_wb, linf_wb, l2_wb = load_velocity_norms(wb_dir)
println(@sprintf("  L∞ in [%.4g, %.4g] m/s", minimum(linf_wb), maximum(linf_wb)))

# Apply floor before log (clips exact-zero WB values to FLOOR)
clamp_fl(v) = max.(v, FLOOR)
linf_std_p = clamp_fl(linf_std)
l2_std_p   = clamp_fl(l2_std)
linf_wb_p  = clamp_fl(linf_wb)
l2_wb_p    = clamp_fl(l2_wb)

println("Norms ready; building figure...")

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(900, 500), fontsize=20)

ax = Axis(fig[1, 1];
    xlabel         = "time t[s]",
    ylabel         = "velocity |v|[m/s]",
    yscale         = log10,
    xticklabelsize = 16,
    yticklabelsize = 16,
)

lines!(ax, t_std, linf_std_p; color=COL_STD, linestyle=:solid, linewidth=2.0,
       label=L"\mathrm{Standard},\;L^\infty")
lines!(ax, t_std, l2_std_p;   color=COL_STD, linestyle=:dash,  linewidth=2.0,
       label=L"\mathrm{Standard},\;L^2\;\mathrm{(RMS)}")
lines!(ax, t_wb,  linf_wb_p;  color=COL_WB,  linestyle=:solid, linewidth=2.0,
       label=L"\mathrm{WB},\;L^\infty")
lines!(ax, t_wb,  l2_wb_p;    color=COL_WB,  linestyle=:dash,  linewidth=2.0,
       label=L"\mathrm{WB},\;L^2\;\mathrm{(RMS)}")

# Legend: upper-left avoids the growing Standard lines.
axislegend(ax; position=:rb, labelsize=16, framevisible=true)

# ── save ──────────────────────────────────────────────────────────────────────
outdir = normpath(joinpath(@__DIR__, "..", "figures", "stationary"))
mkpath(outdir)

save(joinpath(outdir, "$(out_prefix).pdf"), fig)
save(joinpath(outdir, "$(out_prefix).png"), fig; px_per_unit=3)
println("Saved → $(abspath(outdir))/$(out_prefix).{pdf,png}")
