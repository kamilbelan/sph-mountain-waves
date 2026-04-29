# stationary_field.jl — Figure 2: Spurious velocity field (Standard scheme)
#
# Shows the spurious velocity field that develops in the Standard (non-WB) SPH
# scheme when initialised in hydrostatic equilibrium.  The discretisation-level
# pressure-gradient error drives a non-zero circulation that should be
# identically zero in the well-balanced formulation.
#
# Usage:
#   julia --project=. scripts/stationary_field.jl <std_run_dir> [out_prefix]
#
# Output: figures/stationary/<out_prefix>.{pdf,png}

using HDF5
using JLD2
using CairoMakie
using LaTeXStrings
using Colors
using Printf

# ── tuneable parameters ───────────────────────────────────────────────────────

# Frame-selection: blowup threshold as a fraction of the reference sound speed.
const C_REF             = sqrt(65_000.0 * 1.4 / 1.393)   # ≈ 255.7 m/s
const BLOWUP_FRAC       = 0.50                            # 50 % of c_ref
const BLOWUP_THRESHOLD  = BLOWUP_FRAC * C_REF

# SPH interpolation grid resolution (for the zoomed plot window)
const N_GRID_X = 200
const N_GRID_Z = 200

# Target number of arrows across the horizontal domain
const N_ARROWS_X = 35

# ── constants ─────────────────────────────────────────────────────────────────
const FLUID = 0.0

# Sequential colormap: white → sand-light → sand → rose → wine (Paul Tol)
const CMAP_SEQ = cgrad([
    colorant"#FFFFFF",   # white
    colorant"#EEE8C0",   # sand-light
    colorant"#DDCC77",   # sand
    colorant"#CC6677",   # rose
    colorant"#882255",   # wine
])

# ── 2-D Wendland C² kernel ────────────────────────────────────────────────────
@inline function wendland2(h::Float64, r::Float64)::Float64
    q = r / h
    q >= 2.0 && return 0.0
    t = 1.0 - 0.5 * q
    return (7.0 / (4π * h * h)) * t^4 * (1.0 + 2.0 * q)
end

# ── SPH interpolation onto a regular grid ─────────────────────────────────────
function sph_interp_grid(
    xp::Vector{Float64}, zp::Vector{Float64},
    h0::Float64,
    mp::Vector{Float64}, rhop::Vector{Float64},
    fp::Vector{Float64},
    x_grid::AbstractRange, z_grid::AbstractRange,
)
    Nx  = length(x_grid)
    Nz  = length(z_grid)
    dx  = Float64(step(x_grid))
    dz  = Float64(step(z_grid))

    field_sum  = zeros(Nx, Nz)
    weight_sum = zeros(Nx, Nz)

    support = 2.0 * h0
    ix_half = ceil(Int, support / dx) + 1
    iz_half = ceil(Int, support / dz) + 1

    for j in eachindex(xp)
        wvol = mp[j] / rhop[j]   

        # nearest grid indices (1-based)
        i_cen = round(Int, (xp[j] - first(x_grid)) / dx) + 1
        k_cen = round(Int, (zp[j] - first(z_grid)) / dz) + 1

        i_lo = max(1,  i_cen - ix_half)
        i_hi = min(Nx, i_cen + ix_half)
        k_lo = max(1,  k_cen - iz_half)
        k_hi = min(Nz, k_cen + iz_half)

        for i in i_lo:i_hi, k in k_lo:k_hi
            r = hypot(x_grid[i] - xp[j], z_grid[k] - zp[j])
            w = wendland2(h0, r) * wvol
            field_sum[i, k]  += w * fp[j]
            weight_sum[i, k] += w
        end
    end

    # Shepard normalisation; uncovered cells → NaN 
    out = fill(NaN32, Nx, Nz)
    for i in 1:Nx, k in 1:Nz
        if weight_sum[i, k] > 1e-6 # <--- Fixed boundary speckle bug
            out[i, k] = Float32(field_sum[i, k] / weight_sum[i, k])
        end
    end
    return out
end

# ── frame selection ───────────────────────────────────────────────────────────
function select_frame_index(run_dir::String)
    all_files = readdir(run_dir; join=true)
    frames = sort(filter(f -> occursin(r"^frame_\d{4}\.h5$", basename(f)), all_files))
    isempty(frames) && error("No frame_*.h5 files in: $run_dir")

    println("  Scanning $(length(frames)) frames for |v|_max...")
    max_vels = Float64[]
    for fpath in frames
        h5open(fpath, "r") do fid
            vel   = read(fid["velocities"])
            types = read(fid["types"])
            mask  = types .== FLUID
            vx    = vel[1, mask]
            vz    = vel[2, mask]
            vmax  = isempty(vx) ? 0.0 : maximum(hypot.(vx, vz))
            push!(max_vels, vmax)
        end
    end

    idx = findlast(v -> v < BLOWUP_THRESHOLD, max_vels)
    if idx === nothing
        @warn "All frames exceed blowup threshold; using last frame"
        idx = length(frames)
    elseif idx == length(frames)
        println("  Simulation stayed below threshold; using last frame.")
    else
        println(@sprintf("  Blowup at frame %d; using frame %d (|v|_max = %.4g m/s)",
                         idx + 1, idx, max_vels[idx]))
    end

    return frames[idx], max_vels[idx]
end

# ── argument parsing ──────────────────────────────────────────────────────────
length(ARGS) >= 1 ||
    error("Usage: julia --project=. scripts/stationary_field.jl <run_dir> [out_prefix]")

std_dir = ARGS[1]

# Dynamic filename extraction
if length(ARGS) >= 2
    out_prefix = ARGS[2]
else
    parts = splitpath(normpath(std_dir))
    out_prefix = length(parts) >= 2 ? "$(parts[end-1])_$(parts[end])" : "stationary_field"
end

meta = load(joinpath(std_dir, "metadata.jld2"))
dr   = Float64(meta["dr"])
eta  = Float64(meta["η"])
h0   = eta * dr
println(@sprintf("Loaded metadata: dr = %.1f m, η = %.2f, h0 = %.1f m", dr, eta, h0))

println("Selecting analysis frame from: $std_dir")
frame_path, v_max_frame = select_frame_index(std_dir)
println(@sprintf("  Selected: %s   |v|_max = %.4g m/s",
                 basename(frame_path), v_max_frame))

# ── load selected frame ───────────────────────────────────────────────────────
t_snap, x_f, z_f, vx_f, vz_f, rho_f, m_approx = h5open(frame_path, "r") do fid
    t_raw  = read(HDF5.attributes(fid)["time"])
    t      = t_raw isa AbstractArray ? first(t_raw) : t_raw # <--- Fixed attribute array bug
    pos    = read(fid["positions"])
    vel    = read(fid["velocities"])
    dens   = read(fid["densities"])
    types  = read(fid["types"])
    mask   = types .== FLUID

    x   = pos[1, mask]
    z   = pos[2, mask]
    vx  = vel[1, mask]
    vz  = vel[2, mask]
    rho = dens[mask]

    # --- SAFETY CHECK FOR DELETED PARTICLES ---
    if isempty(rho)
        error("ZERO fluid particles (type == $FLUID) found in frame $(basename(frame_path)). Did the simulation explode and delete them, or is the type wrong?")
    end
    # ------------------------------------------

    rho0_approx = maximum(rho)
    m_val       = rho0_approx * dr^2

    t, x, z, vx, vz, rho, m_val
end

N_fluid = length(x_f)
println(@sprintf("  Loaded %d FLUID particles at t = %.4g s", N_fluid, t_snap))

m_f     = fill(m_approx, N_fluid)
vnorm_f = hypot.(vx_f, vz_f)
println(@sprintf("  |v|_max = %.4g m/s  (fluid)", maximum(vnorm_f)))

# ── SPH interpolation onto zoomed grid ────────────────────────────────────────
x_min_d = minimum(x_f);  x_max_d = maximum(x_f)
z_min_d = 0.0;           z_max_d = maximum(z_f)

# CROPPED PLOT WINDOW (Center 100 km)
x_min_plot = -50_000.0
x_max_plot =  50_000.0

x_grid = range(x_min_plot, x_max_plot; length=N_GRID_X)
z_grid = range(z_min_d, z_max_d; length=N_GRID_Z)

println("  Interpolating |v| onto $(N_GRID_X) × $(N_GRID_Z) grid...")
vnorm_grid = sph_interp_grid(x_f, z_f, h0, m_f, rho_f, vnorm_f, x_grid, z_grid)
println("  Done.")

# Axes in km (based on plotted region)
x_grid_km = collect(x_grid) ./ 1e3
z_grid_km = collect(z_grid) ./ 1e3
dom_len_km = (x_max_plot - x_min_plot) / 1e3
dom_hgt_km = (z_max_d - z_min_d) / 1e3

# Colour range: capped slightly higher so the gradient is smooth
valid_vals = filter(!isnan, vec(vnorm_grid))
max_val    = isempty(valid_vals) ? 0.0 : Float64(maximum(valid_vals))
vmax_grid  = min(max(1e-6, max_val), 80.0) # <--- Raised cap to 80 m/s

# ── velocity arrows on a coarse spatial grid ──────────────────────────────────
N_arrows_z = max(2, round(Int, N_ARROWS_X * dom_hgt_km / dom_len_km))
dx_a = (x_max_plot - x_min_plot) / N_ARROWS_X
dz_a = (z_max_d - z_min_d) / N_arrows_z

arrow_x  = Float64[]
arrow_z  = Float64[]
arrow_vx = Float64[]
arrow_vz = Float64[]

for ix in 0:(N_ARROWS_X - 1), iz in 0:(N_arrows_z - 1)
    xc = x_min_plot + (ix + 0.5) * dx_a
    zc = z_min_d + (iz + 0.5) * dz_a
    
    sq_dists = @. (x_f - xc)^2 + (z_f - zc)^2
    best = argmin(sq_dists)
    
    # Distance threshold to avoid drawing arrows in empty sky
    if sq_dists[best] < (1.5 * dx_a)^2
        push!(arrow_x,  x_f[best])
        push!(arrow_z,  z_f[best])
        push!(arrow_vx, vx_f[best])
        push!(arrow_vz, vz_f[best])
    end
end

v_arrow_max    = max(maximum(hypot.(arrow_vx, arrow_vz)), 1e-10)
arrow_scale_km = 0.03 * dom_len_km / v_arrow_max
println(@sprintf("  Arrow scale: %.3g km/(m/s)  (v_max for arrows = %.4g m/s)",
                 arrow_scale_km, v_arrow_max))

# ── figure ────────────────────────────────────────────────────────────────────
t_label = @sprintf("%.4g", t_snap)

# Adjusted figure size for the ~4:1 cropped aspect ratio
fig = Figure(size=(1000, 350), fontsize=20)

ax = Axis(fig[1, 1];
    xlabel         = "distance x[km]",
    ylabel         = "height z[km]",
    xticklabelsize = 16,
    yticklabelsize = 16,
    aspect         = DataAspect(), 
)

# Heatmap with color clipping
hm = heatmap!(ax, x_grid_km, z_grid_km, vnorm_grid;
    colormap   = CMAP_SEQ,
    colorrange = (0.0, vmax_grid),
)

Colorbar(fig[1, 2], hm;
    label         = "velocity |v|[m/s]",
    labelsize     = 18,
    ticklabelsize = 16,
    width         = 15,
)

# Velocity arrows overlay with increased size
arrows2d!(ax,
    arrow_x  ./ 1e3,
    arrow_z  ./ 1e3,
    arrow_vx .* arrow_scale_km,
    arrow_vz .* arrow_scale_km;
    color      = RGBf(0.15, 0.15, 0.15),
    tipwidth   = 8,
    tiplength  = 8,
    shaftwidth = 0.8,  
    align      = :center,
)

# ── save ──────────────────────────────────────────────────────────────────────
outdir = normpath(joinpath(@__DIR__, "..", "figures", "stationary"))
mkpath(outdir)

save(joinpath(outdir, "$(out_prefix).pdf"), fig)
save(joinpath(outdir, "$(out_prefix).png"), fig; px_per_unit=3)
println("Saved → $(abspath(outdir))/$(out_prefix).{pdf,png}")
