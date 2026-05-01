using DrWatson
@quickactivate "SPH"
using HDF5
using JLD2
using CairoMakie
using Colors
using Printf

const FLUID = 0.0
const MOUNTAIN = 1.0 

# ── 2-D Wendland C² kernel ────────────────────────────────────────────────────
@inline function wendland2(h::Float64, r::Float64)::Float64
    q = r / h
    q >= 2.0 && return 0.0
    t = 1.0 - 0.5 * q
    return (7.0 / (4π * h * h)) * t^4 * (1.0 + 2.0 * q)
end

# ── SPH interpolation onto a regular grid ─────────────────────────────────────
function sph_interp_grid(xp, zp, h0, mp, rhop, fp, x_grid, z_grid)
    Nx, Nz = length(x_grid), length(z_grid)
    dx, dz = Float64(step(x_grid)), Float64(step(z_grid))
    field_sum, weight_sum = zeros(Nx, Nz), zeros(Nx, Nz)
    support = 2.0 * h0
    ix_half, iz_half = ceil(Int, support / dx) + 1, ceil(Int, support / dz) + 1

    for j in eachindex(xp)
        wvol = mp[j] / rhop[j]   
        i_cen = round(Int, (xp[j] - first(x_grid)) / dx) + 1
        k_cen = round(Int, (zp[j] - first(z_grid)) / dz) + 1

        for i in max(1, i_cen - ix_half):min(Nx, i_cen + ix_half)
            for k in max(1, k_cen - iz_half):min(Nz, k_cen + iz_half)
                r = hypot(x_grid[i] - xp[j], z_grid[k] - zp[j])
                w = wendland2(h0, r) * wvol
                field_sum[i, k]  += w * fp[j]
                weight_sum[i, k] += w
            end
        end
    end
    out = fill(NaN32, Nx, Nz)
    for i in 1:Nx, k in 1:Nz
        if weight_sum[i, k] > 1e-6
            out[i, k] = Float32(field_sum[i, k] / weight_sum[i, k])
        end
    end
    return out
end

# ── Updated Argument Parsing ──────────────────────────────────────────────────
# Usage: julia plot_emerging_wave.jl <run_dir> <frame_file.h5> <formulation>
length(ARGS) >= 3 || error("Usage: julia plot_emerging_wave.jl <run_dir> <frame_file.h5> <formulation>")
run_dir     = ARGS[1]
frame_file  = ARGS[2]
formulation = ARGS[3] # SPH or PTH

# Extract frame number from filename (e.g., "frame_0100.h5" -> "0100")
frame_num = replace(basename(frame_file), r"frame_(\d+)\.h5" => s"\1")

# Load metadata
meta = load(joinpath(run_dir, "metadata.jld2"))
dr, eta = Float64(meta["dr"]), Float64(meta["η"])
h0 = eta * dr

# Load frame
filepath = joinpath(run_dir, frame_file)
t_snap, x_f, z_f, vz_f, rho_f, m_approx, x_b, z_b, theta_f = h5open(filepath, "r") do fid
    pos   = read(fid["positions"])
    vel   = read(fid["velocities"])
    dens  = read(fid["densities"])
    types = read(fid["types"])
    theta = read(fid["pot_temperatures"]) 
    
    mask_fluid = types .== FLUID
    mask_bnd   = types .== MOUNTAIN  
    
    rho0 = maximum(dens[mask_fluid])
    
    (read(HDF5.attributes(fid)["time"]), 
     pos[1, mask_fluid], pos[2, mask_fluid], 
     vel[2, mask_fluid], dens[mask_fluid], rho0 * dr^2,
     pos[1, mask_bnd], pos[2, mask_bnd],
     theta[mask_fluid])
end

m_f = fill(m_approx, length(x_f))

# Define plotting grid 
x_min, x_max = -40_000.0, 50_000.0  
z_min, z_max = 0.0, 14_000.0
x_grid = range(x_min, x_max; length=300)
z_grid = range(z_min, z_max; length=150)

println("Interpolating vertical velocity (w)...")
w_grid = sph_interp_grid(x_f, z_f, h0, m_f, rho_f, vz_f, x_grid, z_grid)

println("Interpolating potential temperature (theta)...")
theta_grid = sph_interp_grid(x_f, z_f, h0, m_f, rho_f, theta_f, x_grid, z_grid)

# ── Aesthetic Settings ────────────────────────
tol_diverging = cgrad([
    colorant"#882255", colorant"#AA4466", colorant"#CC6677", colorant"#E8A5AE",
    colorant"#FFFFFF", colorant"#B8DCF0", colorant"#88CCEE", colorant"#5577BB", colorant"#332288"
])
wmax = 0.8 

# ── Figure ────────────────────────────────────
fig = Figure(size=(1000, 450), fontsize=22)
ax = Axis(fig[1, 1], xlabel="distance x [km]", ylabel="height z [km]")

# 1. Plot the interpolated fluid field (w)
hm = heatmap!(ax, collect(x_grid)./1e3, collect(z_grid)./1e3, w_grid; 
              colormap=tol_diverging, colorrange=(-wmax, wmax), interpolate=true)

# 2. Scatter the mountain bedrock
scatter!(ax, x_b ./ 1e3, z_b ./ 1e3, color=(:darkgray, 1.0), markersize=4.5)

# 3. Add Isentropes (Constant Theta)
contour!(ax, collect(x_grid)./1e3, collect(z_grid)./1e3, theta_grid; 
         levels=25, color=(:black, 0.5), linewidth=0.8)

xlims!(ax, x_min/1e3, x_max/1e3)
ylims!(ax, -0.4, z_max/1e3) 

Colorbar(fig[1, 2], hm, label="w [m/s]")

# ── Updated Save Logic ────────────────────────────────────────────────────────
out_dir = plotsdir("evolutionary")
mkpath(out_dir)

# Filename pattern: wave_FORMULATION_FRAME.pdf
out_filename = "wave_$(formulation)_$(frame_num).pdf"
save(joinpath(out_dir, out_filename), fig)

println("Saved -> $(joinpath(out_dir, out_filename))")
