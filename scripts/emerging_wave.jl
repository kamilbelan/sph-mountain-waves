using HDF5
using JLD2
using CairoMakie
using Colors
using Printf

const FLUID = 0.0

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

length(ARGS) >= 2 || error("Usage: julia plot_emerging_wave.jl <run_dir> <frame_file.h5>")
run_dir = ARGS[1]
frame_file = ARGS[2]

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
    
    # NOTE: Change "theta" if your HDF5 key is named differently!
    theta = read(fid["pot_temperatures"]) 
    
    mask_fluid = types .== FLUID
    mask_bnd   = types .!= FLUID  
    
    rho0 = maximum(dens[mask_fluid])
    
    (read(HDF5.attributes(fid)["time"]), 
     pos[1, mask_fluid], pos[2, mask_fluid], 
     vel[2, mask_fluid], dens[mask_fluid], rho0 * dr^2,
     pos[1, mask_bnd], pos[2, mask_bnd],
     theta[mask_fluid])
end

m_f = fill(m_approx, length(x_f))

# Define plotting grid 
x_min, x_max = -50_000.0, 150_000.0  
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

# Removed the title!
ax = Axis(fig[1, 1], xlabel="Distance x [km]", ylabel="Height z [km]")

# 1. Plot the interpolated fluid field (w)
hm = heatmap!(ax, collect(x_grid)./1e3, collect(z_grid)./1e3, w_grid; 
              colormap=tol_diverging, colorrange=(-wmax, wmax), interpolate=true)

# 2. Scatter the boundary particles (Mountain)
# Increased markersize and darkened so it fuses into a solid foundation
scatter!(ax, x_b ./ 1e3, z_b ./ 1e3, color=(:black, 0.7), markersize=8)

# 3. Add Isentropes (Constant Theta)
contour!(ax, collect(x_grid)./1e3, collect(z_grid)./1e3, theta_grid; 
         levels=25, color=:black, linewidth=1.2)

# Strict axis limits
xlims!(ax, x_min/1e3, x_max/1e3)
ylims!(ax, -1.0, z_max/1e3) 

Colorbar(fig[1, 2], hm, label="w [m/s]")

mkpath("figures/evolutionary")
save("figures/evolutionary/emerging_wave.pdf", fig)
println("Saved -> figures/evolutionary/emerging_wave.pdf")
