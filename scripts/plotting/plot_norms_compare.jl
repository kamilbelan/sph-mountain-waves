# plot_norms_compare.jl — Standard vs well-balanced: Mach number, linear axes.
#
# Companion to stationary_field.jl. Together they make one argument:
#   stationary_field.jl  — WHAT goes wrong (the spurious field, Standard scheme)
#   this script          — HOW BADLY, and that WB does not do it at all.
#
# Single linear axes comparing the Standard SPH and well-balanced hydrostatic
# tests. The ordinate is the Mach number built on the reference sound speed,
#
#     Ma = ‖v‖_∞ / c,     c = sqrt(65e3 * γ / ρ0).
#
# LINEAR, not logarithmic. check_exact_zero.jl establishes that the WB run is
# at *bitwise* zero velocity — every particle, every frame. Zero has no place
# on a log axis, and clamping it to a decade floor draws a line that does not
# exist. On a linear axis Ma = 0 is a real, plottable point sitting on the
# abscissa, which is exactly the claim being made.
#
# For the same reason there is no ε_mach reference line: the WB result is not
# round-off-small, it is exactly zero, and marking machine epsilon would imply
# the error merely sits below it.
#
# ‖v‖_∞ is read from maximum_velocities.csv (written by finalize_diagnostics);
# the HDF5 frames are read only to certify the exact-zero claim and to count
# the particles and frames quoted in the annotation.
#
# Usage:
#   julia --project=@. scripts/plotting/plot_norms_compare.jl \
#         <std_run_dir> <wb_run_dir> [out_prefix]
#
# Output: plots/stationary/<out_prefix>.pdf   (vector)
#
# SIZING: the figure is emitted at 560 × 380 pt with a 13 pt base font; the
# smallest label is 11 pt. Include it at its natural size. Scaling below ≈0.9
# pushes the smallest labels under 10 pt.

using DrWatson
@quickactivate "SPH"
using CSV
using DataFrames
using HDF5
using TOML
using CairoMakie
using LaTeXStrings
using Printf

# ── constants ─────────────────────────────────────────────────────────────────
const FLUID = 0.0

const X_LO, X_HI = 0.0, 55.0
const Y_LO, Y_HI = 0.0, 1.05

# stationary_field.jl selects its snapshot as the last frame below
# BLOWUP_FRAC * c. Mirroring that rule here means the marker always points at
# the instant the companion field figure actually shows — no hardcoded time,
# and the two figures cannot drift apart.
const SNAPSHOT_FRAC = 0.5

# Paul Tol colourblind-safe palette. Wine goes to the well-balanced curve: it
# is the darker, more saturated of the two, and the WB result is the claim the
# figure exists to make, so it must carry the stronger colour.
const COL_STD = "#4477AA"   # blue — Standard SPH
const COL_WB  = "#882255"   # wine — well-balanced p-θ
const COL_REF = "#BBBBBB"   # light grey, reference line

# ── configs ───────────────────────────────────────────────────────────────────
function sound_speed(cfg_path::String)
    isfile(cfg_path) || error("global params TOML not found: $cfg_path")
    cfg = TOML.parsefile(cfg_path)
    γ  = Float64(cfg["gamma"])
    ρ0 = Float64(cfg["rho0"])
    return sqrt(65e3 * γ / ρ0), γ, ρ0
end

# ── loading ───────────────────────────────────────────────────────────────────
function read_diag_csv(run_dir::String, fname::String)
    path = joinpath(run_dir, fname)
    isfile(path) || error("missing diagnostics CSV: $path")
    df = CSV.read(path, DataFrame)
    ("t" in names(df) && "u" in names(df)) ||
        error("$path does not have the expected (t, u) columns")
    return Float64.(df.t), Float64.(df.u)
end

"""
Certify the exact-zero claim from the HDF5 frames.

Returns (n_frames, n_fluid, n_nonzero): the frame count, the per-frame FLUID
particle count, and how many particle-frames carry a velocity component that is
not bitwise +0.0. The annotation may only be drawn when n_nonzero == 0.
"""
function exact_zero_census(run_dir::String)
    frames = sort(filter(f -> occursin(r"^frame_\d{4}\.h5$", basename(f)),
                         readdir(run_dir; join=true)))
    isempty(frames) && error("No frame_*.h5 files found in: $run_dir")

    n_fluid   = -1
    n_nonzero = 0
    for fpath in frames
        h5open(fpath, "r") do fid
            vel   = read(fid["velocities"])
            types = read(fid["types"])
            mask  = types .== FLUID
            vx, vz = vel[1, mask], vel[2, mask]
            n = length(vx)

            if n_fluid < 0
                n_fluid = n
            elseif n != n_fluid
                error("FLUID count changes between frames ($n_fluid → $n); " *
                      "the per-frame count cannot be quoted as a single number.")
            end

            @inbounds for i in 1:n
                if reinterpret(UInt64, vx[i]) != 0x0000000000000000 ||
                   reinterpret(UInt64, vz[i]) != 0x0000000000000000
                    n_nonzero += 1
                end
            end
        end
    end
    return length(frames), n_fluid, n_nonzero
end

"""Mantissa-and-exponent form, e.g. 1.00347e25 → "1.0e25"."""
function fmt_ma(x::Float64)
    e = floor(Int, log10(x))
    m = x / 10.0^e
    return @sprintf("%.1fe%d", m, e)
end

# ── args ──────────────────────────────────────────────────────────────────────
length(ARGS) >= 2 || error("Usage: julia --project=@. scripts/plotting/plot_norms_compare.jl " *
                           "<std_run_dir> <wb_run_dir> [out_prefix]")

std_dir = ARGS[1]
wb_dir  = ARGS[2]
out_prefix = length(ARGS) >= 3 ? ARGS[3] : "stationary_norms_compare"

c, γ, ρ0 = sound_speed(projectdir("config", "global_params.toml"))
@printf("c = %.6f m/s   (γ = %.4f, ρ0 = %.4f)\n", c, γ, ρ0)

t_std, linf_std = read_diag_csv(std_dir, "maximum_velocities.csv")
t_wb,  linf_wb  = read_diag_csv(wb_dir,  "maximum_velocities.csv")

ma_std = linf_std ./ c
ma_wb  = linf_wb  ./ c

# ── certify the WB exact-zero claim before annotating it ──────────────────────
n_frames, n_fluid, n_nonzero = exact_zero_census(wb_dir)
@printf("WB census: %d frames × %d FLUID particles, %d non-(+0.0) components\n",
        n_frames, n_fluid, n_nonzero)
n_nonzero == 0 ||
    error("WB run is NOT at bitwise zero ($n_nonzero components differ); " *
          "the exact-zero annotation would be false — revise the figure.")
all(==(0.0), ma_wb) ||
    error("WB maximum_velocities.csv contains non-zero entries; " *
          "it disagrees with the HDF5 frames.")

# ── numbers quoted in the annotations, all taken from the data ────────────────
i_peak = argmax(ma_std)
ma_peak, t_peak = ma_std[i_peak], t_std[i_peak]
@printf("Standard peak: Ma = %.6e at t = %.4f s\n", ma_peak, t_peak)

i_exit = findfirst(>(Y_HI), ma_std)
i_exit === nothing && error("Standard curve never exceeds the ylims ceiling $Y_HI; " *
                            "the exit arrow has nothing to point at.")
t_exit = t_std[i_exit]
@printf("Standard leaves the axis (Ma > %.2f) at t = %.4f s\n", Y_HI, t_exit)

n_clipped = count(>(Y_HI), ma_std)
@printf("note: %d of %d Standard samples lie above the ceiling and are off-scale.\n",
        n_clipped, length(ma_std))

i_snap = findlast(<(SNAPSHOT_FRAC), ma_std)
if i_snap === nothing
    @warn "Standard run is never below Ma = $SNAPSHOT_FRAC; no snapshot marker drawn."
else
    @printf("field snapshot (companion figure): t = %.3f s, Ma = %.4f\n",
            t_std[i_snap], ma_std[i_snap])
end

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(560, 380), fontsize=13)

ax = Axis(fig[1, 1];
          xlabel = L"t\;[\mathrm{s}]",
          ylabel = L"\mathrm{Ma} = \|v\|_\infty / c",
          xticks = 0:10:50,
          yticks = 0:0.2:1.0,
          xticklabelsize = 12,
          yticklabelsize = 12)

xlims!(ax, X_LO, X_HI)
ylims!(ax, Y_LO, Y_HI)

# Ma = 0.5 reference
hlines!(ax, [0.5]; color=COL_REF, linestyle=:dot, linewidth=1.5)
text!(ax, 1.0, 0.515; text=L"\mathrm{Ma} = 0.5", fontsize=11, color=:gray40,
      align=(:left, :bottom))

# Instant shown by the companion field figure, drawn under the curves.
if i_snap !== nothing
    vlines!(ax, [t_std[i_snap]]; color=:gray55, linestyle=:dashdot, linewidth=1.2)
    text!(ax, t_std[i_snap] - 0.7, 0.05; text="field snapshot", fontsize=11,
          color=:gray40, rotation=π/2, align=(:left, :baseline))
end

# Plain strings, not LaTeXStrings: MathTeXEngine renders "-" as a minus sign
# even inside \text{}, giving "well − balanced p−θ". Plain strings keep the
# hyphens and match the annotation below, which names the same two schemes.
lines!(ax, t_std, ma_std; color=COL_STD, linestyle=:solid, linewidth=2.0,
       label="Standard SPH")
# Heavier and dashed: the WB curve lies exactly on the abscissa, so it has to
# read against the axis spine it coincides with.
lines!(ax, t_wb, ma_wb; color=COL_WB, linestyle=:dash, linewidth=2.8,
       label="well-balanced p-θ")

# ── annotations ───────────────────────────────────────────────────────────────
# Short upward arrow where the Standard curve leaves the top of the axis.
arrows2d!(ax, [t_exit], [0.80], [0.0], [0.20];
          color      = COL_STD,
          tipwidth   = 8,
          tiplength  = 8,
          shaftwidth = 1.6,
          align      = :tail)

text!(ax, t_exit + 1.2, 0.72;
      text  = @sprintf("→ Ma ≈ %s at t = %.1f s", fmt_ma(ma_peak), t_peak),
      fontsize = 11, color = COL_STD, align = (:left, :center))

# Sits directly above the zero line it describes, so the statement and the
# curve it refers to are read together.
text!(ax, X_HI - 1.0, 0.03;
      text  = @sprintf("well-balanced p-θ:  v ≡ +0.0 exactly\nall %d particles, all %d frames",
                       n_fluid, n_frames),
      fontsize = 11, color = COL_WB, align = (:right, :bottom))

axislegend(ax; position=:lt, labelsize=11, framevisible=false, rowgap=1)

# ── save ──────────────────────────────────────────────────────────────────────
outdir = plotsdir("stationary")
mkpath(outdir)
outfile = joinpath(outdir, "$(out_prefix).pdf")
save(outfile, fig)          # CairoMakie writes true vector PDF
println("Saved → ", abspath(outfile), "  (560 × 380 pt, vector)")
