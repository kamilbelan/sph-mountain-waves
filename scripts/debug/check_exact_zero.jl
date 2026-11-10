# check_exact_zero.jl — is the WB stationary run *exactly* at rest?
#
# Loads every HDF5 frame of a run and reports, for FLUID particles only:
#   • max |v|                        (full 17-digit precision)
#   • how many particles have v bitwise == 0.0, out of the FLUID count
#   • max |v| / c,  c = sqrt(65e3 * γ / ρ0)
#
# and then classifies the run as one of
#   (1) max|v| is exactly 0.0 at every frame
#   (2) max|v| ~ 1e-16 and flat in time
#   (3) max|v| ~ 1e-16 and growing
#
# Usage:
#   julia --project=@. scripts/debug/check_exact_zero.jl <run_dir> [global_params.toml]

using DrWatson
@quickactivate "SPH"
using HDF5
using TOML
using Printf

const FLUID = 0.0

# ── configs ───────────────────────────────────────────────────────────────────
# γ and ρ0 live in the *global* params, not in the per-run metadata (which only
# records sim_params). Prefer an explicit path, else the project's config file.
function load_gamma_rho0(cfg_path::String)
    isfile(cfg_path) || error("global params TOML not found: $cfg_path")
    cfg = TOML.parsefile(cfg_path)
    haskey(cfg, "gamma") || error("no `gamma` in $cfg_path")
    haskey(cfg, "rho0")  || error("no `rho0` in $cfg_path")
    return Float64(cfg["gamma"]), Float64(cfg["rho0"]), cfg_path
end

# ── frame reading ─────────────────────────────────────────────────────────────
"""
Return (t, max|v|, n_bitwise_zero, n_neg_zero, n_fluid) for one frame.

`n_bitwise_zero` counts particles whose vx and vz are both the bit pattern of
+0.0. `n_neg_zero` counts those that are numerically zero but carry a sign bit
(-0.0) in at least one component — numerically identical, bitwise not.
"""
function frame_stats(fpath::String)
    h5open(fpath, "r") do fid
        t_raw = read(HDF5.attributes(fid)["time"])
        t     = t_raw isa AbstractArray ? first(t_raw) : t_raw

        vel   = read(fid["velocities"])   # 3 × N
        types = read(fid["types"])        # N
        mask  = types .== FLUID

        vx = vel[1, mask]
        vz = vel[2, mask]
        n_fluid = length(vx)
        n_fluid == 0 && error("frame has no FLUID particles: $fpath")

        vmax = 0.0
        n_bit = 0
        n_neg = 0
        @inbounds for i in 1:n_fluid
            x, z = vx[i], vz[i]
            m = hypot(x, z)
            m > vmax && (vmax = m)

            bx = reinterpret(UInt64, x)
            bz = reinterpret(UInt64, z)
            if bx == 0x0000000000000000 && bz == 0x0000000000000000
                n_bit += 1
            elseif x == 0.0 && z == 0.0
                n_neg += 1        # ±0.0 mix: numerically zero, bitwise not +0.0
            end
        end
        return t, vmax, n_bit, n_neg, n_fluid
    end
end

# ── args ──────────────────────────────────────────────────────────────────────
length(ARGS) >= 1 ||
    error("Usage: julia --project=@. scripts/debug/check_exact_zero.jl <run_dir> [global_params.toml]")

run_dir  = ARGS[1]
cfg_path = length(ARGS) >= 2 ? ARGS[2] : projectdir("config", "global_params.toml")

isdir(run_dir) || error("run directory does not exist: $run_dir")

frames = sort(filter(f -> occursin(r"^frame_\d{4}\.h5$", basename(f)),
                     readdir(run_dir; join=true)))
isempty(frames) && error("No frame_*.h5 files found in: $run_dir")

γ, ρ0, used_cfg = load_gamma_rho0(cfg_path)
c = sqrt(65e3 * γ / ρ0)

println("run       : ", run_dir)
println("frames    : ", length(frames))
println("configs   : ", used_cfg)
@printf("γ         : %.17e\n", γ)
@printf("ρ0        : %.17e\n", ρ0)
@printf("c         : %.17e  m/s   [c = sqrt(65e3 * γ / ρ0)]\n", c)
@printf("eps(F64)  : %.17e\n", eps(Float64))
println()

# ── per-frame report ──────────────────────────────────────────────────────────
println(rpad("frame", 12), rpad("t [s]", 24), rpad("max|v| [m/s]", 26),
        rpad("v==0 bitwise / N_fluid", 26), "max|v|/c")
println("-"^110)

times   = Float64[]
vmaxes  = Float64[]
neg_tot = 0

for f in frames
    t, vmax, n_bit, n_neg, n_fluid = frame_stats(f)
    push!(times, t); push!(vmaxes, vmax)
    global neg_tot += n_neg

    @printf("%-12s %.17e   %.17e   %9d / %-9d   %.17e\n",
            basename(f)[1:end-3], t, vmax, n_bit, n_fluid, vmax / c)
end

println()
if neg_tot > 0
    @printf("note: %d particle-frames were numerically zero but not bitwise +0.0 (signed zero).\n\n", neg_tot)
end

# ── classification ────────────────────────────────────────────────────────────
n_zero_frames = count(==(0.0), vmaxes)
vmax_overall  = maximum(vmaxes)

@printf("frames with max|v| exactly 0.0 : %d / %d\n", n_zero_frames, length(vmaxes))
@printf("overall max|v|                 : %.17e m/s\n", vmax_overall)
@printf("overall max|v| / c             : %.17e\n", vmax_overall / c)
println()

if n_zero_frames == length(vmaxes)
    verdict = 1
else
    # Growth test on the non-zero tail: compare the last quarter of the record
    # against the first non-zero quarter. A well-balanced scheme sitting on
    # round-off noise stays flat; a leaking one drifts upward.
    n  = length(vmaxes)
    q  = max(1, n ÷ 4)
    early = maximum(vmaxes[1:q])
    late  = maximum(vmaxes[end-q+1:end])
    grew  = early > 0.0 ? (late / early) : (late > 0.0 ? Inf : 1.0)
    @printf("first-quarter max : %.17e\n", early)
    @printf("last-quarter  max : %.17e\n", late)
    @printf("growth ratio      : %.6g\n\n", grew)
    verdict = grew > 10.0 ? 3 : 2
end

println("="^70)
if verdict == 1
    println("VERDICT: case (1) — max|v| is exactly 0.0 at EVERY frame.")
    println("The WB run is at machine-exact rest; zeros must be clamped to the")
    println("axis floor for display on a log scale, and the caption must say so.")
elseif verdict == 2
    println("VERDICT: case (2) — max|v| is at round-off level and FLAT in time.")
    println("Equilibrium is preserved to machine precision; no secular drift.")
else
    println("VERDICT: case (3) — max|v| is at round-off level but GROWING.")
    println("Round-off is being amplified; the balance is not exactly preserved.")
end
println("="^70)
