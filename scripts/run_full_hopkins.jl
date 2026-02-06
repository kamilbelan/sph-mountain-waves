using DrWatson
@quickactivate "SPH"

include(srcdir("formulations", "StaticFullHopkins.jl"))
using .StaticFullHopkins 

# define global constants (physics & geometry)
global_params = Dict(
	# physics
	:g        => 9.81,
	:R_mass   => 287.05,
	:cp       => 7 * 287.05 / 2,
	:cv       => 7 * 287.05 / 2 - 287.05,
	:γ        => 7 / 5,
	:R_gas    => 8.314,
	:T_bg     => 250.0,
	:ρ0       => 1.393,
	:N        => sqrt(0.0196),         

	# numerics
	:rho_floor => 1e-6,
	:P_floor   => 1e-10,
	:ε         => 0.01,
	:α         => 0.1, # usually α = 0.05 - 0.2
	:β         => 0.2, # usually β = 2 α

	# geometry
	:dom_height => 26e3,
	:dom_length => 400e3,
	:bc_width   => 3120.0, 
	:hₘ         => 0.0,
	:a          => 1000.0,
	:zₜ         => 26e3,
	:zᵦ         => 12e3
)

# DrWatson will generate every possible combination of these lists.
sweep_params = Dict(
	:η      => [1.5, 1.8],          
	:dr     => [26e3/25, 26e3/50, 26e3/75],            
	:t_end  => [2.0, 5.0],          
	:γᵣ     => [10 * sqrt(0.0196)]   # N = sqrt(0.0196)
)

# dict_list merges the global dict with the sweep dict and expands the lists.
job_list = dict_list(merge(global_params, sweep_params))

println("Starting StaticFullHopkins batch simulation")
println("Total simulations to run: $(length(job_list))")

for (i, params) in enumerate(job_list)
	println("\n[Job $i/$(length(job_list))] Parameters: η=$(params[:η]), t_end=$(params[:t_end])")
	output_path = StaticFullHopkins.run_sim(params)
	println("FINISHED. Saved to: $output_path")
end

println("\n BATCH COMPLETE.")
