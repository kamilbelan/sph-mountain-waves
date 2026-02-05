using DrWatson
@quickactivate "SPH"

include(scriptsdir("wcsph_perturbed_witch.jl"))

run_dir = PerturbedStaticWitch.main()
