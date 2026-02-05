using DrWatson
@quickactivate "SPH"

include(scriptsdir("full_hopkins_perturbed_witch.jl"))

run_dir = PerturbedStaticWitch.main()
