using DrWatson
@quickactivate "SPH"

include(scriptsdir("hopkins_perturbed_witch.jl"))

run_dir = PerturbedStaticWitch.main()
