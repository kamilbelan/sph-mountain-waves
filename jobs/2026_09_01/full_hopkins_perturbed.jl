using DrWatson
@quickactivate "SPH"

include(scriptsdir("full_hopkins_perturbed.jl"))

run_dir = PerturbedStaticWitch.main()
