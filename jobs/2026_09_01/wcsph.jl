using DrWatson
@quickactivate "SPH"

include(scriptsdir("wcsph.jl"))

run_dir = PerturbedStaticWitch.main()
