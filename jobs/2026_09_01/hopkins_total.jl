using DrWatson
@quickactivate "SPH"

include(scriptsdir("hopkins_total.jl"))

run_dir = AdiabaticStaticWitch.main()
