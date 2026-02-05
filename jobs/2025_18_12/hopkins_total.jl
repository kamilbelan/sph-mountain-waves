using DrWatson
@quickactivate "SPH"

include(scriptsdir("hopkins_total_witch.jl"))

run_dir = AdiabaticStaticWitch.main()
