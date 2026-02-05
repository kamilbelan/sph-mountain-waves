using DrWatson
@quickactivate "SPH"

include(scriptsdir("pavelka_total_witch.jl"))

run_dir = AdiabaticStaticWitch.main()
