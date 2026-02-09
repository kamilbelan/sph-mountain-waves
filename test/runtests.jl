using DrWatson, Test
@quickactivate "SPH"

# Here you include files using `srcdir`
#include(srcdir("full_hopkins_perturbed.jl"))
#include(scriptsdir("full_hopkins_perturbed_witch.jl"))
#include(scriptsdir("hopkins_perturbed_witch.jl"))
#include(scriptsdir("hopkins_total_witch.jl"))
#include(scriptsdir("pavelka_total_witch.jl"))
#include(scriptsdir("wcsph_perturbed_witch.jl"))

# Run test suite
#println("Starting tests")
#ti = time()

#@testset "SPH tests" begin
#    @test 1 == 1
#end

#ti = time() - ti
#println("\nTest took total time of:")
#println(round(ti/60, digits = 3), " minutes")
