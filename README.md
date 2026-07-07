# SPH

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21246806.svg)](https://doi.org/10.5281/zenodo.21246806)

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> SPH

It is authored by Kamil Belan.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "SPH"
```
which auto-activate the project and enable local path handling from DrWatson.

# Running the simulations

To run the simulations, one has to do the following

- edit the global parameters in `global_params.toml`
- edit the simulation paramaters `sim_params.toml`
- run the script `run_sim.jl` e.g. like this

```
julia> using Pkg
julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
julia> Pkg.activate("path/to/this/project")
julia> Pkg.instantiate()
julia> @quickactivate "SPH"
julia> include(scriptsdir("run_sim.jl"))
```

# Results

The data produced by the scripts can be found in the `data` directory. A separate folder is created for each run.

## Acknowledgements

The SPH backend (`lib/SmoothedParticles/9tnSQ/`) is a fork of
[SmoothedParticles.jl](https://github.com/OndrejKincl/SmoothedParticles.jl)
by Ondřej Kincl and Michal Pavelka (Charles University), vendored here as a
local development version. See `lib/SmoothedParticles/9tnSQ/LICENSE` for its
license terms.

## License

This project's own code is released under the MIT License — see
[LICENSE](LICENSE). The vendored SmoothedParticles.jl fork retains its
original MIT license (Ondřej Kincl) in `lib/SmoothedParticles/9tnSQ/LICENSE`.

## Citation

If you use this code, please cite it — see [CITATION.cff](CITATION.cff)
(or use GitHub's "Cite this repository" button). Archived on Zenodo at
[10.5281/zenodo.21246806](https://doi.org/10.5281/zenodo.21246806).


