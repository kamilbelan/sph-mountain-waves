using SmoothedParticles
using Parameters

"""
    make_system(Particle::Type, global_params::Dict, sim_params::Dict)

Creates the system's geometry and places particles of type `Particle in the positions.
"""
function make_system(Particle::Type, global_params::Dict, sim_params::Dict)
		# unpack needed parameters
		@unpack dom_height, dom_length, h_m, a = global_params
		@unpack η, dr, dt_rel, t_end, γ_r_rel = sim_params

		h0 = η * dr
		bc_width = 6*dr

		grid = Grid(dr, :hexagonal)
		domain = Rectangle(-dom_length / 2.0, 0.0, dom_length / 2.0, dom_height)
		fence = BoundaryLayer(domain, grid, bc_width)
		witch_profile(x) = (h_m * a^2) / (x^2 + a^2)
		mountain = Specification(domain, x -> (x[2] <= witch_profile(x[1])))

		sys = ParticleSystem(Particle, domain + fence, h0)

		# passing params to the Particle constructor
		generate_particles!(sys, grid, domain - mountain, x -> Particle(x, VEC0, FLUID, global_params, sim_params))
		generate_particles!(sys, grid, fence, x -> Particle(x, VEC0, WALL, global_params, sim_params))
		generate_particles!(sys, grid, mountain, x -> Particle(x, VEC0, FLUID, global_params, sim_params))

		create_cell_list!(sys)
		return sys
	end


