using SmoothedParticles
using Parameters

import SmoothedParticles: Grid2, covering, is_inside, boundarybox, RealVector, Shape


"""
ExpGrid(dr, K)

Creates an isotropically spaced grid whose spacing decays exponentially.
"""
mutable struct ExpGrid <: Grid2
	dr::Float64
	K::Float64
end

function covering(grid::ExpGrid, s::Shape)::Vector{RealVector}
	rect = boundarybox(s)
	dr_base = grid.dr
	K = grid.K

	# hexagonal factors 
	a_factor = (4/3)^(1/4)
	b_factor = (3/4)^(1/4)

	xs = RealVector[]

	# upper half
	y = 0.0
	layer = 0
	while y <= rect.x2_max
		s_local = dr_base * exp(K * y / 2.0)
		a_local = s_local * a_factor
		b_local = s_local * b_factor

		if y >= rect.x2_min
			i_min = Int(floor(rect.x1_min / a_local))
			i_max = Int(ceil(rect.x1_max / a_local))
			offset = iseven(layer) ? 0.0 : 0.5 * a_local

			for i in i_min:i_max
				pos = RealVector((i * a_local) + offset, y, 0.0)
				if is_inside(pos, s)
					push!(xs, pos)
				end
			end
		end

		y += b_local
		layer += 1
	end
	# lower half
	y = 0.0
	layer = 0
	while true
		# step down into the boundary
		layer -= 1

		s_local = dr_base * exp(K * y / 2.0)
		b_local = s_local * b_factor
		a_local = s_local * a_factor

		y -= b_local 

		if y < rect.x2_min
			break
		end

		if y <= rect.x2_max
			i_min = Int(floor(rect.x1_min / a_local))
			i_max = Int(ceil(rect.x1_max / a_local))
			offset = iseven(layer) ? 0.0 : 0.5 * a_local

			for i in i_min:i_max
				pos = RealVector((i * a_local) + offset, y, 0.0)
				if is_inside(pos, s)
					push!(xs, pos)
				end
			end
		end
	end   
	return xs
end
"""
make_system(Particle::Type, global_params::Dict, sim_params::Dict)

Creates the system's geometry and places particles of type `Particle in the positions.
"""
function make_system(Particle::Type, global_params::Dict, sim_params::Dict)
	# unpack needed parameters

	@unpack dom_height, dom_length, a = global_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel, h_m = sim_params
	@unpack g, R_mass, T_bg = global_params

	K = g / (R_mass * T_bg)
	h0 = η * dr
	bc_width = 12 * dr 

	# place the particles into ExpGrid
	grid = ExpGrid(dr, K)

	domain = Rectangle(-dom_length / 2.0, 0.0, dom_length / 2.0, dom_height)
	fence = BoundaryLayer(domain, grid, bc_width)
	witch_profile(x) = (h_m * a^2) / (x^2 + a^2)
	mountain = Specification(domain, x -> (x[2] <= witch_profile(x[1])))

	# AD HOC 3.0!
	sys = ParticleSystem(Particle, domain + fence, 3.0*h0)

	# passing params to the Particle constructor
	generate_particles!(sys, grid, domain - mountain, x -> Particle(x, VEC0, FLUID, global_params, sim_params))
	generate_particles!(sys, grid, fence, x -> Particle(x, VEC0, WALL, global_params, sim_params))
	generate_particles!(sys, grid, mountain, x -> Particle(x, VEC0, FLUID, global_params, sim_params))

	create_cell_list!(sys)
	return sys
end


