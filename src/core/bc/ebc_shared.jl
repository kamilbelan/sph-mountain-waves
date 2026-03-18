using SmoothedParticles
using Parameters

"""
    reset_ebc_search!(p::AbstractParticle)

Resets the particles fields `best_match_id, min_dist_sq` to perform a fresh new search in next step.
"""
function reset_ebc_search!(p::AbstractParticle)
	if (p.type == INFLOW_GHOST) || (p.type == OUTFLOW_GHOST)
		p.best_match_id = -1
		p.min_dist_sq = Inf
	end

end

"""
search_best_extrapolator!(p::AbstractParticle, q::AbstractParticle, r::Float64, sim_params::Dict)

For the ghost particle p it finds the best (ie closest) fluid (or incoming) particle q to extrapolate some quantities and stores its id.
"""
function search_best_extrapolator!(p::AbstractParticle, q::AbstractParticle, dr::Float64, K::Float64)
	@inbounds if ((p.type == INFLOW_GHOST) || (p.type == OUTFLOW_GHOST)) && ((q.type == FLUID) || (q.type == INFLOW_INCOMING))
		a_factor = (4/3)^(1/4)
		a_local = dr * exp(K * max(0.0, p.x[2]) / 2.0) * a_factor

		# the grid is not regular, we need to determine ideal target point. 
		# inflow ghost (west) targets fluid +3 columns to the right.
		# outflow ghost (east) targets fluid -3 columns to the left.
		direction = p.type == INFLOW_GHOST ? 1.0 : -1.0
		x_target = p.x + RealVector(direction * 3.0 * a_local, 0.0, 0.0)

		# compute the distance of q to x_target
		vec_dist = q.x - x_target
		dist_sq = SmoothedParticles.dot(vec_dist, vec_dist)

		# update argmin
		if dist_sq < p.min_dist_sq
			p.min_dist_sq = dist_sq
			p.best_match_id = q.id
		end
	end
end

"""
    manage_particle_lifecycle!(sys::ParticleSystem, grid::ExpGrid, sim_params::Dict)

Marks inflow particle inside domain as fluid particles and sets the particles exiting the domain as inflow particles.
"""
function manage_particle_lifecycle!(sys::ParticleSystem, dr::Float64, K::Float64, dom_length::Float64)
	a_factor = (4/3)^(1/4)

	for p in sys.particles
		# retype inflow into fluid if it is in the domain
		if p.type == INFLOW_INCOMING && p.x[1] > - dom_length/2.0
			p.type = FLUID
		end

		# take the exiting particle and teleport it to inflow
		if p.type == FLUID && p.x[1] > dom_length/2.0
			# calculate the horizontal spacing at its original spawn altitude
			a_local = dr * exp(K * max(0.0, p.spawn_y) / 2.0) * a_factor

			# teleport it to the back of the incoming inflow layer
			p.x = RealVector(-dom_length/2.0 - 1.5 * a_local, p.spawn_y, 0.0)

			# set the particle type to INFLOW_INCOMING
			p.type = INFLOW_INCOMING
		end
	end
end
