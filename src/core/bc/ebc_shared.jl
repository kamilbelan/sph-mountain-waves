using SmoothedParticles
using Parameters

"""
    reset_ebc_search!(p::AbstractParticle)

Resets the ghost particles fields `best_match_id, min_dist_sq` to perform a fresh new search in next step.
"""
function reset_ghost_ebc_search!(p::AbstractParticle)
	if (p.type == INFLOW_GHOST) || (p.type == OUTFLOW_GHOST)
		p.best_match_id = -1
		p.min_dist_sq = Inf
	end

end

"""
    reset_mountain_ebc_search!(p::AbstractParticle)

Resets the mountain particles fields `best_match_id, min_dist_sq` to perform a fresh new search in next step.
"""
function reset_mountain_ebc_search!(p::AbstractParticle)
	if p.type == MOUNTAIN
		p.best_match_id = -1
		p.min_dist_sq = Inf
	end

end


"""
    search_ghost_extrapolator!(p::AbstractParticle, q::AbstractParticle, r::Float64, sim_params::Dict)

For the ghost particle p it finds the best (ie closest) fluid (or incoming) particle q to extrapolate some quantities and stores its id.
"""
function search_ghost_extrapolator!(p::AbstractParticle, q::AbstractParticle, dr::Float64, K::Float64)
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
    mountain_normal(x, h_m, a)

Computes the outward unit normal of the Witch of Agnesi surface h(x) = h_m*a²/(x²+a²) at horizontal position x.
Points upward into the fluid.
"""
function mountain_normal(x::Float64, h_m::Float64, a::Float64)::RealVector
	dh_dx = -2.0 * h_m * a^2 * x / (x^2 + a^2)^2
	n = RealVector(-dh_dx, 1.0, 0.0)
	return n / SmoothedParticles.norm(n)
end

"""
    search_mountain_extrapolator!(p, q, r)

For each mountain wall particle p, finds the nearest fluid (or incoming) particle q.
Uses pure nearest-neighbor: unlike inflow/outflow ghosts, mountain wall particles
sit directly below well-supported fluid particles, so no offset target is needed.
"""
function search_mountain_extrapolator!(p::AbstractParticle, q::AbstractParticle, r::Float64)
	@inbounds if p.type == MOUNTAIN && (q.type == FLUID || q.type == INFLOW_INCOMING)
		dist_sq = r * r
		if dist_sq < p.min_dist_sq
			p.min_dist_sq = dist_sq
			p.best_match_id = q.id
		end
	end
end

"""
    apply_mountain_freeslip!(sys, h_m, a)

Enforces the free-slip velocity BC at the mountain wall particles by first-order extrapolation
of the nearest fluid particle's velocity, then reflecting the normal component.
"""
function apply_mountain_freeslip!(sys::ParticleSystem, h_m::Float64, a::Float64)
	for p in sys.particles
		if p.type == MOUNTAIN && p.best_match_id > 0
			# retrieve the closest particle
			q = sys.particles[p.best_match_id]
			dx = p.x - q.x

			# first-order (Taylor) extrapolation of velocity
			u = q.v[1] + SmoothedParticles.dot(dx, q.grad_u)
			w = q.v[2] + SmoothedParticles.dot(dx, q.grad_w)
			v_e = RealVector(u, w, 0.0)

			# reflect normal component to enforce v·n = 0 at the surface (zero average above/below the layer)
			n = mountain_normal(p.x[1], h_m, a)
			v_normal = SmoothedParticles.dot(v_e, n)
			p.v = v_e - 2.0 * v_normal * n
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


