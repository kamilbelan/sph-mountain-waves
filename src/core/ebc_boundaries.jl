using SmoothedParticles
using Parameters

"""
    reset_ebc_gradients!(p::Particle)

Sets grad_ρ, grad_u, grad_w to the zero vector.
"""
@inbounds function reset_ebc_gradients!(p::Particle)
	p.grad_ρ = VEC0
	p.grad_u = VEC0
	p.grad_w = VEC0
end

"""
    compute_ebc_gradients!(p::Particle, q::Particle, r::Float64)

Computes the gradient grad_ρ, grad_u, grad_w using a SPH sum
"""
@inbounds function compute_ebc_gradients!(p::Particle, q::Particle, r::Float64)
	if (p.type == FLUID) || (p.type == INFLOW_INCOMING)
		x_pq = p.x - q.x
		ker = q.m * rDwendland2(p.h, r)

		p.grad_ρ += (q.ρ - p.ρ) / p.ρ * ker * x_pq
		p.grad_u += (q.v[1] - p.v[1]) / p.ρ * ker * x_pq
		p.grad_w += (q.v[2] - p.v[2]) / p.ρ * ker * x_pq
	end
end

"""
    reset_ebc_search!(p::Particle)

Resets the particles fields `best_match_id, min_dist_sq` to perform a fresh new search in next step.
"""
@inbounds function reset_ebc_search!(p::Particle)
	if (p.type == INFLOW_GHOST) || (p.type == OUTFLOW_GHOST)
		p.best_match_id = -1
		p.min_dist_sq = Inf
	end

end

"""
search_best_extrapolator!(p::Particle, q::Particle, r::Float64, sim_params::Dict)

For the ghost particle p it finds the best (ie closest) fluid (or incoming) particle q to extrapolate some quantities and stores its id.
"""
@inbounds function search_best_extrapolator!(p::Particle, q::Particle, r::Float64, grid::ExpGrid)
	if ((p.type == INFLOW_GHOST) || (p.type == OUTFLOW_GHOST)) && ((q.type == FLUID) || (q.type == INFLOW_INCOMING))
		dr = grid.dr
		K = grid.K
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
    apply_extrapolation!(sys::ParticleSystem)

Extrapolates the density and velocity of the ghost particles in `sys` using the optimal (ie closest particles).
"""
function apply_extrapolation!(sys::ParticleSystem)
	for p in sys.particles
		if ((p.type == INFLOW_GHOST) || (p.type == OUTFLOW_GHOST)) && p.best_match_id > 0
			# retrieve the closest particle
			q = sys.particles[p.best_match_id]
			dx = p.x - q.x

			# extrapolate the ghost particle values by the closest fluid
			# we are using first order taylor here
			p.ρ = q.ρ + SmoothedParticles.dot(dx, q.grad_ρ)
			u = q.v[1] + SmoothedParticles.dot(dx, q.grad_u)
			w = q.v[2] + SmoothedParticles.dot(dx, q.grad_w)
			p.v = RealVector(u, w, 0.0)
		end
	end
	
end

"""
    manage_particle_lifecycle!(sys::ParticleSystem, grid::ExpGrid, sim_params::Dict)

Marks inflow particle inside domain as fluid particles and sets the particles exiting the domain as inflow particles.
"""
function manage_particle_lifecycle!(sys::ParticleSystem, grid::ExpGrid, dom_length::Float64)
	a_factor = (4/3)^(1/4)
	dr = grid.dr
	K = grid.K

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

"""
    set_inflow_values!(p::Particle, v_initial::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)

Assigns correct field values for particles at the inflow.
"""
function set_inflow_values!(p::Particle, v_initial::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)
	if p.type == INFLOW_INCOMING
			inflow_velocity = v_initial * VECX
			p.v = inflow_velocity
			p.ρ_bg = background_density(p.x[2], ρ0, T_bg, g, R_mass)
			p.ρ = p.ρ_bg
			p.ρ′ = 0.0
			p.Dρ = 0.0
			p.Dv = VEC0
		end
	
end
