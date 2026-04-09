using SmoothedParticles
using Parameters

"""
    reset_ebc_gradients!(p::AbstractParticle)

Sets grad_δρ, grad_u, grad_w to the zero vector.
"""
function reset_ebc_gradients!(p::AbstractParticle)
	p.grad_δρ = VEC0
	p.grad_u = VEC0
	p.grad_w = VEC0
end

"""
    compute_ebc_gradients!(p::AbstractParticle, q::AbstractParticle, r::Float64)

Computes the gradient grad_δρ, grad_u, grad_w using a SPH sum
"""
function compute_ebc_gradients!(p::AbstractParticle, q::AbstractParticle, r::Float64)
	@inbounds if (p.type == FLUID || p.type == INFLOW_INCOMING) &&
	             (q.type == FLUID || q.type == INFLOW_INCOMING)
		x_pq = p.x - q.x
		ker = q.m * rDwendland2(p.h, r)

		p.grad_δρ += (q.δρ - p.δρ) / p.ρ * ker * x_pq
		p.grad_u += (q.v[1] - p.v[1]) / p.ρ * ker * x_pq
		p.grad_w += (q.v[2] - p.v[2]) / p.ρ * ker * x_pq
	end
end

"""
    apply_extrapolation!(sys::ParticleSystem)

Extrapolates the density and velocity of the ghost particles in `sys` using the optimal (ie closest particles).
"""
function apply_extrapolation!(sys::ParticleSystem, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)
	for p in sys.particles

		if ((p.type == INFLOW_GHOST) || (p.type == OUTFLOW_GHOST)) && p.best_match_id > 0 && p.best_match_id <= length(sys.particles)
			# retrieve the closest particle
			q = sys.particles[p.best_match_id]
			dx = p.x - q.x

			# extrapolate the ghost particle values by the closest fluid
			# we are using first order taylor here
			p.δρ = q.δρ + SmoothedParticles.dot(dx, q.grad_δρ)
			p.ρ = background_density(p.x[2], ρ0, T_bg, g, R_mass) + p.δρ
			u = q.v[1] + SmoothedParticles.dot(dx, q.grad_u)
			w = q.v[2] + SmoothedParticles.dot(dx, q.grad_w)
			p.v = RealVector(u, w, 0.0)
		end

		if p.type == MOUNTAIN && p.best_match_id > 0 && p.best_match_id <= length(sys.particles)
			q = sys.particles[p.best_match_id]
			dx = p.x - q.x
			p.δρ = q.δρ + SmoothedParticles.dot(dx, q.grad_δρ)
			p.ρ = background_density(p.x[2], ρ0, T_bg, g, R_mass) + p.δρ
		end
	end
	
end

"""
    set_inflow_values!(p::AbstractParticle, v_initial::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)

Assigns correct field values for particles at the inflow.
"""
function set_inflow_values!(p::AbstractParticle, v_initial::Float64, ρ0::Float64, T_bg::Float64, g::Float64, R_mass::Float64)
	if p.type == INFLOW_INCOMING
			inflow_velocity = v_initial * VECX
			p.v = inflow_velocity
			p.ρ_bg = background_density(p.x[2], ρ0, T_bg, g, R_mass)
			p.δρ = 0.0
			p.Dρ = 0.0
			p.ρ = p.ρ_bg + p.δρ
			p.Dv = VEC0
		end
end
