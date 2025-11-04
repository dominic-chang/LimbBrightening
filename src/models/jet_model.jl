using Krang
using Roots

function magnetic_stream_function(r, θ, p) # 2410.00954v1 Eq (11)
	return (r^p)*(1-abs(cos(θ)))
end

function lapse(metric, r, θ)
	Δt = Krang.Δ(metric, r)
	Ξt = Krang.Ξ(metric, r, θ; Δ = Δt)
	Σt = Krang.Σ(metric, r, θ)

	return sqrt(Σt*Δt/Ξt)#Eq 3 2410.00954v1
end

function normal_observer_BL_d(metric, r, θ)#Eq 3 2410.00954v1
	return [-lapse(metric, r, θ), 0, 0, 0]
end

function dψ_dθ(r, θ, p)
	signθ = θ > π/2 ? -1 : 1
	return signθ*r^p*sin(θ)
end
function dψ_dr(r, θ, p)
	return p*r^(p-1)*(1-abs(cos(θ)))
end
function Br(metric, r, θ, p) #Eq 6 2410.00954v1

	α = lapse(metric, r, θ)
	g = Krang.met_det(metric, r, θ)
	return α*dψ_dθ(r, θ, p)/(√(-g))
end
function Bθ(metric, r, θ, p)#Eq 6 2410.00954v1
	α = lapse(metric, r, θ)
	g = Krang.met_det(metric, r, θ)
	return -α*dψ_dr(r, θ, p)/(√(-g))
end

function magnetic_field_line_angular_rotation(metric, θp, p)#Eq 36 & 37 2410.00954v1
	a = metric.spin
	if p == 0
		return a/8
	else
		return a/(4*(1+sec(θp/2)^2))
	end
end

function magnetic_field_BL_u(metric, r, θ, p)
	α = lapse(metric, r, θ)
	g = Krang.met_det(metric, r, θ)
	cosθ = cos(θ)
	signθ = sign(cosθ)
	ψ = magnetic_stream_function(r, θ, p)
	rp = Krang.horizon(metric)
	θp = acos(1-ψ/rp^p)
	ΩF = magnetic_field_line_angular_rotation(metric, θp, p)
	I = (p > 0 ? -4π*ψ*ΩF*signθ : -2π*ψ*ΩF*(1+abs(cos(θp))))
	Bbar = I/(2π)
	Δt = Krang.Δ(metric, r)
	sinθ2 = 1-cosθ^2

	return [0.0, dψ_dθ(r, θ, p)/(√(-g)), -dψ_dr(r, θ, p)/(√(-g)), Bbar/(Δt*sinθ2)]
end

function Nco(metric, r, ψ, p)#Eq A19 2410.00954v1
	θ = acos(1-ψ/r^p)
	rp = Krang.horizon(metric)
	θp = acos(1-ψ/rp^p)
	ΩF = magnetic_field_line_angular_rotation(metric, θp, p)
	met_dd = metric_dd(metric, r, θ)
	return -(met_dd[1, 1] + 2*met_dd[1, 4]*ΩF + met_dd[4, 4]*ΩF^2)
end

function dNco_dr(metric, r, θp, p)
	a = metric.spin
	rp = Krang.horizon(metric)
	ψ = magnetic_stream_function(rp, θp, p)

	# Common terms
	r_p = r^p
	sqrt_1_minus_a2 = sqrt(1 - a^2)
	one_plus_sqrt = 1 + sqrt_1_minus_a2
	ψ_term = (1 + sqrt_1_minus_a2)^(-p) * ψ
	sec_arccos = abs(sec(1/2 * acos(1 - ψ_term)))
	ΩF = magnetic_field_line_angular_rotation(metric, θp, p)

	denom = (r^(2 + 2p) + a^2 * (r^p - ψ)^2)^2

	# Prefactor r^(-1 - 2p)/denom
	pref = r^(-1 - 2p) / denom

	# Numerator terms (direct translation)
	t1 = -2a * r^(3 + 4p) * ψ * (2*(1 + p)*r^p - (1 + 2p)*ψ) * ΩF

	t2 = -2a^3 * r^(1 + 2p) * (r^p - ψ) * ψ * (2*(-1 + p)*r^(2p) + 3*r^p*ψ - ψ^2) * ΩF

	t3 = a^6 * p * (r^p - ψ)^5 * ψ * ΩF^2

	t4 =
		-a^4 * r * ψ * (-r^p + ψ) *
		(
			(-2 + 3p) * r^(1 + 4p) +
			4 * (-1 + 2p) * r^(3p) * ψ -
			(-7 + 8p) * r^(1 + 3p) * ψ -
			4 * (-2 + 3p) * r^(2p) * ψ^2 +
			(-9 + 8p) * r^(1 + 2p) * ψ^2 +
			(-5 + 8p) * r^p * ψ^3 -
			(-5 + 4p) * r^(1 + p) * ψ^3 +
			(1 - 2p) * ψ^4 +
			(-1 + p) * r * ψ^4
		) * ΩF^2

	t5 = r^(3 + 4p) * (
		r^(2p) + (-2 + p) * r^(3 + p) * ψ * ΩF^2 - (-1 + p) * r^3 * ψ^2 * ΩF^2
	)

	t6 =
		a^2 * r^(1 + 2p) *
		(
			-r^(4p) + 2 * (1 + p) * r^(3p) * ψ - (1 + 2p) * r^(2p) * ψ^2 +
			(-4 + 3p) * r^(3 + 3p) * ψ * ΩF^2 +
			4 * (1 + 2p) * r^(2 + 2p) * ψ^2 * ΩF^2 - (-10 + 7p) * r^(3 + 2p) * ψ^2 * ΩF^2 -
			4 * (1 + 3p) * r^(2 + p) * ψ^3 * ΩF^2 + 2 * (-4 + 3p) * r^(3 + p) * ψ^3 * ΩF^2 +
			(1 + 4p) * r^2 * ψ^4 * ΩF^2 - 2 * (-1 + p) * r^3 * ψ^4 * ΩF^2
		)

	numer = t1 + t2 + t3 + t4 + t5 + t6

	return pref * numer

end

function d2Nco_dr2(metric, r, θp, p)
	a = metric.spin
	rp = Krang.horizon(metric)
	ψ = magnetic_stream_function(rp, θp, p)

	denom = (r^(2 + 2p) + a^2 * (r^p - ψ)^2)^3

	pref = - r^(-2*(1 + p)) / denom

	t1 = 2 * r^(5 + 8p)

	t2 = -4 * a^5 * p * r^(1 + 4p) * (r^p - ψ)^2 * ψ * (((-1 + p) * r^p) + ψ + 2p * ψ) * ΩF

	t3 = -4 * a * (1 + p) * r^(5 + 6p) * ψ * ((2 + p) * r^p - (1 + 2p) * ψ) * ΩF

	t4 = -4 * a^3 * r^(3 + 4p) * ψ * (
			 2 * (-3 + p + p^2) * r^(3p) -
			 (-15 + 2 * (-4 + p) * p) * r^(2p) * ψ -
			 3 * (1 + p) * (4 + p) * r^p * ψ^2 +
			 (1 + p) * (3 + 2p) * ψ^3
		 ) * ΩF

	t5 = (-2 + p) * (-1 + p) * r^(8 + 7p) * ψ * ΩF^2

	t6 = -(-1 + p) * (-1 + 2p) * r^(8 + 6p) * ψ^2 * ΩF^2

	t7 = a^8 * p * (r^p - ψ)^6 * ψ * ((1 + p) * r^p - (1 + 2p) * ψ) * ΩF^2

	t8 =
		a^6 * r * (r^p - ψ)^2 * ψ *
		(
			2 * (1 + 2p^2) * r^(1 + 5p) -
			r^(4p) * (9r + 2p * (4 - 3r + p * (-8 + 9r))) * ψ +
			r^(3p) * (4 * (5 - 7p) * p + (16 + p * (-21 + 29p)) * r) * ψ^2 -
			r^(2p) * (4 * (5 - 8p) * p + (14 + p * (-27 + 22p)) * r) * ψ^3 +
			r^p * (p * (10 + 9p * (-2 + r) - 15r) + 6r) * ψ^4 -
			(-1 + 2p) * (p * (-2 + r) - r) * ψ^5
		) * ΩF^2

	t9 =
		a^2 * r^(3 + 4p) *
		(
			-6 * r^(4p) + 2 * (1 + p) * (6 + p) * r^(3p) * ψ -
			2 * (3 + p) * (1 + 2p) * r^(2p) * ψ^2 +
			2 * (3 + 2 * (-2 + p) * p) * r^(3 + 3p) * ψ * ΩF^2 +
			8 * (1 + p) * (1 + 2p) * r^(2 + 2p) * ψ^2 * ΩF^2 -
			(15 + 2p * (-13 + 7p)) * r^(3 + 2p) * ψ^2 * ΩF^2 -
			4 * (2 + 9p * (1 + p)) * r^(2 + p) * ψ^3 * ΩF^2 +
			3 * (-1 + p) * (-4 + 5p) * r^(3 + p) * ψ^3 * ΩF^2 +
			2 * (1 + 6p + 8p^2) * r^2 * ψ^4 * ΩF^2 -
			3 * (-1 + p) * (-1 + 2p) * r^3 * ψ^4 * ΩF^2
		)

	t10 =
		a^4 * r^(1 + 2p) * ψ *
		(
			3 * r^(2 + p) * (p * (10 + 9p * (-2 + r) - 15r) + 6 * (2 + r)) * ψ^4 * ΩF^2 -
			3 * (-1 + p) * (-2 + 2p * (-2 + r) - r) * r^2 * ψ^5 * ΩF^2 +
			2 * r^(5p) * ((-1 + p) * p + 3 * (1 + (-1 + p) * p) * r^3 * ΩF^2) +
			3 * r^(3p) * ψ^2 * (-2p * (1 + p) + r^2 * (-27p * r + 8 * (3 + 2r) + p^2 * (-32 + 19r)) * ΩF^2) -
			r^(4p) * ψ * (-6p + r^2 * (-4p * (4 + 9r) + 3 * (8 + 9r) + p^2 * (-32 + 30r)) * ΩF^2) -
			r^(2p) * ψ^3 * (-2p * (1 + 2p) + r^2 * (78 + 42r + p * (40 - 87r + 2p * (-52 + 27r))) * ΩF^2)
		)

	numer = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10

	return pref * numer
end

function stagnation_surface(metric, θp, p)
	ro = horizon(metric)
	Roots.find_zero((x->dNco_dr(metric, x, θp, p), x->d2Nco_dr2(metric, x, θp, p)), ro)
end

function Eco(metric, ro, θo, θp, p)#Eq A8 2410.00954v1
	ΩF = magnetic_field_line_angular_rotation(metric, θp, p)
	rp = Krang.horizon(metric)
	ψ = magnetic_stream_function(rp, θp, p)

	normalization = 1/sqrt(Nco(metric, ro, ψ, p))
	uo_u = normalization*[1, 0, 0, ΩF]
	uo_d = Krang.metric_dd(metric, ro, θo) * uo_u
	return -uo_d[1]-ΩF*uo_d[4]
end

function electric_field_BL_u(metric, r, θ, B_BL_u) # Eq D136 2307.06372

	Δt = Krang.Δ(metric, r)
	Σt = Krang.Σ(metric, r, θ)
	Ξt = Krang.Ξ(metric, r, θ; Δ = Δt)

	ωt = Krang.ω(metric, r, θ; Ξ = Ξt)
	_, Br, Bθ, _ = B_BL_u
	signθ = θ > π/2 ? -1 : 1
	return [0, signθ*(ΩF-ωt)*Ξt*sin(θ)*Bθ/Σt, -signθ*(ΩF-ωt)*Ξt*sin(θ)*Br/(Σt*Δt), 0]

end

function velocity_ff_BL_u(metric, r, θp, p, B_BL_u, E_BL_u)
	rp = Krang.horizon(metric)
	ψ = magnetic_stream_function(rp, θp, p)
	θ = acos(1-ψ/r^p)
	met_dd = Krang.metric_dd(metric, r, θ)
	α = lapse(metric, r, θ)
	ΩF = magnetic_field_line_angular_rotation(metric, θp, p)

	B_BL_d = met_dd * B_BL_u
	E_BL_d = met_dd * E_BL_u

	_, Br_BL_d, Bθ_BL_d, Bϕ_BL_d = B_BL_d
	_, Er_BL_d, Eθ_BL_d, Eϕ_BL_d = E_BL_d

	B2 = B_BL_u' * B_BL_d
	E2 = E_BL_u' * E_BL_d

	Bhat_BL_u = B_BL_u ./ sqrt(B2)
	g = √(-Krang.met_det(metric, r, θ))

	n_BL_u = Krang.metric_uu(metric, r, θ) * normal_observer_BL_d(metric, r, θ)

	v_drift_BL_u = (α/(B2*g)) .* [0, (Eθ_BL_d*Bϕ_BL_d - Bθ_BL_d*Eϕ_BL_d), (Eϕ_BL_d*Br_BL_d - Bϕ_BL_d*Er_BL_d), (Er_BL_d*Bθ_BL_d - Br_BL_d*Eθ_BL_d)]
	γo = 1/sqrt(abs(1-E2/B2)) # TODO: Fix the need for abs here

	met_dd = Krang.metric_dd(metric, r, θ)
	h = (met_dd[1, 4] + met_dd[4, 4]*ΩF)*Bhat_BL_u[4]
	f = γo*(α-(met_dd[1, 4] + met_dd[4, 4]*ΩF)*v_drift_BL_u[4])
	ro = stagnation_surface(metric, θp, p)
	θo = acos(1-ψ/ro^p)

	pm = sign(r-ro)
	eco = Eco(metric, ro, θo, θp, p)
	ξ = (h*f + pm*eco*sqrt(abs(eco^2-f^2+h^2)))/(eco^2+h^2)# TODO: fix abs
	γ = γo/√(1-ξ^2)
	return γ*(n_BL_u + v_drift_BL_u + ξ/γo*Bhat_BL_u)
end
