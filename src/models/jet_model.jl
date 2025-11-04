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

function Bϕ(metric, r, θ, p)#Eq 7 2410.00954v1

	cosθ = cos(θ)
	signθ = sign(cosθ)
	ψ = magnetic_stream_function(r, θ, p)
	rp = Krang.horizon(metric)
	θp = acos(1-ψ/rp^p)
	ΩF = magnetic_field_line_angular_rotation(metric, θp, p)
	I = (p > 0 ? -4π*ψ*ΩF*signθ : -2π*ψ*ΩF*(1+abs(cosθ)))#Eq 31 2410.00954v1
	Bbar = I/(2π)#Eq 29 2410.00954v1
	Δt = Krang.Δ(metric, r)

	return Bbar/(Δt*(1-cosθ^2))
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

function Er(metric, ψ, r, θ, p, Bθ)
	rp = Krang.horizon(metric)
	Σt = Krang.Σ(metric, r, θ)
	Ξt = Ξ(metric, r, θ; Δ = Δt)
	ωt = ω(metric, r, θ; Ξ = Ξt)
	θp = acos(1-ψ/rp^p)
	ΩF = magnetic_field_line_angular_rotation(metric, θp, p)
	At = Krang.A(metric, r, θ; Δ = Δ(metric, r))
	return (ΩF-ωt)*At*sin(θ)*Bθ/Σt
end

function Eθ(metric, ψ, r, θ, p, Br)
	rp = Krang.horizon(metric)
	Σt = Krang.Σ(metric, r, θ)
	Ξt = Ξ(metric, r, θ; Δ = Δt)
	ωt = ω(metric, r, θ; Ξ = Ξt)
	θp = acos(1-ψ/rp^p)
	ΩF = magnetic_field_line_angular_rotation(metric, θp, p)
	Δt = Δ(metric, r)
	At = Krang.A(metric, r, θ; Δ = Δt)
	return -(ΩF-ωt)*At*sin(θ)*Br/(Σt*Δt)
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

	# Helper functions
	function sec(x)
		return 1/cos(x)
	end

	# Common terms
	r_p = r^p
	sqrt_1_minus_a2 = sqrt(1 - a^2)
	one_plus_sqrt = 1 + sqrt_1_minus_a2
	ψ_term = (1 + sqrt_1_minus_a2)^(-p) * ψ
	sec_arccos = sec(1/2 * acos(1 - ψ_term))

	# Base denominator terms
	base_denom = (r^(2 + 2p) + a^2 * (r_p - ψ)^2)^2
	sec_denom = (1 + sec_arccos^2)^2

	# Main numerator terms
	term1 = 16r^(3 + 6p)

	term2 = a^8 * p * (r_p - ψ)^5 * ψ

	term3 = -a^2 * r^(1 + 4p) * (
		16r^(2p) - 32(1 + p)r_p * ψ +
		16(1 + p)r^(2 + p) * ψ - (-2 + p)r^(5 + p) * ψ +
		16(1 + 2p)ψ^2 - 8(1 + 2p)r^2 * ψ^2 + (-1 + p)r^5 * ψ^2
	)

	term4 =
		a^4 * r^(1 + 2p) * ψ *
		(
			-16(-1 + p)r^(3p) + (-4 + 3p)r^(3 + 3p) +
			8(-5 + 2p)r^(2p) * ψ + 4(1 + 2p)r^(2 + 2p) * ψ -
			(-10 + 7p)r^(3 + 2p) * ψ + 32r_p * ψ^2 -
			4(1 + 3p)r^(2 + p) * ψ^2 + 2(-4 + 3p)r^(3 + p) * ψ^2 -
			8ψ^3 + (1 + 4p)r^2 * ψ^3 - 2(-1 + p)r^3 * ψ^3
		)

	term5 = -a^6 * r * ψ * (-r_p + ψ) * (
		(-2 + 3p)r^(1 + 4p) + 4(-1 + 2p)r^(3p) * ψ -
		(-7 + 8p)r^(1 + 3p) * ψ - 4(-2 + 3p)r^(2p) * ψ^2 +
		(-9 + 8p)r^(1 + 2p) * ψ^2 + (-5 + 8p)r_p * ψ^3 -
		(-5 + 4p)r^(1 + p) * ψ^3 + (1 - 2p)ψ^4 + (-1 + p)r * ψ^4
	)

	term6 = 8r^(1 + 2p) * (
				4r^(2 + 4p) - a^4 * (r_p - ψ) * ψ * (
					2(-1 + p)r^(2p) + 3r_p * ψ - ψ^2
				) - a^2 * r^(2p) * (
					4r^(2p) - 8(1 + p)r_p * ψ + 2(1 + p)r^(2 + p) * ψ +
					4(1 + 2p)ψ^2 - (1 + 2p)r^2 * ψ^2
				)
			) * sec_arccos^2

	term7 = 16r^(1 + 4p) * (
				r^(2 + 2p) - a^2 * (r_p - ψ) * (r_p - (1 + 2p)ψ)
			) * sec_arccos^4

	# Final assembly
	numerator = r^(-1 - 2p) * (term1 + term2 + term3 + term4 + term5 + term6 + term7)
	denominator = 8 * base_denom * sec_denom

	return numerator / denominator
end

function d2Nco_dr2(metric, r, θp, p)
	a = metric.spin
	rp = Krang.horizon(metric)
	ψ = magnetic_stream_function(rp, θp, p)

	# Helper functions
	function sec(x)
		return 1/cos(x)
	end

	# Common terms
	r_p = r^p
	sqrt_1_minus_a2 = sqrt(1 - a^2)
	one_plus_sqrt = 1 + sqrt_1_minus_a2
	ψ_term = (1 + sqrt_1_minus_a2)^(-p) * ψ
	sec_arccos = sec(1/2 * acos(1 - ψ_term))

	# Base denominator terms
	base_denom = (r^(2 + 2p) + a^2 * (r_p - ψ)^2)^4
	sec_denom = (1 + sec_arccos^2)^2

	# Main numerator terms
	term1 = 96r^(7 + 10p)

	term2 = a^2 * r^(5 + 8p) * (
		-576r^(2p) +
		32(36 + 47p + 12p^2 + p^3)r_p * ψ -
		16(6 + 11p + 6p^2 + p^3)r^(2 + p) * ψ +
		p * (2 - 3p + p^2)r^(5 + p) * ψ -
		32(18 + 47p + 24p^2 + 4p^3)ψ^2 +
		16(3 + 11p + 12p^2 + 4p^3)r^2 * ψ^2 -
		2p * (1 - 3p + 2p^2)r^5 * ψ^2
	)

	term3 =
		a^4 * r^(3 + 6p) *
		(
			96r^(4p) +
			64(-6 - 13p + 6p^2 + p^3)r^(3p) * ψ -
			48(-12 - p + 4p^2 + p^3)r^(2 + 3p) * ψ +
			p * (10 - 9p + 5p^2)r^(5 + 3p) * ψ +
			192(3 + 13p + 4p^2)r^(2p) * ψ^2 +
			16(-90 - 103p - 24p^2 + 4p^3)r^(2 + 2p) * ψ^2 +
			8(3 + 11p + 12p^2 + 4p^3)r^(4 + 2p) * ψ^2 -
			2p * (13 - 21p + 14p^2)r^(5 + 2p) * ψ^2 -
			192(2 + 13p + 14p^2 + 3p^3)r_p * ψ^3 +
			96(12 + 25p + 16p^2 + 3p^3)r^(2 + p) * ψ^3 -
			12(2 + 11p + 18p^2 + 9p^3)r^(4 + p) * ψ^3 +
			12p * (2 - 5p + 3p^2)r^(5 + p) * ψ^3 +
			32(3 + 26p + 48p^2 + 16p^3)ψ^4 -
			32(9 + 25p + 24p^2 + 8p^3)r^2 * ψ^4 +
			2(3 + 22p + 48p^2 + 32p^3)r^4 * ψ^4 -
			8p * (1 - 3p + 2p^2)r^5 * ψ^4
		)

	term4 =
		a^10 * p * r * (r_p - ψ)^3 * ψ *
		(
			(10 + 9p + 5p^2)r^(1 + 6p) +
			8(-1 + 4p^2)r^(5p) * ψ - (44 + 39p + 37p^2)r^(1 + 5p) * ψ -
			4(-7 + 13p^2)r^(4p) * ψ^2 +
			6(13 + 8p + 15p^2)r^(1 + 4p) * ψ^2 +
			20(-2 + 5p^2)r^(3p) * ψ^3 -
			6(12 - p + 17p^2)r^(1 + 3p) * ψ^3 -
			30(-1 + 3p^2)r^(2p) * ψ^4 +
			(38 - 51p + 61p^2)r^(1 + 2p) * ψ^4 +
			6(-2 + 7p^2)r_p * ψ^5 -
			3(4 - 11p + 7p^2)r^(1 + p) * ψ^5 +
			2(1 - 4p^2)ψ^6 + 2(1 - 3p + 2p^2)r * ψ^6
		)

	term5 = a^12 * p * (1 + p) * (r_p - ψ)^8 * ψ * (
		(2 + p)r_p - 2ψ * (1 + 2p)
	)

	term6 =
		2a^6 * r^(1 + 4p) * ψ *
		(
			-6r^2 * (r_p - ψ)^2 * (2r_p - ψ) * (
				4r^(2p) - 8r_p * ψ + 6r^(2 + p) * ψ + 4ψ^2 - 3r^2 * ψ^2
			) +
			p * (r_p - ψ) *
			(
				-16r^(4p) + 120r^(2 + 4p) + 10r^(5 + 4p) +
				64r^(3p) * ψ - 64r^(2 + 3p) * ψ - 12r^(4 + 3p) * ψ -
				32r^(5 + 3p) * ψ - 96r^(2p) * ψ^2 - 320r^(2 + 2p) * ψ^2 +
				218r^(4 + 2p) * ψ^2 + 40r^(5 + 2p) * ψ^2 + 64r_p * ψ^3 +
				352r^(2 + p) * ψ^3 - 212r^(4 + p) * ψ^3 - 24r^(5 + p) * ψ^3 -
				16ψ^4 - 88r^2 * ψ^4 + 53r^4 * ψ^4 + 6r^5 * ψ^4
			) -
			3p^2 * r^2 *
			(
				16r^(5p) + r^(3 + 5p) + 96r^(4p) * ψ - 32r^(2 + 4p) * ψ -
				10r^(3 + 4p) * ψ - 192r^(3p) * ψ^2 + 56r^(2 + 3p) * ψ^2 +
				32r^(3 + 3p) * ψ^2 - 44r^(3 + 2p) * ψ^3 + 112r_p * ψ^4 -
				28r^(2 + p) * ψ^4 + 27r^(3 + p) * ψ^4 - 32ψ^5 + 8r^2 * ψ^5 -
				6r^3 * ψ^5
			) +
			p^3 * (
				16r^(5p) - 24r^(2 + 5p) + 5r^(5 + 5p) + 64r^(4p) * ψ -
				32r^(2 + 4p) * ψ + 48r^(4 + 4p) * ψ - 36r^(5 + 4p) * ψ -
				224r^(3p) * ψ^2 + 256r^(2 + 3p) * ψ^2 - 182r^(4 + 3p) * ψ^2 +
				84r^(5 + 3p) * ψ^2 + 128r^(2p) * ψ^3 - 192r^(2 + 2p) * ψ^3 +
				208r^(4 + 2p) * ψ^3 - 92r^(5 + 2p) * ψ^3 + 80r_p * ψ^4 -
				40r^(2 + p) * ψ^4 - 97r^(4 + p) * ψ^4 + 51r^(5 + p) * ψ^4 -
				64ψ^5 + 32r^2 * ψ^5 + 20r^4 * ψ^5 - 12r^5 * ψ^5
			)
		)

	term7 =
		2a^8 * r^(1 + 2p) * (r_p - ψ) * ψ *
		(
			3r^2 * (r_p - ψ)^3 * ψ * (-2r_p + ψ)^2 +
			2p * (r_p - ψ)^2 * (
				4r^(4p) + 5r^(3 + 4p) - 8r^(3p) * ψ - 30r^(2 + 3p) * ψ -
				14r^(3 + 3p) * ψ + 4r^(2p) * ψ^2 + 23r^(2 + 2p) * ψ^2 +
				15r^(3 + 2p) * ψ^2 - 8r^(2 + p) * ψ^3 - 8r^(3 + p) * ψ^3 +
				2r^2 * ψ^4 + 2r^3 * ψ^4
			) +
			3p^2 * r^2 * (r_p - ψ) * (
				r^(1 + 5p) + 16r^(4p) * ψ + 12r^(3p) * ψ^2 -
				13r^(1 + 3p) * ψ^2 - 8r^(2p) * ψ^3 + 26r^(1 + 2p) * ψ^3 -
				18r^(1 + p) * ψ^4 + 4r * ψ^5
			) +
			p^3 * (
				-8r^(6p) + 5r^(3 + 6p) - 40r^(5p) * ψ + 48r^(2 + 5p) * ψ -
				39r^(3 + 5p) * ψ + 72r^(4p) * ψ^2 - 154r^(2 + 4p) * ψ^2 +
				105r^(3 + 4p) * ψ^2 + 8r^(3p) * ψ^3 + 198r^(2 + 3p) * ψ^3 -
				139r^(3 + 3p) * ψ^3 - 32r^(2p) * ψ^4 - 172r^(2 + 2p) * ψ^4 +
				102r^(3 + 2p) * ψ^4 + 84r^(2 + p) * ψ^5 - 42r^(3 + p) * ψ^5 -
				16r^2 * ψ^6 + 8r^3 * ψ^6
			)
		)

	term8 =
		16r^(1 + 4p) *
		(
			12r^(6 + 6p) +
			a^8 * p * ψ * (-r_p + ψ)^3 * (
				(-1 + p^2)r^(2p) + (2 + 7p^2)r_p * ψ + (-1 + 4p^2)ψ^2
			) -
			a^2 * r^(4 + 4p) * (
				72r^(2p) - 4(36 + 47p + 12p^2 + p^3)r_p * ψ +
				(6 + 11p + 6p^2 + p^3)r^(2 + p) * ψ +
				4(18 + 47p + 24p^2 + 4p^3)ψ^2 -
				(3 + 11p + 12p^2 + 4p^3)r^2 * ψ^2
			) +
			a^4 * r^(2 + 2p) *
			(
				12r^(4p) + 8(-6 - 13p + 6p^2 + p^3)r^(3p) * ψ -
				3(-12 - p + 4p^2 + p^3)r^(2 + 3p) * ψ +
				24(3 + 13p + 4p^2)r^(2p) * ψ^2 +
				(-90 - 103p - 24p^2 + 4p^3)r^(2 + 2p) * ψ^2 -
				24(2 + 13p + 14p^2 + 3p^3)r_p * ψ^3 +
				6(12 + 25p + 16p^2 + 3p^3)r^(2 + p) * ψ^3 +
				4(3 + 26p + 48p^2 + 16p^3)ψ^4 -
				2(9 + 25p + 24p^2 + 8p^3)r^2 * ψ^4
			) -
			a^6 * (r_p - ψ) * ψ *
			(
				3r^2 * (r_p - ψ)^3 * (2r_p - ψ) -
				p * (r_p - ψ)^2 * (
					-4r^(2p) + 15r^(2 + 2p) + 8r_p * ψ + 22r^(2 + p) * ψ -
					4ψ^2 - 11r^2 * ψ^2
				) +
				6p^2 * r^2 * (
					r^(4p) + 7r^(3p) * ψ - 5r^(2p) * ψ^2 - 5r_p * ψ^3 + 2ψ^4
				) +
				p^3 * (
					-4r^(4p) + 3r^(2 + 4p) - 20r^(3p) * ψ + 7r^(2 + 3p) * ψ +
					36r^(2p) * ψ^2 - 25r^(2 + 2p) * ψ^2 + 4r_p * ψ^3 -
					r^(2 + p) * ψ^3 - 16ψ^4 + 4r^2 * ψ^4
				)
			)
		) * sec_arccos^2

	term9 =
		32r^(1 + 4p) *
		(
			3r^(6 + 6p) +
			a^6 * p * (r_p - ψ)^3 * ψ * (
				(-1 + p^2)r^(2p) + (2 + 7p^2)r_p * ψ + (-1 + 4p^2)ψ^2
			) -
			a^2 * r^(4 + 4p) * (
				18r^(2p) - (36 + 47p + 12p^2 + p^3)r_p * ψ +
				(18 + 47p + 24p^2 + 4p^3)ψ^2
			) +
			a^4 * r^(2 + 2p) * (r_p - ψ) * (
				3r^(3p) + (-9 - 26p + 12p^2 + 2p^3)r^(2p) * ψ +
				(9 + 52p + 36p^2 + 2p^3)r_p * ψ^2 -
				(3 + 26p + 48p^2 + 16p^3)ψ^3
			)
		) * sec_arccos^4

	# Final assembly
	numerator = r^(-3 - 2p) * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9)
	denominator = 8 * base_denom * sec_denom

	return numerator / denominator
end

function stagnation_surface(metric, θp, p)
	ro = 4.0
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


function drift_velocity(metric, r, θ, θp, p)
	met_dd = Krang.metric_dd(metric, r, θ)
	α = lapse(metric, r, θ)
	#rp = Krang.horizon(metric)
	Σt = Krang.Σ(metric, r, θ)
	Δt = Krang.Δ(metric, r)
	Ξt = Krang.Ξ(metric, r, θ; Δ = Δt)
	ωt = Krang.ω(metric, r, θ; Ξ = Ξt)
	#ψ = magnetic_stream_function(rp, θp, p)
	ΩF = magnetic_field_line_angular_rotation(metric, θp, p)
	At = Krang.A(metric, r, θ; Δ = Δt)

	B_BL_u = magnetic_field_BL_u(metric, r, θ, p)
	B_BL_d = met_dd * B_BL_u
	E_BL_u = [0, (ΩF-ωt)*At*sin(θ)*Bθ/Σt, -(ΩF-ωt)*At*sin(θ)*Br/(Σt*Δt), 0]
	E_BL_d = met_dd * E_BL_u

	Er_BL_d, Eθ_BL_d, Eϕ_BL_d = E_BL_d
	Br_BL_d, Bθ_BL_d, Bϕ_BL_d = B_BL_d

	B2 = B_BL_u * B_BL_d
	g = sqrt(-Krang.metric_det(metric, r, θ))


	return (α/(B2*g)) .* [(Eθ_BL_d*Bϕ_BL_d - Bθ_BL_d*Eϕ_BL_d), (Eϕ_BL_d*Br_BL_d - Bϕ_BL_d*Er_BL_d), (Er_BL_d*Bθ_BL_d - Br_BL_d*Eθ_BL_d)]

end

function electric_field_BL_u(metric, r, θ, B_BL_u)

	Δt = Krang.Δ(metric, r)
	Σt = Krang.Σ(metric, r, θ)
	Ξt = Krang.Ξ(metric, r, θ; Δ = Δt)

	ωt = Krang.ω(metric, r, θ; Ξ = Ξt)
	_, Br, Bθ, _ = B_BL_u
	return [0, (ΩF-ωt)*Ξt*sin(θ)*Bθ/Σt, -(ΩF-ωt)*Ξt*sin(θ)*Br/(Σt*Δt), 0]

end
function velocity_ff_BL_u(metric, r, θp, p, B_BL_u, E_BL_u)
	rp = Krang.horizon(metric)
	ψ = rp^p*(1-cos(θp))
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
	ξ = (h*f + pm*eco*sqrt(eco^2-f^2+h^2))/(eco^2+h^2)
	γ = γo/√(1-ξ^2)
	return γ*(n_BL_u + v_drift_BL_u + ξ/γo*Bhat_BL_u)
end
