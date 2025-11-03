function magnetic_stream_function(r, θ, p) # 2410.00954v1 Eq (11)
    return r^p*(1-abs(cos(θ)))
end

function lapse(metric, r, θ)
    Δt = Krang.Δ(metric, r)
    Ξt = Krang.Ξ(metric, r, θ; Δ = Δt)
    Σt = Krang.Σ(metric, r, θ)

    return -Ξt/(Σt*Δt)#Eq 3 2410.00954v1
end
function dψ_dθ(r, θ, p)
    signθ = θ > π/2 ? -1 : 1
    return signθ*r^p*sin(θ)
end
function dψ_dr(r, θ, p)
    return p*magnetic_stream_function(r, θ, p)/r
end
function Br(metric, r, θ, p)
    α = lapse(metric, r, θ)
    g = Krang.met_det(metric, r, θ)
    return α*p*dψ_dθ(r,θ,p)/(√(-g)*r)
end
function Bθ(metric, r, θ, p)
    α = lapse(metric, r, θ)
    g = Krang.met_det(metric, r, θ)
    return -α*p*dψ_dr(r,θ,p)/(√(-g)*r)
end