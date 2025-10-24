function Krang._raytrace(
    observation,
    pix::Krang.AbstractPixel,
    mesh::Krang.Mesh{<:Krang.ConeGeometry{T,A},<:Krang.AbstractMaterial};
    res,
) where {T,A}
    geometry = mesh.geometry
    material = mesh.material
    θs = geometry.opening_angle
    subimgs = material.subimgs

    isindir = false
    opticaldepth = 0.0
    for _ = 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        for n in subimgs
            #νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
            intersection, issuccess = if Krang.isFastLight(material)
                if Krang.isAxisymmetric(material)
                    rs, νr, νθ, _, issuccess = Krang.emission_radius(pix, θs, isindir, n)
                    Krang.Intersection(zero(rs), rs, θs, zero(rs), νr, νθ), issuccess
                else
                    rs, ϕs, νr, νθ, issuccess =
                        Krang.emission_coordinates_fast_light(pix, θs, isindir, n)
                    Krang.Intersection(zero(rs), rs, θs, ϕs, νr, νθ), issuccess
                end
            else
                ts, rs, ϕs, νr, νθ, issuccess =
                    Krang.emission_coordinates(pix, θs, isindir, n)
                Krang.Intersection(ts, rs, θs, ϕs, νr, νθ), issuccess
            end
            if issuccess && (Krang.horizon(Krang.metric(pix)) < rs < T(Inf))
                source, curr_depth = upreferred.(material(pix, intersection))
                observation += source * (1 - exp(-curr_depth)) * (exp(-opticaldepth))
                opticaldepth += curr_depth
            end
        end
    end
    return observation
end
