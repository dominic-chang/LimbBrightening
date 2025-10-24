#Constants
Unitful.@unit c "c" SpeedofLight (2.998 * 10^8 * Unitful.m * Unitful.s^-1) false true
const G = u"G"
const me = u"me"
const mp = u"mp"
const e = u"q"
const kB = u"k"
const h = u"h"
const rg = upreferred(6.5 * 10^9 * UnitfulAstro.Msun * G / c^2)
const ke = upreferred(8.988 * 10^9 * u"N" * u"m^2/C^2")


function b_func_power(r, mag_field_magnitude, rb, p_b)
    """
        Magnetic field as a function of distance from the black hole
    """
    return mag_field_magnitude * (r / rb)^p_b
end


function te_func(r, t_e0, rb, p_temp)
    """
        Temperature as a function of distance from the black hole
    """
    return t_e0 * (r / rb)^p_temp
end


function theta_e_func(temp)
    """
        Dimensionless temperature value
    """
    return (kB * temp / (me * c^2))
end


function nth_func(r, n_th0, rb, p_dens)
    """
        Density at as a function of distance from the black hole:
    """
    return n_th0 * (r / rb)^p_dens
end


function nu_c_func(bsinθ_b, theta_e)
    """ 
        Frequnecy scalar
    """

    return 3 / 4π * e * bsinθ_b * theta_e^2 / (me)
end


function synchrotron_func_I(ν_νc)
    """
        Synchrotron emission function (A18 from MNRAS 462, 115–136 (2016))

        # Params
        -   `ν_νc` : Frequency ratio
    """
    return 2.5651 *
           (1 + 1.92 * (ν_νc^(-1 / 3)) + (0.9977 * ν_νc^(-2 / 3))) *
           exp(-1.8899 * ν_νc^(1 / 3))
end

function black_body_func(ν, temp)
    """
        Black body emission function
    """
    return 2 * h * ν^3 / (c^2) * 1 / (exp(upreferred(h * ν / (kB * temp))) - 1)
end

println("""

    Loaded Constants
    ----------------
    -   c : Speed of light ($(1c |> upreferred))
    -   G : Gravitational constant ($(1G |> upreferred))
    -   me : Electron mass ($(1me |> upreferred))
    -   mp : Proton mass ($(1mp |> upreferred))
    -   e : Elementary charge ($(1e |> upreferred))
    -   kB : Boltzmann constant ($(1kB |> upreferred))
    -   h : Planck constant ($(1h |> upreferred))
    -   rg : Schwarzschild radius ($(1rg |> upreferred))
    -   ke : Coulomb constant ($(1ke |> upreferred))
""")
