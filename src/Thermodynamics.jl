export ρref, pres, q_sat, q_v, q_l, temp, rho
export incloud_LWP, calc_LCL
export Γs, moist_adiabat, temp_ft

"""
    ρref(T)    

    calculate reference density given temperature
    and reference pressure pref
"""
function ρref(T)
    return pref / (Rd * T)
end

"""
    pres(z, T)

    assumes hydrostatic balance
"""
function pres(z, T)
    return psurf * exp((-g * z) / (Rd * T));
end

"""
    rho(z, T)    

    calculate density given altitude and temperature
"""
function rho(z, T)
    return pres(z,T) / (Rd * T)
end

"""
    q_sat(z, T)

    uses Clasius-Clapeyron relation with assumed constant 
    latent heat of vaporization term L0=2.5e6
"""
function q_sat(z, T)
    psat = e0 * exp(-L0/Rv * (1 / T - 1/T0));
    return Rd/Rv * psat / (pres(z,T) - psat)
end

"""
    q_l(z, T, qt)

    as difference between total specific humidity and saturation specific humidity
    if undersaturated, then ql=0
"""
function q_l(z, T, qt)
    return max(qt - q_sat(z,T), 0.0);
end

"""
    temp(z, s, qt)

    uses saturation adjustment on the liquid water static energy (s)
"""
function temp(z, s, qt)
    s_act(T) = Cp*T + g*z - L0*q_l(z,T,qt);
    f(x) = s - s_act(x);
    Tguess = (s - g*z) / Cp;
    return find_zero(f, eltype(s)(Tguess), Order1(), atol=0.1)
end

"""
    calc_LCL(u)

    calculate the lifiting condensation level
"""
function calc_LCL(u)
    zi, sM, qtM, SST, CF = u;

    f(z) = qtM - q_sat(z,temp(z, sM, qtM));
    if f(0) > 0
        zb = 0.0;
    elseif f(zi) < 0
        zb = zi;
    else
        zb = find_zero(f, (0.0,zi), Bisection());
    end

    return zb
end

"""
    incloud_LWP(u)

    calulcate the in-cloud liquid water path
"""
function incloud_LWP(u, zb)
    zi, sM, qtM, SST, CF = u;
    if zb >= zi
        icLWP = 0.0
    else
        dz = 1.0;
        z = zb:dz:zi;
        T = temp.(z,sM,qtM);
        ρ = rho.(z,T);
        ql = q_l.(z,T,qtM);
        icLWP = sum(ρ .* ql .* dz);
    end
    return icLWP
end

"""
    Γs - saturated adiabatic lapse rate
"""
function Γs(z, T)
    rv = q_sat(z, T) / (1 - q_sat(z, T));
    return Γd * (1 + (L0*rv)/(Rd*T)) / (1 + (L0^2*rv)/(Rv*Cp*T^2))
end

"""
    moist_adiabat(Tsurf, zft, p)
    calculate moist adiabat given a surface temperature (Tsurf),
    up to an altitude zft, with the parameters p
    - first calculates the zLCL
    - then calculates the moist adiabatic profile with dz=10m up to zft
    returns (T,z) profile
"""
function moist_adiabat(Tsurf, zft, p)
    qsurf = p.RHtrop0 * q_sat(0, Tsurf); # surface humidity
    # find LCL
    f(x) = q_sat(x, Tsurf - x * Γd) - qsurf;
    if f(0) < 0
        zLCL = 0;
    elseif f(zft) > 0
        zLCL = zft;
    else 
        zLCL = find_zero(f, (0, zft), Bisection()); # LCL
    end

    # calculate moist adiabat
    dz = 10.0;
    z = range(0, zft, step=dz);
    T = zeros(length(z));
    for (i,zi) in enumerate(z)
        if zi <= zLCL # below LCL follow dry adiabat
            T[i] = Tsurf - zi*Γd;
        else # above LCL follow saturated adiabat
            T[i] = T[i-1] - Γs(zi, T[i-1])*dz;
        end
    end
    return T, z
end

"""
    calculate actual moist adiabat by integrating
    go up dry adiabat to LCL and then saturated adiabat
"""
function temp_ft(Tsurf, zft, p)
    if zft < 10
        println(zft)
    end
    T, z = moist_adiabat(Tsurf, zft, p);
    return T[end]
end
