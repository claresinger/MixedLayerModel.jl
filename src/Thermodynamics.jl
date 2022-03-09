export ρref, pres, q_sat, q_v, q_l, temp, rho
export incloud_LWP, calc_LCL
export Γs, moist_adiabat, temp_ft
export calc_qft0

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
    q_v(z, T, qt)

    as the minimum between the total specific humidty and saturation specific humidity
"""
function q_v(z, T, qt)
    return min(qt, q_sat(z,T));
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
    temp(z, h, qt)

    uses saturation adjustment on the enthalpy
"""
function temp(z, h, qt)
    h_act(T) = Cp*T + g*z + L0*q_v(z,T,qt);
    f(x) = h - h_act(x);
    Tqt = (h - g*z - L0*qt) / Cp;
    return find_zero(f, eltype(h)(Tqt), Order1(), atol=0.1)
end

"""
    calc_LCL(u)

    calculate the lifiting condensation level
"""
function calc_LCL(u)
    zi, hM, qtM, SST, CF = u;

    f(z) = qtM - q_sat(z,temp(z, hM, qtM));
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
    zi, hM, qtM, SST, CF = u;

    dz = 1.0;
    z = zb:dz:zi;
    T = temp.(z,hM,qtM);
    ρ = rho.(z,T);
    ql = q_l.(z,T,qtM);
    return sum(ρ .* ql .* dz)
end

"""
    calc_qft0(RHft, Gamma_q, sft0, Gamma_s)

    calculate the initial free-tropospheric humidity given
    ft RH, ft humidity lapse rate, initial ft dry static energy, 
    and ft dry static energy lapse rate
"""
function calc_qft0(RHft, Gamma_q, sft0, Gamma_s)
    zft = 900.0;
    qft(x) = x + Gamma_q * zft;
    hft(x) = Cp * (sft0 + Gamma_s * zft) + L0 * qft(x);
    Tft(x) = temp(zft, hft(x), qft(x));
    f(x) = x - q_sat(zft, Tft(x)) * RHft;
    x0 = find_zero(f, (0.0,0.1), Bisection());
    return x0 - Gamma_q * zft;
end

"""
    Γs - saturated adiabatic lapse rate
"""
function Γs(z, T)
    rv = q_sat(z, T) / (1 - q_sat(z, T));
    return Γd * (1 + (L0*rv)/(Rd*T)) / (1 + (L0^2*rv)/(Rv*Cp*T^2))
end

"""
"""
function moist_adiabat(Tsurf, zft, p)
    qsurf = p.RHtrop * q_sat(0, Tsurf); # surface humidity
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
    return T
end

"""
    calculate actual moist adiabat by integrating
    go up dry adiabat to LCL and then saturated adiabat
"""
function temp_ft(Tsurf, zft, p)
    T = moist_adiabat(Tsurf, zft, p);
    return T[end]
end
