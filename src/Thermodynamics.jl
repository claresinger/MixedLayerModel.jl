export ρref, pres, q_sat, q_v, q_l, temp, rho
export RH, theta
export incloud_LWP, calc_LCL
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
    qsat = Rd/Rv * psat / (pres(z,T) - psat);
    return qsat
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
    Tguess = eltype(h)((h - g*z - L0*q_v(z, 300.0, qt)) / Cp);
    T = find_zero(f, Tguess, Order1());

    # T = Tguess;
    # try
    #     T = find_zero(f, Tguess, Order1());#, xatol=0.1, atol=1.0);
    # catch
    #     println(z)
    #     println(h)
    #     println(qt)
    #     println("guess: ", Tguess)
    # end
    # println(z)
    # println("guess: ", Tguess)
    # println("success: ", T)
    return T
end

"""
    RH(z, h, qt)

    relative humidity is the ratio of 
    total specific humidity to saturation
    max value is 1
"""
function RH(z, h, qt)
    qsat = q_sat(z, temp(z, h, qt));
    x = min(qt / qsat, 1.0);
    return x
end

"""
    theta(z,h,qt)

    calculate the potential temperature
"""
function theta(z,h,qt)
    T = temp(z,h,qt);
    p = pres(z,T);
    θ = T * (pref / p) ^ (Rd/Cp);
    return θ
end

"""
    calc_LCL(u)

    calculate the lifiting condensation level
"""
function calc_LCL(u)
    zi, hM, qtM, SST, CF = u;

    f(z) = qtM - q_sat(z,temp(z,hM,qtM));
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
    liq_wat_path = sum(ρ .* ql .* dz);

    return liq_wat_path

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
    qft0 = find_zero(f, (0.0,0.1), Bisection());
    qft0 = qft0 - Gamma_q * zft;
    return qft0
end

"""
    Γm - moist adiabatic lapse rate calculation
"""
function Γm(Tsurf, RHsurf)
    qv = q_sat(0.0, Tsurf) * RHsurf;
    Γ = g * (1 + (L0*qv)/(Rd*Tsurf)) / (Cp + (L0^2*qv)/(Rv*Tsurf^2));
    return Γ
end