export ρref, pres, q_sat, q_v, q_l, temp, rho
export RH, theta
export calc_LWP, calc_LCL
export calc_qft0

"""
    ρref(T)    

    calculate reference density given temperature
    and reference pressure pref
"""
function ρref(T)
    return pref ./ (Rd .* T)
end

"""
    pres(z, T)

    assumes hydrostatic balance
"""
function pres(z, T)
    return psurf .* exp.((-g .* z) ./ (Rd .* T));
end

"""
    rho(z, T)    

    calculate density given altitude and temperature
"""
function rho(z, T)
    return pres(z,T) ./ (Rd .* T)
end

"""
    q_sat(z, T)

    uses Clasius-Clapeyron relation with assumed constant 
    latent heat of vaporization term L0=2.5e6
"""
function q_sat(z, T)
    psat = e0 .* exp.(-L0/Rv .* (1 ./ T .- 1/T0));
    qsat = Rd/Rv .* psat ./ (pres(z,T) .- psat);
    return qsat
end

"""
    q_v(z, T, qt)

    as the minimum between the total specific humidty and saturation specific humidity
"""
function q_v(z, T, qt)
    return min.(qt, q_sat(z,T));
end

"""
    q_l(z, T, qt)

    as difference between total specific humidity and saturation specific humidity
    if undersaturated, then ql=0
"""
function q_l(z, T, qt)
    return max.(qt - q_sat(z,T), 0.0);
end

"""
    temp(z, h, qt)

    uses saturation adjustment on the enthalpy
    if no zero is found, then set temp = 0°C
"""
function temp(z, h, qt)
    h_act(T) = Cp .* T .+ g .* z .+ L0 .* q_v(z,T,qt);
    f(x) = h - h_act(x);
    
    T_guess = 300.0;
    T = T0;
    try
        T = find_zero(f, T_guess);
    catch
        T = T0;
    end
        
    return T
end

"""
    RH(z, h, qt)

    relative humidity is the ratio of 
    total specific humidity to saturation
    max value is 1
"""
function RH(z, h, qt)
    qsat = q_sat(z, temp.(z, h, qt));
    x = min.(qt ./ qsat, 1.0);
    return x
end

"""
    theta(z,h,qt)

    calculate the potential temperature
"""
function theta(z,h,qt)
    T = temp.(z,h,qt);
    p = pres(z,T);
    θ = T .* (pref ./ p) .^ (Rd/Cp);
    return θ
end

"""
    calc_LCL(zi, hM, qtM)

    calculate the lifiting condensation level
"""
function calc_LCL(zi, hM, qtM)
    zb = zi;
    try
        f(z) = qtM - q_sat(z,temp(z,hM,qtM));
        zb = find_zero(f, (0.0,zi), Bisection());
    catch
        z0 = 0.0;
        ql0 = q_l(z0,temp(z0, hM, qtM),qtM);
        if ql0 > 0
            zb = z0;
        else
            zb = zi;
        end
    end
    return zb
end

"""
    calc_LWP(zi, hM, qtM)

    calulcate the liquid water path
"""
function calc_LWP(zi, hM, qtM)
    zb = calc_LCL(zi, hM, qtM);
    x = 0.0;
    dz = 1.0;
    for z in collect(zb:dz:zi)
        T = temp.(z,hM,qtM);
        ρ = rho(z,T);
        ql = q_l(z,T,qtM);
        x += ρ * ql * dz;
    end
    liq_wat_path = x;

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
    f(x) = x .- q_sat(zft, Tft(x)) .* RHft;
    qft0 = find_zero(f, (0.0,0.1), Bisection());
    qft0 = qft0 - Gamma_q * zft;
    return qft0
end