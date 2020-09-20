module MLMparams

export interact_surf_params

@with_kw mutable struct interact_surf_params
    SST0::Real = 290.0; # (K)
    dSST::Real = 0.0; # (K/day)
    
    D::Real = 6.0e-6; # (1/s)
    CO2::Real = 400; # (ppm)
    
    RHsurf::Real = 0.80;
    RHft::Real = 0.25;
    
    Gamma_q::Real = -3e-6; # (kg/kg/m)
    sft0::Real = 297; # (K)
    Gamma_s::Real = 5e-3; # (K/m)
    qft0::Real = calc_qft0(RHft, Gamma_q, sft0, Gamma_s); # (kg/kg)
    
    A::Real = 2.0;
    
    V::Real = 10.0; # m/s
    CTh::Real = 8e-4;
    CTq::Real = 8e-4;
    
    Hw::Real = 1.0;
    OHU::Real = 10.0;
    
    etype::ent_type = bflux();
    ftype::flux_type = varFlux();
end

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

end