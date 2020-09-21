using Parameters
using Roots
using MixedLayerModel.Thermodynamics
using MixedLayerModel.Entrainment
using MixedLayerModel.SurfaceFluxes

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

# @with_kw mutable struct dycoms_params
#     SST0::Real = 292.5; # (K)
#     dSST::Real = 0.0; # (K/day)
    
#     D::Real = 3.75e-6; # (1/s)
#     DFR::Real = 48.0; # (W/m^2)
    
#     CTh::Real = 0.0008; #0.000945;
#     CTq::Real = 0.0008; #0.0013;
#     V::Real = 5.25; # (m/s)
    
#     RHsurf::Real = 0.70;
#     RHft::Real = 0.25;
    
#     Gamma_q::Real = -2e-6; # (kg/kg/m)
#     sft0::Real = 296; #294.75; # (K)
#     Gamma_s::Real = 5e-3; # (K/m)
#     qft0::Real = calc_qft0(RHft, Gamma_q, sft0, Gamma_s); # (kg/kg)
    
#     A::Real = 2.0;
#     a::Real = 1.0;

#     LHF::Real = 93.0; # (W/m^2) 93
#     SHF::Real = 16.0; # (W/m^2) 16, LHF*(μ-1.0), LHF*(μ-1.0) + DFR*(1-a)
    
#     etype::ent_type = bflux();
#     rtype::rad_type = direct();
#     ftype::flux_type = fixFlux();
# end

# @with_kw mutable struct cgils_params
#     SST0::Real = 290.0; # (K)
#     dSST::Real = 0.0; # (K/day)
    
#     D::Real = 6.0e-6; # (1/s)
#     DFR::Real = 60.0; # (W/m^2)
    
#     RHsurf::Real = 0.70;
#     RHft::Real = 0.25;
    
#     Gamma_q::Real = -2e-6; # (kg/kg/m)
#     sft0::Real = 296; #294.75; # (K)
#     Gamma_s::Real = 5e-3; # (K/m)
#     qft0::Real = calc_qft0(RHft, Gamma_q, sft0, Gamma_s); # (kg/kg)
    
#     A::Real = 2.0;

#     LHF::Real = 60.0; # (W/m^2)
#     SHF::Real = 10.0; # (W/m^2)
    
#     etype::ent_type = bflux();
#     rtype::rad_type = direct();
#     ftype::flux_type = fixFlux();
# end

# @with_kw mutable struct doomsday_params
#     SST0::Real = 290.0; # (K)
#     dSST::Real = 0.0; # (K/day)
    
#     D::Real = 3.0e-6; # (1/s)
#     DFR::Real = 74.0; # (W/m^2)
    
#     CTh::Real = 0.0008;
#     CTq::Real = 0.0008;
#     V::Real = 5.25; # (m/s)
    
#     RHsurf::Real = 0.70;
#     RHft::Real = 0.25;
    
#     Gamma_q::Real = -2e-6; # (kg/kg/m)
#     sft0::Real = 296; #294.75; # (K)
#     Gamma_s::Real = 5e-3; # (K/m)
#     qft0::Real = calc_qft0(RHft, Gamma_q, sft0, Gamma_s); # (kg/kg)
    
#     A::Real = 2.0;
#     a::Real = 1.0;

#     LHF::Real = 74.0; # (W/m^2)
#     SHF::Real = 2.0; # (W/m^2)
    
#     etype::ent_type = bflux();
#     rtype::rad_type = direct();
#     ftype::flux_type = fixFlux();
# end

# @with_kw mutable struct CPQ_params
#     SST0::Real = 290.0; # (K)
#     dSST::Real = 0.0; # (K/day)
    
#     D::Real = 6.0e-6; # (1/s)
#     DFR::Real = 74.0; # (W/m^2)
    
#     RHsurf::Real = 0.80;
#     RHft::Real = 0.25;
    
#     Gamma_q::Real = -3e-6; # (kg/kg/m)
#     sft0::Real = 297; # (K)
#     Gamma_s::Real = 5e-3; # (K/m)
#     qft0::Real = calc_qft0(RHft, Gamma_q, sft0, Gamma_s); # (kg/kg)
    
#     A::Real = 2.0;

#     LHF::Real = 107.0; # (W/m^2)
#     SHF::Real = 2.3; # (W/m^2)
    
#     etype::ent_type = bflux();
#     rtype::rad_type = direct();
#     ftype::flux_type = fixFlux();
# end