export basic_params

@with_kw mutable struct basic_params
    # fixed SST for 400ppm
    SST0::Real = 290.0; # (K)

    # baseline CO2
    CO2::Real = 400; # (ppm)
    
    # baseline radiative cooling for non-interactive radiation
    ΔR::Real = 80; # (W/m2)
    
    # subsidence strength
    D::Real = 6.0e-6; # (1/s)
    
    # params for fixedFT inverson specification
    RHft::Real = 0.25;
    Gamma_q::Real = -3e-6; # (kg/kg/m)
    sft0::Real = 297; # (K) #297
    Gamma_s::Real = 5e-3; # (K/m)
    qft0::Real = 0.0; # (kg/kg)

    # params for SST dependent inversion specification
    hj_m::Real = -1.15*1e3 # J/kg/K
    hj_b::Real = -2.26*1e3 # J/kg
    qj_m::Real = -0.42*1e-3 # kg/kg/K
    qj_b::Real = -4.92*1e-3 # kg/kg
    
    a::Real = 0.8;
    A::Real = 2.0;
    
    # params for interactive surface fluxes
    RHsurf::Real = 0.70;
    V::Real = 10.0; # m/s
    CTh::Real = 8e-4;
    CTq::Real = 8e-4;

    # params for fixed surface fluxes
    LHF::Real = 90.0; # (W/m^2)
    SHF::Real = 5.0; # (W/m^2)
    
    # slab ocean params
    Hw::Real = 0.1;
    OHU::Real = 10.0;
    
    # default types
    etype::ent_type = enBal();
    ftype::flux_type = varFlux();
    rtype::rad_type = varRad();
    stype::sst_type = fixSST();
    fttype::ft_type = sstdep();
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
#     qft0::Real = 0.0; # (kg/kg)
    
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
#     qft0::Real = 0.0; # (kg/kg)
    
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
#     qft0::Real = 0.0; # (kg/kg)
    
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
#     qft0::Real = 0.0; # (kg/kg)
    
#     A::Real = 2.0;

#     LHF::Real = 107.0; # (W/m^2)
#     SHF::Real = 2.3; # (W/m^2)
    
#     etype::ent_type = bflux();
#     rtype::rad_type = direct();
#     ftype::flux_type = fixFlux();
# end