export upCO2, climatology

@with_kw mutable struct upCO2
    # fixed SST for 400ppm
    SST0::Real = 290.0; # (K)

    # baseline CO2
    CO2::Real = 400; # (ppm)
    ECS::Real = 3.0; # (K / co2 doubling)
    
    # subsidence strength
    D::Real = 6.0e-6; # (1/s)
    
    # params for interactive surface fluxes
    V::Real = 10.0; # (m/s)
    CTh::Real = 8e-4;
    CTq::Real = 8e-4;
    
    # slab ocean params
    Hw::Real = 1.0; # (m)
    OHU::Real = 10.0; # (W/m2)

    # tropical column params
    # AreaFrac::Real = 0.065;
    RHtrop0::Real = 0.8;
    Ts400::Real = 300.0; # (K)
    RHtropft::Real = 0.2;

    # fixed radiation
    ΔR::Real = 80.0; # (W/m2)

    # fixed surface fluxes
    SHF::Real = 10.0; # (W/m2)
    LHF::Real = 80.0; # (W/m2)

    # default types
    etype::ent_type = enBal();
    ftype::flux_type = varFlux();
    rtype::rad_type = varRad();
    stype::sst_type = fixSST();
    fttype::ft_type = twocol();
end

@with_kw mutable struct climatology
    # fixed SST for 400ppm
    SST0::Real = 290.0; # (K)

    # baseline CO2
    CO2::Real = 400; # (ppm)
    
    # subsidence strength
    D::Real = 6.0e-6; # (1/s)
    
    # params for fixedFT inverson specification
    RHft::Real = 0.25;
    Gamma_q::Real = -3e-6; # (kg/kg/m)
    sft0::Real = 300*Cp; # (K * J/K/kg), Tft0 = 300 (K)
    Gamma_s::Real = Cp*-5e-3 + g; # (K/m * J/K/kg), dT/dz=-5 K/km
    qft0::Real = 0.0; # (kg/kg)
    
    # params for interactive surface fluxes
    V::Real = 10.0; # (m/s)
    CTh::Real = 8e-4;
    CTq::Real = 8e-4;

    # fixed radiation
    ΔR::Real = 80.0;

    # fixed surface fluxes
    SHF::Real = 10.0;
    LHF::Real = 80.0;

    # default types
    etype::ent_type = enBal();
    ftype::flux_type = varFlux();
    rtype::rad_type = varRad();
    stype::sst_type = fixSST();
    fttype::ft_type = fixedFT();
end
