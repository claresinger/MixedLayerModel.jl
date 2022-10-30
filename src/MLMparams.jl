export upCO2, climatology

@with_kw mutable struct upCO2
    # radiation
    SW_a::Real = 120; # (W/m2)
    SW_b::Real = 150; # (W/m2)

    # cloud fraction params
    CFmax::Real = 1.0;       # maximum cloud fraction
    CFmin::Real = 0.2;       # minimum cloud fraction
    decoup_slope::Real = 8; # decoupling slope, m

    # fixed SST for 400ppm
    SST0::Real = 290.0; # (K)

    # fix WV in radiation
    qft_rad::Real = 2e-3;

    # baseline CO2
    CO2::Real = 400; # (ppm)
    τCF::Real = 5; # (days)
    
    # subsidence strength
    D::Real = 6.0e-6; # (1/s)
    α_vent::Real = 2.0e-3; # (m/s)
    
    # params for interactive surface fluxes
    V::Real = 10.0; # (m/s)
    Cd::Real = 1e-3;
    
    # slab ocean params
    Hw::Real = 1.0; # (m)
    OHU::Real = 10.0; # (W/m2)

    # tropical column params
    RHtrop0::Real = 0.8;
    Ts400::Real = 300.0; # (K)
    RHft::Real = 0.2;
    ECS::Real = 4.0; # (K / co2 doubling)
    Eexport::Real = 15.0; # (K)

    EIS0::Real = 10.0; # (K)
    dTdz::Real = -5e-3; # (decreases by 5 K/km)

    # default types
    etype::ent_type = enBal();
    ftype::flux_type = varFlux();
    rtype::rad_type = varRad();
    stype::sst_type = fixSST();
    fttype::ft_type = co2EIS();
    wvtype::wvtype = wvON();
end

@with_kw mutable struct climatology
    # cloud fraction params
    CFmax::Real = 0.8;       # maximum cloud fraction
    CFmin::Real = 0.1;       # minimum cloud fraction
    decoup_slope::Real = 8; # decoupling slope (m)

    # baseline CO2
    CO2::Real = 400; # (ppm)
    
    # subsidence strength
    D::Real = 6.0e-6; # (1/s)
    α_vent::Real = 2.0e-3; # (m/s)
    
    # params for fixedFT inverson specification
    sft0::Real = 300*Cp; # (K * J/K/kg), Tft0 = 300 (K)
    Gamma_s::Real = Cp*-5e-3 + g; # (K/m * J/K/kg), dT/dz=-5 K/km
    RHft::Real = 0.25;
    
    # for fixEIS inversion
    EIS0::Real = 10.0; # (K)
    dTdz::Real = -5e-3; # (decreases by 5 K/km)
    
    # params for interactive surface fluxes
    SST0::Real = 290.; # (K)
    V::Real = 10.0; # (m/s)
    Cd::Real = 1e-3;
    Hw::Real = 1.0; # (m)
    OHU::Real = 10.0; # (W/m2)

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
    fttype::ft_type = fixedFT();
    wvtype::wvtype = wvON();
end
