export upCO2, climatology

@with_kw mutable struct upCO2
    # radiation
    SW_a::Real = 120; # (W/m2)
    SW_b::Real = 150; # (W/m2)

    # cloud fraction params
    CFmax::Real = 1.0;      # maximum cloud fraction
    CFmin::Real = 0.2;      # minimum cloud fraction
    Dslope::Real = 8;       # decoupling slope
    Dcrit::Real = 0.9;      # critical decoupling parameter        

    # fixed SST for 400ppm
    SST0::Real = 290.0; # (K)

    # fix WV in radiation
    qft_rad::Real = 2e-3;

    # baseline CO2
    CO2::Real = 400; # (ppm)
    τCF::Real = 2; # (days)
    
    # subsidence strength
    D::Real = 6.0e-6; # (1/s)
    α_vent::Real = 1.0e-3; # (m/s)
    
    # params for interactive surface fluxes
    V::Real = 10.0; # (m/s)
    Cd::Real = 8e-4;
    
    # slab ocean params
    Hw::Real = 0.2; # (m), 0.2m ~ 10 days
    OHU::Real = 10.0; # (W/m2)

    # tropical column params
    RHtrop0::Real = 0.8;
    RHft::Real = 0.2;
    dTdz::Real = -5e-3; # (decreases by 5 K/km)

    EIS0::Real = 8.0; # (K)
    ECS::Real = 1.5; # (K / co2 doubling)
    Eexport::Real = 10.0; # (K)

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
    CFmax::Real = 0.8;          # maximum cloud fraction
    CFmin::Real = 0.05;         # minimum cloud fraction
    Dslope::Real = 8;       # decoupling slope
    Dcrit::Real = 0.9;      # critical decoupling parameter 

    # baseline CO2
    CO2::Real = 400; # (ppm)
    τCF::Real = 2; # (days)
    
    # subsidence strength
    D::Real = 6.0e-6; # (1/s)
    α_vent::Real = 1.0e-3; # (m/s)
    
    # params for fixedFT inverson specification
    sft0::Real = 300*Cp; # (K * J/K/kg), Tft0 = 300 (K)
    Gamma_s::Real = Cp*-5e-3 + g; # (K/m * J/K/kg), dT/dz=-5 K/km
    RHft::Real = 0.25;
    
    # for fixEIS inversion
    EIS0::Real = 8.0; # (K)
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
