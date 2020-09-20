module Radiation

using ..Thermodynamics

export calc_surf_RAD, calc_cloud_RAD

"""
    calculate net SW and LW radiation at the surface
"""
function calc_surf_RAD(u,p)
    zi, hM, qM, SST = u;

    # shortwave calculation
    LWP = calc_LWP(zi, hM, qM)*1000.0; # kg/m^2 \to g/m^2
    αc = cloud_albedo(LWP);
    αs = 0.1; # surface albedo of ocean water
    SW_net = (1-αc) * (1-αs) * S_subtr;

    # longwave calculation
    ϵc_down = cloud_emissivity_down(LWP);
    zb = calc_LCL(zi,hM,qM);
    zc = (zi+zb)/2.0;
    Tc = temp(zc,hM,qM);
    LW_down = σ_SB * Tc^4;
    LW_up = σ_SB * SST^4;

    # sum them up
    RAD = SW_net + LW_down - LW_up;
    return RAD
end

"""
    calculate the net OLR at cloud-top based on CO2
    this is the 3-layer atmosphere model
"""
function calc_cloud_RAD(u,p)
    zi, hM, qM, SST = u;
    Tct = temp(zi,hM,qM);
    Ta = Tct - 5.0;
    ϵa = atmos_emissivity(p);
    
    LWP = calc_LWP(zi, hM, qM)*1000.0; # kg/m^2 \to g/m^2
    ϵc_up = cloud_emissivity_up(LWP);
    
    ΔR = ϵc_up * σ_SB * Tct^4 - ϵa * σ_SB * Ta^4;
    return ΔR
end

"""
    albedo of the cloud given LWP in g/m^2
    cloud albedo from Stephens 1978 part 2. eq 1 and 7.
    backscatter β = 0.07, looked up from table 2.
"""
function cloud_albedo(LWP)
    sza = 60.0; # degrees
    r_e = 10.0; # um
    backscatter = 0.07; 
    A = (2 * cosd(sza) * r_e)/ (3 * backscatter);
    αc = LWP / (A + LWP); 
    return αc
end

"""
    emissivity of the atmosphere given a concentration of CO2 in ppm
    function chosen to have eps=0.8 for 400ppm and 0.9 for 800ppm
    this yields reasonable ΔR values of order 50 W/m^2 for current temps
"""
function atmos_emissivity(p)
    #ϵa = exp(-1/p.CO2+1/400.0)*0.8;
    ϵa = log(p.CO2)/log(400.0)*0.8;
    return ϵa
end

"""
    emissivity of the cloud given LWP in g/m^2
    emissivity up from Stephens 1978 part 2. eq 15 and 16.
"""
function cloud_emissivity_up(LWP)
    a0_up = 0.13; # m^2/g
    ϵc_up = 1 - exp(-a0_up * LWP); 
    return ϵc_up
end

"""
    emissivity of the cloud given LWP in g/m^2
    emissivity down from Stephens 1978 part 2. eq 15 and 16.
"""
function cloud_emissivity_down(LWP)
    a0_down = 0.158; # m^2/g
    ϵc_down = 1 - exp(-a0_down * LWP); 
    return ϵc_down
end

end