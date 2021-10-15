export rad_type, varRad, fixRad
export calc_surf_RAD, calc_cloudtop_RAD

## create type for radiation
## one where ΔR is prescribed
## one where ΔR is calculated from LWP and cloud top temperature
abstract type rad_type end
struct varRad <: rad_type end
struct fixRad <: rad_type end

"""
    calculate net SW and LW radiation at the surface
"""
function calc_surf_RAD(u,p)
    zi, hM, qM, SST, CF = u;

    # shortwave calculation
    LWP = incloud_LWP(u)*1000.0; # kg/m^2 \to g/m^2
    αc = cloud_albedo(LWP, CF);
    αs = 0.1; # surface albedo of ocean water
    SW_net = (1-αc) * (1-αs) * S_subtr;

    # longwave calculation
    ϵc_down = cloud_emissivity(LWP, CF);
    zb = calc_LCL(u);
    zc = (zi+zb)/2.0;
    Tc = temp(zc,hM,qM);
    LW_down = σ_SB * Tc^4;
    LW_up = σ_SB * SST^4;

    # sum them up
    RAD = SW_net + LW_down - LW_up;
    return RAD
end

"""
    returns the prescribed cloud-top radiative cooling ΔR
"""
function calc_cloudtop_RAD(u,p,rtype::fixRad)
    return p.ΔR
end

"""
    calculate the net ΔR at cloud-top based on CO2
    balance between upwelling and downwelling longwave
    downwelling longwave is based on an effective temperature
    which is empirically fit to LES

    gives ΔR ≈ 80 W/m2 for 400 ppm CO2
"""
function calc_cloudtop_RAD(u,p,rtype::varRad)
    zi, hM, qM, SST, CF = u;
    Tct = temp(zi,hM,qM);
    LWP = incloud_LWP(u)*1e3; # kg/m^2 \to g/m^2
    ϵc_up = cloud_emissivity(LWP, CF);
    
    Teff = 263.5 + 10.8*log(p.CO2/400.0);
    
    ΔR = ϵc_up * σ_SB * Tct^4 - σ_SB * Teff^4;
    return ΔR
end

"""
    albedo of the cloud given LWP in g/m^2
    cloud albedo from Stephens 1978 part 2. eq 1 and 7.
    backscatter β = 0.07, looked up from table 2.
    zenith angle θ = 60° and effective radius r_e = 10 um.
"""
function cloud_albedo(LWP, CF)
    # θ = 60.0; # degrees
    # r_e = 10.0; # um
    # β = 0.07; 
    # A = (2 * cosd(θ) * r_e)/ (3 * β);
    # αc = LWP / (A + LWP); 

    m = 0.795
    Lx = 19.136
    αc = m * (1 - Lx/(Lx+LWP));
    αc = αc * CF;
    return αc
end

"""
    cloud_emissivity(LWP)

    emissivity of the cloud as a function of LWP
    ϵ = 1 - exp(-a0 * LWP) with a0 = 0.15 m^2/given

    based on Stephens 1978 part II: eq 15 and 16
"""
function cloud_emissivity(LWP, CF)
    a0 = 0.15; # m^2/g
    ϵc = 1 - exp(-a0 * LWP); 
    ϵc = ϵc * CF;
    return ϵc
end