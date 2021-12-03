export rad_type, varRad, fixRad
export calc_surf_RAD, calc_cloudtop_RAD

## create type for radiation
## one where ΔR is prescribed
## one where ΔR is calculated from LWP and cloud top temperature
abstract type rad_type end
struct varRad <: rad_type end
struct fixRad <: rad_type end

"""
    free-tropospheric temperature
    parameterized as a function of CO2 
    fit to LES results
    this is the effective emissions temperature
    it gets warmer with increasing CO2 as the troposphere
    gets optically thicker and the effective level of emissions
    gets lower / closer to the cloud top 
"""
function Tatmos(p)
    return 263.5 + 10.8*log(p.CO2/400.0);
end

"""
    calculate net SW and LW radiation at the surface
"""
function calc_surf_RAD(u, p, LWP)
    zi, hM, qM, SST, CF = u;

    # shortwave calculation
    αc = cloud_albedo(LWP, CF);
    SW_net = (1-αc) * (1-α_ocean) * S_subtr;

    # longwave calculation
    # ϵc_down = cloud_emissivity(LWP);
    # Tc = temp(zi,hM,qM);
    # Teff = Tatmos(p);
    # Ta = temp(zi/2.0,hM,qM);
    # LW_down = CF * (ϵc_down * σ_SB * Tc^4) + (1-CF) * (σ_SB * Teff^4);
    # LW_down = CF * (ϵc_down * σ_SB * Tc^4) + (1-CF) * (σ_SB * Ta^4);
    # LW_up = σ_SB * SST^4;
    # LW_net = LW_down - LW_up;
    
    # fit from LES experiments
    LW_net = -30.0;

    return SW_net + LW_net
end

"""
    returns the prescribed cloud-top radiative cooling ΔR
"""
function calc_cloudtop_RAD(u, p, LWP, rtype::fixRad)
    return p.ΔR
end

"""
    calculate the net ΔR at cloud-top based on CO2
    balance between upwelling and downwelling longwave
    downwelling longwave is based on an effective temperature
    which is empirically fit to LES

    gives ΔR ≈ 80 W/m2 for 400 ppm CO2
"""
function calc_cloudtop_RAD(u, p, LWP, rtype::varRad)
    zi, hM, qM, SST, CF = u;
    Tct = temp(zi,hM,qM);
    ϵc_up = cloud_emissivity(LWP);
    Teff = Tatmos(p);
    # ΔR = CF * σ_SB * ϵc_up * (Tct^4 - Teff^4);
    ΔR = CF * σ_SB * (ϵc_up * Tct^4 - Teff^4);
    return ΔR
end

"""
    albedo of the cloud given LWP in kg/m^2
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

    m = 0.795;
    Lx = 19.136*1e-3;
    αc = m * (1 - Lx/(Lx+LWP));
    αc = αc * CF;
    return αc
end

"""
    cloud_emissivity(LWP)

    emissivity of the cloud as a function of LWP in kg/m2
    ϵ = 1 - exp(-a0 * LWP) with a0 = 0.15 m^2/given

    based on Stephens 1978 part II: eq 15 and 16
"""
function cloud_emissivity(LWP)
    a0 = 0.15; # m^2/g
    ϵc = 1 - exp(-a0 * LWP * 1e3); 
    return ϵc
end