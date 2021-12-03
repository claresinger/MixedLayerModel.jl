export rad_type, varRad, fixRad
export calc_surf_RAD, calc_cloudtop_RAD
export toa_net_rad, trop_sst

export cloud_emissivity, Tatmos
export test_trop_sst

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

# function ΔTa(p)
#     return 22.5 - 0.008*p.CO2
# end

"""
    calculate net SW and LW radiation at the surface
"""
function calc_surf_RAD(u, p, LWP)
    zi, hM, qM, SST, CF = u;

    # shortwave calculation
    αc = cloud_albedo(LWP, CF);
    SW_net = (1-αc) * (1-α_ocean) * S_subtr;

    # longwave calculation
    # zi, hM, qM, SST, CF = u;
    # ϵc_down = cloud_emissivity(LWP);
    # Tc = temp(zi,hM,qM,SST);
    # Teff = Tatmos(p);
    # Ta = temp(zi/2.0,hM,qM,SST);
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
    Tct = temp(zi,hM,qM,SST);
    ϵc_up = cloud_emissivity(LWP);
    Teff = Tatmos(p);
    #Teff = Tct - ΔTa(p);
    ΔR = CF * σ_SB * ϵc_up * (Tct^4 - Teff^4);
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
    ϵ = 1 - exp(-a0 * LWP) with a0 = 0.15 m^2/g

    based on Stephens 1978 part II: eq 15 and 16
"""
function cloud_emissivity(LWP)
    a0 = 0.15 * 1e3; # m^2/kg
    ϵc = 1 - exp(-a0 * LWP); 
    return ϵc
end

"""
    toa_net_rad(u, LWP)
    calculates net sw and OLR at TOA
    OLR is linear func of SST based on LES
"""
function toa_net_rad(u, LWP)
    zi, hM, qM, SST, CF = u;

    T = 0.8;
    αc = cloud_albedo(LWP, CF);
    α = T*CF*αc + (1-CF)*α_ocean;
    SW_net = (1-α)*S_subtr;
    
    OLR = -491.0 + 2.57*SST;
    R_s = SW_net - OLR;
    return R_s
end

"""
    trop_sst(u, p, LWP)
    1. calculates subtropical TOA net SW and OLR
    2. calculates subtropical TOA radiative imbalance (relative to 400ppm)
    3. translates that to tropical TOA imbalance
    4. calculates tropical emission temp from tropical OLR
    5. calculates emission height as func of CO2 and H2O concentration
"""
function trop_sst(u, p, LWP)
    zi, hM, qM, SST, CF = u;

    # net TOA imbalance
    if p.CO2 == 400.0
        ΔR_s = 0.0;
    else
        ΔR_s = toa_net_rad(u, LWP) - p.R_s_400;
    end
    ΔR_t = -p.AreaFrac/(1-p.AreaFrac)*ΔR_s;
    
    # emission height parameterization
    thermo_x = log((Rd/Rv)*(e0/psurf)*p.RHtrop) + (L0/Rv)*(1/T0);
    A = 1292.5;
    B = 886.8;
    He(Ts, CO2) = A * log(CO2) + B * thermo_x + B * (-L0/Rv) / Ts;
    
    # change in emission temperature
    Ts400 = 300.0; # K 
    Te400 = Ts400 - Γm(Ts400, p.RHtrop)*He(Ts400, 400);
    ΔTe = -ΔR_t / (4*σ_SB*Te400^3);
    
    # # water vapor and lapse rate feedback on
    # ΔHe(ΔTs) = A * log(p.CO2/400) + B * (L0/Rv) / Ts400^2 * ΔTs;
    # ΔΓ(ΔTs) = Γm(Ts400+ΔTs, p.RHtrop) - Γm(Ts400, p.RHtrop);
    # f(ΔTs) = ΔTs - (ΔTe + ΔΓ(ΔTs)*He(Ts400, 400) + Γm(Ts400, p.RHtrop)*ΔHe(ΔTs));
    # ΔTs = find_zero(f, (-20, 20), Bisection());

    # no water vapor feedback or lapse rate feedback
    ΔHe = A * log(p.CO2/400.0);
    ΔTs = ΔTe + Γm(Ts400, p.RHtrop)*ΔHe;
    
    # absolute tropical SST
    #sst_t = Ts400 + ΔTs;

    # fix tropical SST right now 
    sst_t = SST + 10.0;

    return sst_t
end

# function test_trop_sst(ΔR_s, p)
#     ΔR_t = -p.AreaFrac/(1-p.AreaFrac)*ΔR_s;
#     println(p.CO2)
    
#     # emission height parameterization
#     thermo_x = log((Rd/Rv)*(e0/psurf)*p.RHtrop) + (L0/Rv)*(1/T0);
#     A = 1292.5;
#     B = 886.8;
#     He(Ts, CO2) = A * log(CO2) + B * thermo_x + B * (-L0/Rv) / Ts;
    
#     # change in emission temperature
#     Ts400 = 300.0; # K 
#     Te400 = Ts400 - Γm(Ts400, p.RHtrop)*He(Ts400, 400);
#     ΔTe = -ΔR_t / (4*σ_SB*Te400^3);
    
#     # water vapor and lapse rate feedback on
#     ΔHe(ΔTs) = A * log(p.CO2/400) + B * (L0/Rv) / Ts400^2 * ΔTs;
#     ΔΓ(ΔTs) = Γm(Ts400+ΔTs, p.RHtrop) - Γm(Ts400, p.RHtrop);

#     # # no water vapor feedback or lapse rate feedback
#     # ΔHe(ΔTs) = A * log(p.CO2/400);
#     # ΔΓ(ΔTs) = 0.0;
    
#     f(ΔTs) = ΔTs - (ΔTe + ΔΓ(ΔTs)*He(Ts400, 400) + Γm(Ts400, p.RHtrop)*ΔHe(ΔTs));
#     ΔTs = find_zero(f, (-20, 20), Bisection());
    
#     # absolute tropical SST
#     sst_t = Ts400 + ΔTs;
#     println(sst_t)
# end