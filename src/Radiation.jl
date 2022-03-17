export rad_type, varRad, fixRad
export calc_surf_RAD, calc_cloudtop_RAD
export trop_sst

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
function Tatmos(u, p, LWP)
    zi, hM, qM, SST, CF = u;
    SST_trop = trop_sst(u, p, LWP);
    zft = 3000; # 3km
    Tft = temp_ft(SST_trop, zft, p);
    qft = p.RHft * q_sat(zft, Tft);
    Ta = 375 + 2.8*log(p.CO2) + 21.1*log(qft); # co2 and qft
    
    # Ta = 167 + 16.7*log(p.CO2); # just co2
    
    return Ta
end

function ΔTa(u, p, LWP)
    zi, hM, qM, SST, CF = u;
    SST_trop = trop_sst(u, p, LWP);
    zft = 1000;
    Tft = temp_ft(SST_trop, zft, p);
    qft = p.RHft * q_sat(zft, Tft);
    ΔT = -7.92 + 2.51*log(p.CO2) + 4.97*log(qft); # co2 and qft   

    # ΔT = -22.5 + 0.008*p.CO2
    return ΔT
end

"""
    calculate net SW and LW radiation at the surface
"""
function calc_surf_RAD(u, p, LWP)
    zi, hM, qM, SST, CF = u;

    # shortwave calculation
    αc = cloud_albedo(LWP);
    WVtrans = 1/exp(10*qM);
    SW_net = WVtrans * (1-CF*αc) * (1-α_ocean) * S_subtr;

    # LW_net linear with SST with coefficient dependent on log(CO2)
    # greenhouse effect
    a0, b0, a1, b1 = [-1537.0,  4.88, -276.0, 0.88];
    LW_net = CF * (a0 + b0*SST) + (1-CF) * (a1 + b1*SST);

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
    Teff = Tatmos(u, p, LWP);
    # Teff = Tct + ΔTa(u, p, LWP);
    ΔR = CF * σ_SB * ϵc_up * (Tct^4 - Teff^4);
    ΔR += 1;
    return ΔR
end

"""
    albedo of the cloud given LWP in kg/m^2
    fit from LES experiments
    goes between 0 and 0.8
"""
function cloud_albedo(LWP)
    m = 0.795;
    Lx = 19.136*1e-3;
    αc = m * (1 - Lx/(Lx+LWP));
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

function trop_sst(u, p, LWP)
    zi, hM, qM, SST, CF = u;

    # proportionality factor for 
    # tropical temperature increase 
    # relative to albedo decrease
    a = 0.1;
    αc0 = cloud_albedo(0.1);
    CF0 = 1.0;
    Δαc = cloud_albedo(LWP) - αc0;
    ΔCF = CF - CF0;
    ΔT_export = -a * (1-α_ocean) * S_trop/4 * (αc0 * ΔCF + CF0 * Δαc);

    # increase in tropical temperature from
    # direct greenhouse warming in tropics
    # ECS = 3 K per doubling 
    ΔT_greenhouse = 3 / log(2) * log(p.CO2 / 400);

    return p.Ts400 + ΔT_export + ΔT_greenhouse
end

# """
#     toa_net_rad(u, LWP)
#     calculates net sw and OLR at TOA
#     OLR is linear func of SST based on LES
# """
# function toa_net_rad(u, LWP)
#     zi, hM, qM, SST, CF = u;

#     T = 0.8;
#     αc = cloud_albedo(LWP);
#     α = T*CF*αc + (1-CF)*α_ocean;
#     SW_net = (1-α)*S_subtr;
    
#     OLR = -491.0 + 2.57*SST;
#     R_s = SW_net - OLR;
#     return R_s
# end

# """
#     trop_sst(u, p, LWP)
#     1. calculates subtropical TOA net SW and OLR
#     2. calculates subtropical TOA radiative imbalance (relative to 400ppm)
#     3. translates that to tropical TOA imbalance
#     4. calculates tropical emission temp from tropical OLR
#     5. calculates emission height as func of CO2 and H2O concentration
# """
# function trop_sst(u, p, LWP)
#     zi, hM, qM, SST, CF = u;

#     # net TOA imbalance
#     if p.CO2 == 400.0
#         ΔR_s = 0.0;
#     else
#         ΔR_s = toa_net_rad(u, LWP) - p.R_s_400;
#     end
#     ΔR_t = -p.AreaFrac/(1-p.AreaFrac)*ΔR_s;
    
#     # emission height parameterization
#     thermo_x = log((Rd/Rv)*(e0/psurf)*p.RHtrop) + (L0/Rv)*(1/T0);
#     A = 1292.5;
#     B = 886.8;
#     He(Ts, CO2) = A * log(CO2) + B * thermo_x + B * (-L0/Rv) / Ts;
    
#     # change in emission temperature
#     Te400 = p.Ts400 - Γm(p.Ts400)*He(p.Ts400, 400);
#     #σTe3 = 4 * σ_SB * Te400^3;
#     σTe3 = 2.57;
#     ΔTe = -ΔR_t / σTe3;
    
#     # water vapor and lapse rate feedback on
#     ΔHe(ΔTs) = A * log(p.CO2/400) + B * (L0/Rv) / p.Ts400^2 * ΔTs;
#     ΔΓ(ΔTs) = Γm(p.Ts400+ΔTs) - Γm(p.Ts400);
#     f(ΔTs) = ΔTs - (ΔTe + ΔΓ(ΔTs)*He(p.Ts400, 400) + Γm(p.Ts400)*ΔHe(ΔTs));
#     ΔTs = find_zero(f, (-30, 30), Bisection());

#     # # no water vapor feedback or lapse rate feedback
#     # ΔHe = A * log(p.CO2/400.0);
#     # ΔTs = ΔTe + Γm(p.Ts400)*ΔHe;

#     # # only water vapor feedback on, lapse rate off
#     # ΔHe(ΔTs) = A * log(p.CO2/400) + B * (L0/Rv) / p.Ts400^2 * ΔTs;
#     # f(ΔTs) = ΔTs - (ΔTe + Γm(p.Ts400)*ΔHe(ΔTs));
#     # ΔTs = find_zero(f, (-30, 30), Bisection());
    
#     # absolute tropical SST
#     sst_t = p.Ts400 + ΔTs;

#     return sst_t
# end
