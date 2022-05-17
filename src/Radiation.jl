export rad_type, varRad, fixRad
export calc_surf_RAD, calc_cloudtop_RAD
export trop_sst

## create type for radiation
## one where ΔR is prescribed
## one where ΔR is calculated from LWP and cloud top temperature
abstract type rad_type end
struct varRad <: rad_type end
struct fixRad <: rad_type end

# """
#     free-tropospheric temperature
#     parameterized as a function of CO2 
#     fit to LES results
#     this is the effective emissions temperature
#     it gets warmer with increasing CO2 as the troposphere
#     gets optically thicker and the effective level of emissions
#     gets lower / closer to the cloud top 
# """
# function Tatmos(u, p, LWP)
#     zi, hM, qM, SST, CF = u;
#     SST_trop = trop_sst(u, p, LWP);
#     zft = 2000; # 2km
#     Tft = temp_ft(SST_trop, zft, p);
#     qft = p.RHft * q_sat(zft, Tft);
#     Ta = 375 + 2.8*log(p.CO2) + 21.1*log(qft); # co2 and qft
    
#     # Ta = 167 + 16.7*log(p.CO2); # just co2
    
#     return Ta
# end

# """
#     ΔTa(u, p, LWP, fttype::twocol)

#     writes the difference between cloud-top temperature and downwelling emission
#     temperature as a function of CO2 and H2O above-cloud
    
#     qft is taken at zft=1500m to be consistent with TopFluxes in the twocol framework
# """
# function ΔTa(u, p, LWP, fttype::twocol)
#     zi, hM, qM, SST, CF = u;
#     SST_trop = trop_sst(u, p, LWP);
#     zft = 1500;
#     Tft = temp_ft(SST_trop, zft, p);
#     qft = p.RHft * q_sat(zft, Tft);
#     ΔT = 16.0 + 3.0*log(p.CO2) + 8.9*log(qft); # co2 and qft
#     return ΔT
# end

"""
    ΔTa(u, p, LWP)

    writes the difference between cloud-top temperature and downwelling emission
    temperature as a function of CO2 and H2O above-cloud
"""
function ΔTa(u, p, LWP)
    zi, hM, qM, SST, CF = u;
    qft = qjump(u, p, LWP, p.fttype) + qM;
    qft = qft/2;
    ΔT = 16.0 + 3.0*log(p.CO2) + 8.9*log(qft); # co2 and qft
    return ΔT
end

"""
    calc_surf_RAD(u, p, LWP)

    calculate net SW and LW radiation at the surface
"""
function calc_surf_RAD(u, p, LWP)
    zi, hM, qM, SST, CF = u;

    # shortwave calculation
    αc = cloud_albedo(LWP);
    WVtrans = 1/exp(10*qM);
    SW_net = WVtrans * (1-CF*αc) * (1-α_ocean) * S_subtr;

    # direct greenhouse effect in subtropical clear-sky
    a0, a1, a2, b1, b2 = [12.4, -1020, 3.1, -270, 0.86];
    LW_net = (1-CF)*(a0*log(p.CO2/400) + a1 + a2*SST) + CF*(b1 + b2*SST);

    return SW_net + LW_net
end

"""
    calc_cloudtop_RAD(u, p, LWP, rtype::fixRad)

    returns the prescribed cloud-top radiative cooling ΔR
"""
function calc_cloudtop_RAD(u, p, LWP, rtype::fixRad)
    return p.ΔR
end

"""
    calc_cloudtop_RAD(u, p, LWP, rtype::varRad)

    calculate the net ΔR at cloud-top based on CO2
    balance between upwelling and downwelling longwave
    downwelling longwave is based on an effective temperature
    which is empirically fit to LES

    gives ΔR ≈ 75 W/m2 for 400 ppm CO2
"""
function calc_cloudtop_RAD(u, p, LWP, rtype::varRad)
    zi, hM, qM, SST, CF = u;
    Tct = temp(zi,hM,qM);
    #ϵc_up = cloud_emissivity(LWP);
    ϵc_up = 1.0;
    Teff = Tct + ΔTa(u, p, LWP);
    ΔR = CF * σ_SB * ϵc_up * (Tct^4 - Teff^4);
    ΔR = max(ΔR, 1.0);
    return ΔR
end

"""
    cloud_albedo(LWP)

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

"""
    trop_sst(u, p, LWP)

    tropical SST specified as a deviation from a base state p.Ts400
    with warming from export from the subtropics (proportional to all-sky albedo)
    and warming directly from GHG that depends on the ECS parameter
"""
function trop_sst(u, p, LWP)
    zi, hM, qM, SST, CF = u;

    # proportionality factor for 
    # tropical temperature increase 
    # relative to albedo decrease
    a_export = -0.1;
    αc0 = cloud_albedo(50e-3);
    CF0 = 1.0;
    Δαc = cloud_albedo(LWP) - αc0;
    ΔCF = CF - CF0;
    ΔT_export = a_export * (1-α_ocean) * S_trop/4 * (αc0 * ΔCF + CF0 * Δαc);

    # increase in tropical temperature from
    # direct greenhouse warming in tropics
    # ECS = °C per CO2 doubling 
    ΔT_greenhouse = p.ECS / log(2) * log(p.CO2 / 400);

    # T_trop = p.Ts400 + ΔT_export + ΔT_greenhouse;
    # T_trop = (p.Ts400 - 5) + ΔT_export + ΔT_greenhouse;
    # T_trop = (p.Ts400 - 5 - p.SST0) + SST;
    T_trop = p.Ts400 - 5;
    return T_trop
end
