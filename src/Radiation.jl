export rad_type, varRad, fixRad
export surf_SW, surf_LW
export calc_surf_RAD, calc_cloudtop_RAD
export toa_net_rad, trop_sst
export test_trop_sst

export cloud_emissivity, Tatmos

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
    net downward shortwave radiation at surface
"""
function surf_SW(u,p)
    zi, hM, qM, SST, CF = u;

    # shortwave calculation
    LWP = incloud_LWP(u)*1e3; # kg/m^2 --> g/m^2
    αc = cloud_albedo(LWP, CF);
    SW_net = (1-αc) * (1-α_ocean) * S_subtr;
    return SW_net
end

"""
    net downward LW radiation at surface

    formerly written as the difference of two planck terms
    up = σ * SST^4
    down = CF(ϵ * σ * Tc^4) + (1-CF)(σ * Ta^4)

    now written as just a linear function of SST
"""
function surf_LW(u,p)
    zi, hM, qM, SST, CF = u;
    
    # longwave calculation
    # ϵc_down = cloud_emissivity(LWP);
    # zb = calc_LCL(u);
    # Tc = temp(zb,hM,qM);
    # Teff = Tatmos(p);
    # Ta = temp(zb/2.0,hM,qM);
    # LW_down = CF * (ϵc_down * σ_SB * Tc^4) + (1-CF) * (σ_SB * Teff^4);
    # LW_down = CF * (ϵc_down * σ_SB * Tc^4) + (1-CF) * (σ_SB * Ta^4);
    # LW_up = σ_SB * SST^4;
    # LW_net = LW_down - LW_up;
    
    # fit from LES experiments
    a = -30.0;
    b = 0.0;
    LW_net = a - b * (SST-300.0);
    return LW_net
end

"""
    calculate net SW and LW radiation at the surface
"""
function calc_surf_RAD(u,p)
    SW_net = surf_SW(u,p);
    LW_net = surf_LW(u,p);
    RAD = SW_net + LW_net;
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
    LWP = incloud_LWP(u)*1e3; # kg/m^2 --> g/m^2
    ϵc_up = cloud_emissivity(LWP);
    Teff = Tatmos(p);
    ΔR = CF * σ_SB * (ϵc_up * Tct^4 - Teff^4);
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
function cloud_emissivity(LWP)
    a0 = 0.15; # m^2/g
    ϵc = 1 - exp(-a0 * LWP); 
    return ϵc
end

"""
    toa_net_rad(u)
    calculates net sw and OLR at TOA
    OLR is linear func of SST based on LES
"""
function toa_net_rad(u)
    zi, hM, qM, SST, CF = u;

    T = 0.8;
    LWP = incloud_LWP(u)*1e3; # kg/m^2 \to g/m^2
    αc = cloud_albedo(LWP, CF);
    α = T*CF*αc + (1-CF)*α_ocean;
    SW_net = (1-α)*S_subtr;
    
    OLR = -491.0 + 2.57*SST;
    R_s = SW_net - OLR;
    return R_s
end

"""
    calc_R_s_400(u,p)
    if CO2 = 400, returns toa_net_rad
    else returns the saved parameter value
"""
function calc_R_s_400(u, p)
    if p.CO2 == 400.0
        x = toa_net_rad(u);
    else
        x = p.R_s_400;
    end
    return x
end

"""
    trop_sst(u,p)
    1. calculates subtropical TOA net SW and OLR
    2. calculates subtropical TOA radiative imbalance (relative to 400ppm)
    3. translates that to tropical TOA imbalance
    4. calculates tropical emission temp from tropical OLR
    5. calculates emission height as func of CO2 and H2O concentration
"""
function trop_sst(u, p)
    # net TOA imbalance
    R_s = toa_net_rad(u);
    ΔR_s = R_s - calc_R_s_400(u, p);
    ΔR_t = -p.AreaFrac/(1-p.AreaFrac)*ΔR_s;
    # emission height parameterization
    thermo_x = log((Rd/Rv)*(e0/psurf)*p.RHtrop) + (L0/Rv)*(1/T0);
    A = 1292.5;
    B = 886.8;
    He(Ts, CO2) = A * log(CO2) + B * thermo_x + B * (-L0/Rv) / Ts;
    # 400ppm baseline temperatures, lapse rates, emission height
    Ts400 = 300.0; # K 
    Γm400 = Γm(Ts400, p.RHtrop);
    He400 = He(Ts400, 400);
    Te400 = Ts400 - Γm400*He400;
    # change in emission temperature
    ΔTe = -ΔR_t / (4*σ_SB*Te400^3);
    # # change in surface temeprature
    # ΔHe(ΔTs) = A * log(p.CO2/400) + B * (L0/Rv) / Ts400^2 * ΔTs;
    # ΔΓ(ΔTs) = Γm(Ts400+ΔTs, p.RHtrop) - Γm400;

    # # no water vapor feedback or lapse rate feedback
    # ΔHe(ΔTs) = A * log(p.CO2/400);
    # ΔΓ(ΔTs) = 0.0;
    
    # f(ΔTs) = ΔTs - (ΔTe + ΔΓ(ΔTs)*He400 + Γm400*ΔHe(ΔTs));
    # ΔTs = find_zero(f, (eltype(u)(-20.0), eltype(u)(20.0)), Bisection());
    
    ΔTs = 0.0;

    # ΔTs = 0.0;
    # if sign(f(-20)) == sign(f(20))
    #     println(f(-20), "\t", f(20));
    #     println("out of bounds");
    #     println(u);
    #     println(ΔTe);
    # else
    #     try
    #         ΔTs = find_zero(f, (eltype(u)(-20.0), eltype(u)(20.0)), Bisection());
    #     catch
    #         println(f(-20), "\t", f(20));
    #     end
    # end

    # absolute tropical SST
    sst_t = Ts400 + ΔTs;
    return sst_t
end

function test_trop_sst(ΔR_s, p)
    ΔR_t = -p.AreaFrac/(1-p.AreaFrac)*ΔR_s;
    # println(ΔR_t);
    println(p.CO2);

    RHtrop = 0.8;
    thermo_x = log((Rd/Rv)*(e0/psurf)*RHtrop) + (L0/Rv)*(1/T0);
    A = 1292.5;
    B = 886.8;
    He(Ts, CO2) = A * log(CO2) + B * thermo_x + B * (-L0/Rv) / Ts;

    # 400ppm baseline temperatures, lapse rates, emission height
    Ts400 = 300.0; # K 
    Γm400 = Γm(Ts400, RHtrop);
    He400 = He(Ts400, 400);
    Te400 = Ts400 - Γm400*He400;
    #println(Ts400, "\t", Te400, "\t", Γm400, "\t", He400);

    # change in emission temperature
    ΔTe = -ΔR_t / (4*σ_SB*Te400^3);
    println(ΔTe);
    
    # change in surface temeprature
    ΔHe(ΔTs) = A * log(p.CO2/400) + B * (L0/Rv) / Ts400^2 * ΔTs;
    ΔΓ(ΔTs) = Γm(Ts400+ΔTs, RHtrop) - Γm400;
    f(ΔTs) = ΔTs - (ΔTe + ΔΓ(ΔTs)*He400 + Γm400*ΔHe(ΔTs));
    ΔTs = find_zero(f, 1.0);
    println(ΔTs);

    #println(ΔΓ(ΔTs)*He400);

    # absolute tropical SST
    sst_t = Ts400 + ΔTs;
    println(sst_t);
end