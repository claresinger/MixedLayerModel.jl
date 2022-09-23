export rad_type, varRad, fixRad
export calc_surf_RAD, calc_cloudtop_RAD

## create type for radiation
## one where ΔR is prescribed
## one where ΔR is calculated from LWP and cloud top temperature
abstract type rad_type end
struct varRad <: rad_type end
struct fixRad <: rad_type end

"""
    ΔTa(u, p, LWP)
    writes the difference between cloud-top temperature and downwelling emission
    temperature as a function of CO2 and H2O above-cloud
"""
function ΔTa(u, p, LWP)
    zi, sM, qM, SST, CF = u;
    qft = qjump(u, p, LWP, p.fttype) + qM;
    ΔT = -10.1 + 3.1*log(p.CO2) + 5.3*log(qft);
    
    # TODO without proper twocol FT do this:
    # ΔT = -22.5 + 0.008*p.CO2;
    return ΔT
end

"""
    calc_surf_RAD(u, p, LWP)

    calculate net SW and LW radiation at the surface
"""
function calc_surf_RAD(u, p, LWP)
    zi, sM, qM, SST, CF = u;

    # shortwave calculation
    # WVtrans = exp(-10*qM); # TODO: WV abs coefficient
    # αc = cloud_albedo(LWP);
    # SW_net = WVtrans * (1 - (1-CF)*α_ocean - CF*αc) * S_subtr;
    # SW_net = (1 - (1-CF)*0.4 - CF*αc) * S_subtr;

    # LW_net linear with SST with coefficient dependent on log(CO2)
    # direct greenhouse effect in subtropical clear-sky
    # a0, a1, a2, b1, b2 = [12.4, -1020, 3.1, -270, 0.86];
    # LW_net = (1-CF)*(a0*log(p.CO2/400) + a1 + a2*SST) + CF*(b1 + b2*SST);
    
    # TODO simplify radiative fluxes and keep constant
    SW_net = p.SW_a + p.SW_b*(p.CFmax - CF);
    LW_net = -30;

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
    zi, sM, qM, SST, CF = u;
    Tct = temp(zi,sM,qM);
    ϵc_up = cloud_emissivity(LWP);
    Teff = Tct + ΔTa(u, p, LWP);
    ΔR = CF * σ_SB * ϵc_up * (Tct^4 - Teff^4);
    ΔR = max(ΔR, 1)
    return ΔR
end

"""
    cloud_albedo(LWP)

    albedo of the cloud given LWP in kg/m^2
    fit from LES experiments
"""
function cloud_albedo(LWP)
    αmax = 0.98;
    Lx = 36e-3;
    αc = αmax * (LWP)/(Lx + LWP);

    αc = 0.75
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

    ϵc = 0.9
    return ϵc
end
