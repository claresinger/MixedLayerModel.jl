export sst_type, fixSST, varSST
export mlm

## create type for surface energy balance
## one where SST is fixed
## one where SST evolves according to an energy budget
abstract type sst_type end
struct fixSST <: sst_type end
struct varSST <: sst_type end

"""
    dzidt(u, p, ent)

    evolution of inversion height, zi
    balance between entrainment and subsidence
"""
function dzidt(u, p, ent)
    zi, sM, qM, SST, CF = u;
    dzidt = ent - p.D*zi;
    return dzidt
end

"""
    dsMdt(u, p, ent, LWP) 

    evolution of mixed-layer enthalpy, sM
    negative of the vertical energy flux
"""
function dsMdt(u, p, ent, LWP)
    zi, sM, qM, SST, CF = u;
    ΔR = calc_cloudtop_RAD(u, p, LWP, p.rtype);
    S0 = S_0(u, p, p.ftype);
    Szi = S_zi(u, p, ent, LWP);
    dsMdt = -(1/zi) * (Szi - S0 + ΔR/ρref(SST));
    return dsMdt
end

"""
    dqMdt(u, p, ent)

    evolution of mixed-layer total water specific humidity, qM
    negative of the vertical water flux
"""
function dqMdt(u, p, ent, LWP)
    zi, sM, qM, SST, CF = u;
    Q0 = Q_0(u, p, p.ftype);
    Qzi = Q_zi(u, p, ent, LWP);
    dqMdt = -(1/zi) * (Qzi - Q0);
    return dqMdt
end

"""
    dSSTdt(u, p, LWP)

    defined as 0 for fixSST
"""
function dSSTdt(u, p, LWP, stype::fixSST)
    return 0
end

"""
    dSSTdt(u, p, LWP)

    close surface energy budge for varSST
"""
function dSSTdt(u, p, LWP, stype::varSST)
    zi, sM, qM, SST, CF = u;
    RAD = calc_surf_RAD(u, p, LWP);
    SHF = calc_SHF(u, p);
    LHF = calc_LHF(u, p);   
    c = ρw * Cw * p.Hw;
    return (1/c) * (RAD - SHF - LHF - p.OHU)
end

"""
    calculation cloud fraction 
    determined as a logistic function of S, the stability parameter
"""
function dCFdt(u, p, zb, LWP)
    zi, sM, qM, SST, CF = u;
    CFnew = cloud_fraction(u, p, zb, LWP);
    τ_CF = 3600.0*24.0*1.0; # 1 days; cloud fraction adjustment timescale [seconds]
    dCFdt = (CFnew - CF) / τ_CF;
    return dCFdt
end

"""
    mlm(du, u, p, t)    

    define the coupled ODE
      dzi/dt = we - D*zi
      dsM/dt = -dE/dz = 1/zi * (Szi - S0 + dR/rho)
      dqM/dt = -dW/dz = 1/zi * (Qzi - Q0)
      dSST/dt = 1/c * (SWnet - LWnet - SHF - LHF - OHU)
"""
function mlm(du, u, p, t)
    # println(t/3600/24)
    # println(u)
    zb = calc_LCL(u);
    LWP = incloud_LWP(u, zb);
    ent = we(u, p, zb, LWP, p.etype);
    du[1] = dzidt(u, p, ent)
    du[2] = dsMdt(u, p, ent, LWP)
    du[3] = dqMdt(u, p, ent, LWP)
    du[4] = dSSTdt(u, p, LWP, p.stype)
    du[5] = dCFdt(u, p, zb, LWP)
end