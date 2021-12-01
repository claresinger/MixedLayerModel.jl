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
    dzidt = ent - p.D*u[1]
    return dzidt
end

"""
    dhMdt(u, p, ent) 

    evolution of mixed-layer enthalpy, hM
    negative of the vertical energy flux
"""
function dhMdt(u, p, ent, zb)
    zi, hM, qM, SST, CF = u;
    ΔR = calc_cloudtop_RAD(u, p, zb, p.rtype);
    H0 = H_0(u, p, p.ftype);
    Hzi = H_zi(u, p, ent, zb);
    dhMdt = -(1/zi) * (Hzi - H0 + ΔR/ρref(SST));
    return dhMdt
end

"""
    dqMdt(u, p, ent)

    evolution of mixed-layer total water specific humidity, qM
    negative of the vertical water flux
"""
function dqMdt(u, p, ent, zb)
    zi, hM, qM, SST, CF = u;
    Q0 = Q_0(u, p, p.ftype);
    Qzi = Q_zi(u, p, ent, zb);
    dqMdt = -(1/zi) * (Qzi - Q0);
    return dqMdt
end

"""
    dSSTdt(u, p)

    defined as 0 for fixSST
"""
function dSSTdt(u, p, zb, stype::fixSST)
    return 0
end

"""
    dSSTdt(u, p)

    close surface energy budge for varSST
"""
function dSSTdt(u, p, zb, stype::varSST)
    RAD = calc_surf_RAD(u, p, zb);
    SHF = calc_SHF(u, p);
    LHF = calc_LHF(u, p);   
    c = ρw * Cw * p.Hw;
    x = (1/c) * (RAD - SHF - LHF - p.OHU);
    
    zi, hM, qM, SST, CF = u;
    τ_SST = 3600.0*24.0*1.0;
    y = (p.SST0 - SST) / τ_SST;
    return (x + y)
end

"""
    calculation cloud fraction 
    determined as a logistic function of S, the stability parameter
"""
function dCFdt(u, p, zb)
    zi, hM, qM, SST, CF = u;
    CFnew = cloud_fraction(u, p, zb);
    τ_CF = 3600.0*24.0*1.0; # 1 days; cloud fraction adjustment timescale [seconds]
    dCFdt = (CFnew - CF) / τ_CF;
    return dCFdt
end

"""
    mlm(du, u, p, t)    

    define the coupled ODE
      dzi/dt = we - D*zi
      dhM/dt = -dE/dz = 1/zi * (Hzi - H0 + dR/rho)
      dqM/dt = -dW/dz = 1/zi * (Qzi - Q0)
      dSST/dt = 1/c * (SWnet - LWnet - SHF - LHF - OHU)
"""
function mlm(du, u, p, t)
    zb = calc_LCL(u);
    ent = we(u, p, zb, p.etype);
    du[1] = dzidt(u, p, ent)
    du[2] = dhMdt(u, p, ent, zb)
    du[3] = dqMdt(u, p, ent, zb)
    du[4] = dSSTdt(u, p, zb, p.stype)
    du[5] = dCFdt(u, p, zb)
end