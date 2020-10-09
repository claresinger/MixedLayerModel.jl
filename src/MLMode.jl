export sst_type, fixSST, varSST
export mlm

## create type for surface energy balance
## one where SST is fixed
## one where SST evolves according to an energy budget
abstract type sst_type end
struct fixSST <: sst_type end
struct varSST <: sst_type end

"""
    dzidt(u, p)

    evolution of inversion height, zi
    balance between entrainment and subsidence
"""
function dzidt(u, p)
    dzidt = we(u, p, p.etype) - p.D*u[1]
    return dzidt
end

"""
    dhMdt(u, p) 

    evolution of mixed-layer enthalpy, hM
    negative of the vertical energy flux
"""
function dhMdt(u, p)
    zi, hM, qM, SST = u;
    ΔR = calc_cloudtop_RAD(u,p,p.rtype);
    H0 = H_0(u, p, p.ftype);
    Hzi = H_zi(u, p);
    dEdz = 1/zi * (Hzi - H0 + ΔR/ρref(SST));
    dhMdt = - dEdz
    return dhMdt
end

"""
    dqMdt(u, p)

    evolution of mixed-layer total water specific humidity, qM
    negative of the vertical water flux
"""
function dqMdt(u, p)
    zi, hM, qM, SST = u;
    Q0 = Q_0(u, p, p.ftype);
    Qzi = Q_zi(u, p);
    dWdz = (1/zi) * (Qzi - Q0);
    dqMdt = - dWdz
    return dqMdt
end

"""
    dSSTdt(u, p)

    define dSSTdz() function
    defined as 0 for fixSST
"""
function dSSTdt(u, p, stype::fixSST)
    return 0
end

"""
    dSSTdt(u, p)

    define dSSTdz() function
    close surface energy budge for varSST
"""
function dSSTdt(u, p, stype::varSST)
    RAD = calc_surf_RAD(u,p);
    SHF = calc_SHF(u, p);
    LHF = calc_LHF(u, p);   
    c = ρw * Cw * p.Hw;
    x = (1/c) * (RAD - SHF - LHF - p.OHU);
    return x
end

"""
    mlm(du, u, p, t)    

    define the coupled ODE
"""
function mlm(du, u, p, t)
    du[1] = dzidt(u, p)
    du[2] = dhMdt(u, p)
    du[3] = dqMdt(u, p)
    du[4] = dSSTdt(u, p, p.stype)
end