export flux_type, varFlux, fixFlux
export S_0, Q_0, calc_SHF, calc_LHF

## create type for surface fluxes
## one where SHF and LHF are calculated interactively from bulk aerodynamic formula
## one where SHF and LHF are prescribed (fixed)
abstract type flux_type end
struct varFlux <: flux_type end
struct fixFlux <: flux_type end

# adjustment factor based on degree of decoupling (or CF)
# if CF = p.CFmax then flux_adj = 1, if CF = p.CFmin then flux_adj = 1 + p.flux_α
function flux_adj(u,p)
    zi, sM, qM, SST, CF = u;
    return 1 + p.flux_α * ((p.CFmax - CF) / (p.CFmax - p.CFmin));
end

"""
    define surface liquid static energy flux, S_surf
    using bulk aerodynamic formula

    S_surf = C * V * (s0 - s)
"""
function S_0(u, p, ftype::varFlux)
    zi, sM, qM, SST, CF = u;
    s0 = Cp * SST;
    # S0 = p.Cd * p.V * (s0 - sM);
    S0 = p.Cd * p.V * flux_adj(u,p) * (s0 - sM);
    return S0
end

"""
    define the surface liquid static energy flux, S_surf
    given prescribed sensible and latent heat fluxes

    S_surf = SHF / ρref
"""
function S_0(u, p, ftype::fixFlux)
    zi, sM, qM, SST, CF = u;
    return p.SHF / ρref(SST)
end

"""
    define surface moisture flux, Q_surf
    using bulk aerodynamic formula

    Q_surf = C * V * (q0 - q)
"""
function Q_0(u, p, ftype::varFlux)
    zi, sM, qM, SST, CF = u;
    qs = q_sat(0.0,SST);
    # Q0 = p.Cd * p.V * (qs - qM);
    Q0 = p.Cd * p.V * flux_adj(u,p) * (qs - qM);
    return Q0
end

"""
    define the surface moisture flux, Q_surf
    given prescribed latent heat flux

    Q_surf = LHF / (Lv * ρref)
"""
function Q_0(u, p, ftype::fixFlux)
    zi, sM, qM, SST, CF = u;
    Q0 = p.LHF / (ρref(SST) * L0);
    return Q0
end

"""
    calculate the latent heat flux
"""
function calc_LHF(u, p)
    zi, sM, qM, SST, CF = u;
    return ρref(SST) * L0 * Q_0(u, p, p.ftype)
end

"""
    calculate the sensible heat flux
"""
function calc_SHF(u, p)
    zi, sM, qM, SST, CF = u;
    return ρref(SST) * S_0(u, p, p.ftype)
end