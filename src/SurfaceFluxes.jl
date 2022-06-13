export flux_type, varFlux, fixFlux
export S_0, Q_0, calc_SHF, calc_LHF

## create type for surface fluxes
## one where SHF and LHF are calculated interactively from bulk aerodynamic formula
## one where SHF and LHF are prescribed (fixed)
abstract type flux_type end
struct varFlux <: flux_type end
struct fixFlux <: flux_type end

"""
    define surface liquid static energy flux, S_surf
    using bulk aerodynamic formula

    S_surf = C * V * (s0 - s)
"""
function S_0(u, p, ftype::varFlux)
    zi, sM, qM, SST, CF = u;
    s0 = Cp * SST;
    return p.CTh * p.V * (s0 - sM)
end

"""
    define the surface enthalpy flux, H_surf
    given prescribed sensible and latent heat fluxes

    H_surf = (SHF + LHF) / ρref
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
    Q0 = p.CTq * p.V * (qs - qM);
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