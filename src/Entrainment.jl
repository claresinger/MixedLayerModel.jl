module Entrainment

include("Definitions.jl")
using ..Thermodynamics
using ..SurfaceFluxes
using ..Radiation

export we

###########
# create structure for ent_type
# can add as many definitions for entrainment as you like!
# multiple dispatch is such a great thing. yay Julia!
###########
abstract type ent_type end
struct fixed <: ent_type end
struct enBal <: ent_type end
struct bflux <: ent_type end

# fixed entrainment velocity
function we(u, p, etype::fixed)
    w = 0.0015; # 1.5 mm/s
    return w
end

# entrainment velocity obtained via energy balance requirement
function we(u, p, etype::enBal)
    zi, hM, qM, SST = u;
    ΔR = calc_cloud_RAD(u,p);

    # calculate change in s_vl across inversion
    hft = h_ft(zi, p);
    qft = q_ft(zi, p);
    Ds_vli = (hft - hM) - μ*L0*(qft - qM);

    f = ΔR / rho_ref(SST);

    # calculate entrianment rate
    w = (f/Ds_vli);
    return w
end

# buoyancy flux entrainment velocity
# without radiation
function we(u, p, etype::bflux)
    zi, hM, qM, SST = u;
    
    # calculate change in s_v across inversion
    hft = h_ft(zi, p);
    qft = q_ft(zi, p);

    svlft = hft - μ*L0*qft;
    svlM = hM - μ*L0*qM;
    Ds_vli = svlft - svlM;

    # calculate sv flux <w'sv'>(z) = f0(z) + we*f1(z)
    # split into two terms I0, I1 where Ii = \integral fi(z) dz
    H0 = H_0(u, p, p.ftype);
    Q0 = Q_0(u, p, p.ftype);
    zb = calc_LCL(zi, hM, qM);

    A0 = H0 - μ*L0*Q0;
    B0 = β*H0 - ϵ*L0*Q0;
    I0 = A0 * (zb - (zb^2)/(2*zi)) + B0 * ((zi-zb) + (zi^2 - zb^2)/(2*zi));

    A1 = (hft - hM) - μ*L0*(qft - qM);
    B1 = β*(hft - hM) - ϵ*L0*(qft - qM);
    I1 = A1 * (-(zb^2)/(2*zi)) + B1 * ((-zi^2 + zb^2)/(2*zi));
    
    α = (2.5 * p.A) / (zi * Ds_vli);
    w = α*I0 / (1 - α*I1)
    return w
end

end