export ent_type, fixed, enBal, bflux
export we

###########
# create structure for ent_type
###########
abstract type ent_type end
struct fixed <: ent_type end
struct enBal <: ent_type end
struct bflux <: ent_type end

"""
    we(u, p, etype::fixed)

    fixed entrainment velocity of 1.5 mm/s
""" 
function we(u, p, etype::fixed)
    w = 0.0015; # 1.5 mm/s
    return w
end

"""
    we(u, p, etype::enBal)

    entrainment velocity obtained via energy balance requirement
    w = ΔR / (Δs_vli * ρref)
"""
function we(u, p, etype::enBal)
    zi, hM, qM, SST = u;
    ΔR = calc_cloudtop_RAD(u,p,p.rtype);

    # calculate change in s_vl across inversion
    hft = h_ft(zi, p);
    qft = q_ft(zi, p);
    Δs_vli = (hft - hM) - μ*L0*(qft - qM);

    f = ΔR / rho_ref(SST);

    # calculate entrianment rate
    w = (f/Δs_vli);
    return w
end

"""
    we(u, p, etype::bflux)

    entrainment velocity based on buoyancy flux
    without radiation
    
    integral is calculated analytically
"""
function we(u, p, etype::bflux)
    zi, hM, qM, SST = u;
    
    # calculate change in s_vl across inversion
    hft = h_ft(zi, p);
    qft = q_ft(zi, p);
    Δs_vli = (hft - hM) - μ*L0*(qft - qM);

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
    
    α = (2.5 * p.A) / (zi * Δs_vli);
    w = α*I0 / (1 - α*I1)
    return w
end