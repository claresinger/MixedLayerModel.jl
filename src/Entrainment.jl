export ent_type, fixed, enBal, Sally, bflux
export we

###########
# create structure for ent_type
###########
abstract type ent_type end
struct fixed <: ent_type end
struct enBal <: ent_type end
struct Sally <: ent_type end
struct bflux <: ent_type end

"""
    we(u, p, etype::fixed)

    fixed entrainment velocity of 7 mm/s
""" 
function we(u, p, zb, etype::fixed)
    w = 0.007; # 7 mm/s
    return w
end

"""
    we(u, p, etype::enBal)

    entrainment velocity obtained via energy balance requirement
    w = ΔR / (Δs_vli * ρref)
"""
function we(u, p, zb, etype::enBal)
    zi, hM, qM, SST, CF = u;
    ΔR = calc_cloudtop_RAD(u, p, zb, p.rtype);

    # calculate change in s_vl across inversion
    hj = hjump(u, p, zb, p.fttype);
    qj = qjump(u, p, zb, p.fttype);
    Δs_vli = hj - μ*L0*qj;

    f = ΔR / ρref(SST);

    # calculate entrianment rate
    w = (f/Δs_vli);
    return w
end

"""
    we(u, p, etype::Sally)

    entrainment velocity obtained via energy balance requirement
    w = a * ΔR / (Δs_vli * ρref)
"""
function we(u, p, zb, etype::Sally)
    zi, hM, qM, SST, CF = u;
    ΔR = calc_cloudtop_RAD(u, p, zb, p.rtype);

    # calculate change in s_vl across inversion
    hj = hjump(u, p, zb, p.fttype);
    qj = qjump(u, p, zb, p.fttype);
    Δs_vli = hj - μ*L0*qj;

    f = ΔR / ρref(SST);

    # calculate entrianment rate
    w = p.a * (f/Δs_vli);
    return w
end

"""
    we(u, p, etype::bflux)

    entrainment velocity based on buoyancy flux
    without radiation
    
    integral is calculated analytically
"""
function we(u, p, zb, etype::bflux)
    zi, hM, qM, SST = u;
    
    # calculate change in s_vl across inversion
    hj = hjump(u, p, zb, p.fttype);
    qj = qjump(u, p, zb, p.fttype);
    Δs_vli = hj - μ*L0*qj;

    # calculate sv flux <w'sv'>(z) = f0(z) + we*f1(z)
    # split into two terms I0, I1 where Ii = \integral fi(z) dz
    H0 = H_0(u, p, p.ftype);
    Q0 = Q_0(u, p, p.ftype);
    zb = calc_LCL(u);

    A0 = H0 - μ*L0*Q0;
    B0 = β*H0 - ϵ*L0*Q0;
    I0 = A0 * (zb - (zb^2)/(2*zi)) + B0 * ((zi-zb) + (zi^2 - zb^2)/(2*zi));
    
    A1 = hj - μ*L0*qj;
    B1 = β*hj- ϵ*L0*qj;
    I1 = A1 * (-(zb^2)/(2*zi)) + B1 * ((-zi^2 + zb^2)/(2*zi));
    
    α = (2.5 * p.A) / (zi * Δs_vli);
    w = α*I0 / (1 - α*I1)
    return w
end