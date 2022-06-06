export ent_type
export fixed, enBal, bflux
export Δs,  we

###########
# create structure for ent_type
###########
abstract type ent_type end
struct fixed <: ent_type end
struct enBal <: ent_type end
struct bflux <: ent_type end

"""
    Δs(u, p, LWP)

    calculate the jump in dry virtual liquid static energy (s_vl)
    across the inversion
    Δs = Δh - μL Δq
"""
function Δs(u, p, LWP)
    # calculate change in s_vl across inversion
    hj = hjump(u, p, LWP, p.fttype);
    qj = qjump(u, p, LWP, p.fttype);
    sj = hj - μ*L0*qj
    return sj
end

"""
    we(u, p, zb, LWP, etype::fixed)

    fixed entrainment velocity of 7 mm/s
""" 
function we(u, p, zb, LWP, etype::fixed)
    w = 0.007; # 7 mm/s
    return w
end

"""
    we(u, p, zb, LWP, etype::enBal)

    entrainment velocity obtained via energy balance requirement
    w = ΔR / (Δs_vli * ρref)
"""
function we(u, p, zb, LWP, etype::enBal)
    zi, hM, qM, SST, CF = u;
    ΔR = calc_cloudtop_RAD(u, p, LWP, p.rtype);
    w = (ΔR / ρref(SST)) / Δs(u, p, LWP);
    return w
end

"""
    we(u, p, zb, LWP, etype::bflux)

    entrainment velocity based on buoyancy flux
    without radiation
    
    integral is calculated analytically
"""
function we(u, p, zb, LWP, etype::bflux)
    zi, hM, qM, SST = u;

    # calculate sv flux <w'sv'>(z) = f0(z) + we*f1(z)
    # split into two terms I0, I1 where Ii = \integral fi(z) dz
    H0 = H_0(u, p, p.ftype);
    Q0 = Q_0(u, p, p.ftype);

    A0 = H0 - μ*L0*Q0;
    B0 = β*H0 - ϵ*L0*Q0;
    I0 = A0 * (zb - (zb^2)/(2*zi)) + B0 * ((zi-zb) + (zi^2 - zb^2)/(2*zi));
    
    hj = hjump(u, p, LWP, p.fttype);
    qj = qjump(u, p, LWP, p.fttype);
    A1 = hj - μ*L0*qj;
    B1 = β*hj- ϵ*L0*qj;
    I1 = A1 * (-(zb^2)/(2*zi)) + B1 * ((-zi^2 + zb^2)/(2*zi));
    
    A = 2
    α = (2.5 * A) / (zi * Δs(u, p, LWP));
    w = α*I0 / (1 - α*I1)
    return w
end