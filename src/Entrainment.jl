export ent_type
export fixed, enBal, bflux
export sv_jump,  we

###########
# create structure for ent_type
###########
abstract type ent_type end
struct fixed <: ent_type end
struct enBal <: ent_type end
struct bflux <: ent_type end

"""
    sv_jump(u, p, LWP)

    calculate the jump in virtual liquid static energy (sv)
    across the inversion
    Δsv = Δs + CpΔTv - CpΔT
    Tv = (Rd*(1-qt) + Rv*qv)/Rd * T
"""
function sv_jump(u, p, LWP)
    zi, hM, qM, SST, CF = u;
    T_zi = (hM - g*zi - L0*qM)/Cp;
    ql_zi = q_l(zi, T_zi, qM);
    
    hft = hjump(u, p, LWP, p.fttype) + hM;
    qft = qjump(u, p, LWP, p.fttype) + qM;
    Tft = (hft - g*zi - L0*qft)/Cp;

    Tv_ft = (Rd*(1-qft) + Rv*qft)/Rd * Tft;
    Tv_M = (Rd*(1-qM) + Rv*(qM-ql_zi))/Rd * T_zi;
    Δsv = Cp*(Tv_ft - Tv_M);

    return Δsv
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
    w = ΔR / (Δsv * ρref)
"""
function we(u, p, zb, LWP, etype::enBal)
    zi, hM, qM, SST, CF = u;
    ΔR = calc_cloudtop_RAD(u, p, LWP, p.rtype);
    w = (ΔR / ρref(SST)) / sv_jump(u, p, LWP);
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
    α = (2.5 * A) / (zi * sv_jump(u, p, LWP));
    w = α*I0 / (1 - α*I1)
    return w
end