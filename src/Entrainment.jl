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
    we(u, p, zb, LWP, etype::fixed)

    fixed entrainment velocity of 7 mm/s
""" 
function we(u, p, zb, LWP, etype::fixed)
    w = 0.007; # 7 mm/s
    return w
end

"""
    sv_jump(u, p, LWP)

    Δsv = Δs + CpΔTv - CpΔT
        = Δs + Cp(Rv/Rd - 1)(Tft*qft - T_zi*qM) + Cp(Rv/Rd)(T_zi*ql_zi)

    jump in virtual liquid static energy across inversion
    proportional to buoyancy jump
    used in energy balance entrainment
"""
function sv_jump(u, p, LWP)
    zi, sM, qM, SST, CF = u;
    T_zi = temp(zi, sM, qM);
    ql_zi = q_l(zi, T_zi, qM);
    sft = sjump(u, p, LWP, p.fttype) + sM;
    Tft = (sft-g*zi)/Cp;
    qft = qjump(u, p, LWP, p.fttype) + qM;
    Δsv = (sft-sM) + Cp*(Rv/Rd-1)*(Tft*qft - T_zi*qM) + Cp*(Rv/Rd)*(T_zi*ql_zi);
    return Δsv
end

"""
    we(u, p, zb, LWP, etype::enBal)

    entrainment velocity obtained via energy balance requirement
    w = ΔR / (Δsv * ρref)
    Δsv = the jump in virtual liquid static energy
    sv = Cp*Tv + g*z - Lv*ql = s + Cp(Tv - T)
"""
function we(u, p, zb, LWP, etype::enBal)
    zi, sM, qM, SST, CF = u;
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
    zi, sM, qM, SST = u;

    # calculate sv flux <w'sv'>(z) = f0(z) + we*f1(z)
    # split into two terms I0, I1 where Ii = \integral fi(z) dz
    S0 = S_0(u, p, p.ftype);
    Q0 = Q_0(u, p, p.ftype);
    H0 = S0 + L0*Q0; # h = s + Lv*qt

    A0 = H0 - μ*L0*Q0;
    B0 = β*H0 - ϵ*L0*Q0;
    I0 = A0 * (zb - (zb^2)/(2*zi)) + B0 * ((zi-zb) + (zi^2 - zb^2)/(2*zi));
    
    sj = sjump(u, p, LWP, p.fttype);
    qj = qjump(u, p, LWP, p.fttype);
    hj = sj + L0*qj; # h = s + Lv*qt
    A1 = hj - μ*L0*qj;
    B1 = β*hj- ϵ*L0*qj;
    I1 = A1 * (-(zb^2)/(2*zi)) + B1 * ((-zi^2 + zb^2)/(2*zi));
    
    A = 2.0;
    α = (2.5 * A) / (zi * sv_jump(u, p, LWP));
    w = α*I0 / (1 - α*I1)
    return w
end
