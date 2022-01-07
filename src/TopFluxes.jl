export sstdep, co2dep, fixedFT, twocol
export hjump, qjump, H_zi, Q_zi

###########
# create structure for fttype
###########
abstract type ft_type end
struct sstdep <: ft_type end
struct co2dep <: ft_type end
struct fixedFT <: ft_type end
struct twocol <: ft_type end

"""
    qjump(u, p, p.fttype::sstdep)
    defines the inversion jump for qt 
        via linear regression to LES results
"""
function qjump(u, p, zb, fttype::sstdep)
    zi, hM, qM, SST, CF = u;
    qj = p.qj_m * (SST-290) + p.qj_b; # kg/kg
    return qj
end

"""
    hjump(u, p, p.fttype::sstdep)
    defines the inversion jump for h 
        via linear regression to LES results
"""
function hjump(u, p, zb, fttype::sstdep)
    zi, hM, qM, SST, CF = u;
    hj = p.hj_m * (SST-290) + p.hj_b; # m^2/s^2 = J/kg
    return hj
end

"""
    qjump(u, p, p.fttype::co2dep)
    defines the inversion jump for qt 
        via linear regression to LES results
"""
function qjump(u, p, zb, fttype::co2dep)
    zi, hM, qM, SST, CF = u;
    qj = -2.19e-6 * p.CO2 - 4.04e-3; # kg/kg
    qj /= min(1, CF*2.5)
    return qj
end

"""
    hjump(u, p, p.fttype::co2dep)
    defines the inversion jump for h 
        via linear regression to LES results
"""
function hjump(u, p, zb, fttype::co2dep)
    zi, hM, qM, SST, CF = u;
    hj = -6.07 * p.CO2 + 157.0; # m^2/s^2 = J/kg
    hj /= min(1, CF*1.5)
    return hj
end

"""
    qjump(u, p, p.fttype::fixedFT)
    defines qt+(z) in free troposphere -- given Gamma_q
    minimum value of qft of 2 g/kg
"""
function qjump(u, p, zb, fttype::fixedFT)
    zi, hM, qM, SST, CF = u;
    qft = p.qft0 + p.Gamma_q * zi;
    qft = max(qft, 2e-3);
    qj = qft - qM;
    return qj
end

"""
    hjump(u, p, p.fttype::fixedFT)
    defines h+(z) in free troposphere -- given Gamma_s and Gamma_q
"""
function hjump(u, p, zb, fttype::fixedFT)
    zi, hM, qM, SST, CF = u;
    sft = p.sft0 + p.Gamma_s * zi;
    qft = qjump(u, p, zb, p.fttype) + qM;
    hft = Cp * sft + L0 * qft;
    hj = hft - hM;
    return hj
end

"""
    temp_ft(u, p)

    free-tropospheric temperature given tropical sst
"""
function temp_ft(u, p, zb, zft)
    SST_trop = trop_sst(u, p, zb);
    Tft = SST_trop - zft*Γm(SST_trop);
    return Tft
end

"""
    qjump(u, p, p.fttype::twocol)

    specific humidity above cloud given fixed RH=0.2
    and saturation calculated at Tft and fixed 1500 m 
"""
function qjump(u, p, zb, fttype::twocol)
    zi, hM, qM, SST, CF = u;
    zft = 1500.0;
    qft = 0.2 * q_sat(zft, temp_ft(u, p, zb, zft));
    qj = qft - qM;
    return qj
end

"""
    hjump(u, p, p.fttype::twocol)
"""
function hjump(u, p, zb, fttype::twocol)
    zi, hM, qM, SST, CF = u;
    qft = qjump(u, p, zb, p.fttype) + qM;
    zft = zi;
    Tft = temp_ft(u, p, zb, zft);
    hft = Cp*Tft + g*zft + L0*qft;
    hj = hft - hM;
    return hj
end

"""
    H_zi(u, p, ent)

    moist enthalpy flux into the mixed-layer from above at z=zi
    H_zi = -we * (hft - hM)
"""
function H_zi(u, p, ent, zb)
    hj = hjump(u, p, zb, p.fttype);
    Hzi = -ent * hj;
    return Hzi
end

"""
    Q_zi(u, p, ent)

    moisture flux into the mixed-layer from above at z=zi
    Q_zi = -we * (qft - qM)
"""
function Q_zi(u, p, ent, zb)
    qj = qjump(u, p, zb, p.fttype);
    Qzi = -ent * qj;
    return Qzi
end