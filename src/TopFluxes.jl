export sstdep, co2dep, fixedFT, twocol
export hjump, qjump, H_zi, Q_zi

###########
# create structure for fttype
###########
abstract type ft_type end
# struct sstdep <: ft_type end
struct co2dep <: ft_type end
struct fixedFT <: ft_type end
struct twocol <: ft_type end

zft = 1500;

# """
#     qjump(u, p, LWP, p.fttype::sstdep)
#     defines the inversion jump for qt 
#         via linear regression to LES results
# """
# function qjump(u, p, LWP, fttype::sstdep)
#     zi, hM, qM, SST, CF = u;
#     qj = p.qj_m * (SST-290) + p.qj_b; # kg/kg
#     return qj
# end

# """
#     hjump(u, p, LWP, p.fttype::sstdep)
#     defines the inversion jump for h 
#         via linear regression to LES results
# """
# function hjump(u, p, LWP, fttype::sstdep)
#     zi, hM, qM, SST, CF = u;
#     hj = p.hj_m * (SST-290) + p.hj_b; # m^2/s^2 = J/kg
#     return hj
# end

"""
    qjump(u, p, LWP, p.fttype::co2dep)
    defines the inversion jump for qt 
        via linear regression to LES results
"""
function qjump(u, p, LWP, fttype::co2dep)
    zi, hM, qM, SST, CF = u;
    qj = -2.19e-6 * p.CO2 - 4.04e-3; # kg/kg
    qj /= min(1, CF*2.5);
    qj = max(qj, 2e-3 - qM);
    return qj
end

"""
    hjump(u, p, LWP, p.fttype::co2dep)
    defines the inversion jump for h 
        via linear regression to LES results
"""
function hjump(u, p, LWP, fttype::co2dep)
    zi, hM, qM, SST, CF = u;
    hj = -6.07 * p.CO2 + 157.0; # m^2/s^2 = J/kg
    hj /= min(1, CF*1.5);
    hj = max(hj, Cp*SST+g*zi+L0*(qM+qjump(u,p,LWP,p.fttype)) - hM);
    return hj
end

"""
    qjump(u, p, LWP, p.fttype::fixedFT)
    defines qt+(z) in free troposphere -- given Gamma_q
    minimum value of qft of 2 g/kg
"""
function qjump(u, p, LWP, fttype::fixedFT)
    zi, hM, qM, SST, CF = u;
    qft = p.qft0 + p.Gamma_q * zi;
    # qft = max(qft, 2e-3);
    qj = qft - qM;
    return qj
end

"""
    hjump(u, p, LWP, p.fttype::fixedFT)
    defines h+(z) in free troposphere -- given Gamma_s and Gamma_q
"""
function hjump(u, p, LWP, fttype::fixedFT)
    zi, hM, qM, SST, CF = u;
    sft = p.sft0 + p.Gamma_s * zi;
    qft = qjump(u, p, LWP, p.fttype) + qM;
    hft = sft + L0 * qft;
    hj = hft - hM;
    return hj
end

# """
#     qjump(u, p, LWP, p.fttype::fixedFT)
#     defines qt+(z) in free troposphere -- given Gamma_q
#     minimum value of qft of 2 g/kg
# """
# function qjump(u, p, LWP, fttype::fixedFT)
#     zi, hM, qM, SST, CF = u;
#     Tft = temp(zi, hM, qM) + p.ΔTft
#     qft = p.RHft * q_sat(zi, Tft);
#     qj = qft - qM;
#     return qj
# end

# """
#     hjump(u, p, LWP, p.fttype::fixedFT)
#     defines h+(z) in free troposphere -- given Gamma_s and Gamma_q
# """
# function hjump(u, p, LWP, fttype::fixedFT)
#     zi, hM, qM, SST, CF = u;
#     Tzi = temp(zi, hM, qM);
#     hj = Cp*p.ΔTft + L0*(qjump(u, p, LWP, p.fttype) - q_l(zi, Tzi, qM));
#     return hj
# end

"""
    qjump(u, p, LWP, p.fttype::twocol)

    specific humidity above cloud given fixed RH=0.2
    and saturation calculated at Tft
"""
function qjump(u, p, LWP, fttype::twocol)
    zi, hM, qM, SST, CF = u;
    Tft = temp_ft(trop_sst(u, p, LWP), zft, p);
    qft = p.RHtropft * q_sat(zft, Tft);
    qj = qft - qM;
    return qj
end

"""
    hjump(u, p, LWP, p.fttype::twocol)
"""
function hjump(u, p, LWP, fttype::twocol)
    zi, hM, qM, SST, CF = u;
    qft = qjump(u, p, LWP, p.fttype) + qM;
    Tft = temp_ft(trop_sst(u, p, LWP), zft, p);
    hft = Cp*Tft + g*zft + L0*qft;
    hj = hft - hM;
    return hj
end

"""
    H_zi(u, p, ent, LWP)

    moist enthalpy flux into the mixed-layer from above at z=zi
    H_zi = -we * (hft - hM)
"""
function H_zi(u, p, ent, LWP)
    hj = hjump(u, p, LWP, p.fttype);
    Hzi = -ent * hj;
    return Hzi
end

"""
    Q_zi(u, p, ent, LWP)

    moisture flux into the mixed-layer from above at z=zi
    Q_zi = -we * (qft - qM)
"""
function Q_zi(u, p, ent, LWP)
    qj = qjump(u, p, LWP, p.fttype);
    Qzi = -ent * qj;
    return Qzi
end