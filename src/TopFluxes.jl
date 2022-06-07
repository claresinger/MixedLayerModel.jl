export sstdep, co2dep, fixedFT, twocol
export sjump, qjump, S_zi, Q_zi

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
#     zi, sM, qM, SST, CF = u;
#     qj = p.qj_m * (SST-p.SST0) + p.qj_b; # kg/kg
#     return qj
# end

# """
#     sjump(u, p, LWP, p.fttype::sstdep)
#     defines the inversion jump for s
#         via linear regression to LES results
# """
# function sjump(u, p, LWP, fttype::sstdep)
#     zi, sM, qM, SST, CF = u;
#     sj = p.sj_m * (SST-p.SST0) + p.sj_b; # m^2/s^2 = J/kg
#     return sj
# end

"""
    qjump(u, p, LWP, p.fttype::co2dep)
    defines the inversion jump for qt 
        via linear regression to LES results
"""
function qjump(u, p, LWP, fttype::co2dep)
    zi, sM, qM, SST, CF = u;
    qj = -2.19e-6 * p.CO2 - 4.04e-3; # kg/kg
    qj /= min(1, CF*2.5);
    qj = max(qj, 2e-3 - qM);
    return qj
end

"""
    sjump(u, p, LWP, p.fttype::co2dep)
    defines the inversion jump for s
        via linear regression to LES results
"""
function sjump(u, p, LWP, fttype::co2dep)
    zi, sM, qM, SST, CF = u;
    sj = (-6.07 + 2.5*2.19) * p.CO2 + (157.0 + 2.5e3*4.04); # m^2/s^2 = J/kg
    sj /= min(1, CF*1.5)
    return sj
end

"""
    qjump(u, p, LWP, p.fttype::fixedFT)
    defines qt+(z) in free troposphere -- given Gamma_q
    minimum value of qft of 2 g/kg
"""
function qjump(u, p, LWP, fttype::fixedFT)
    zi, sM, qM, SST, CF = u;
    qft = p.qft0 + p.Gamma_q * zi;
    qft = max(qft, 2e-3);
    qj = qft - qM;
    return qj
end

"""
    sjump(u, p, LWP, p.fttype::fixedFT)
    defines s+(z) in free troposphere -- given Gamma_s and Gamma_q
"""
function sjump(u, p, LWP, fttype::fixedFT)
    zi, sM, qM, SST, CF = u;
    sft = p.sft0 + p.Gamma_s * zi;
    sj = sft - sM;
    return sj
end

"""
    qjump(u, p, LWP, p.fttype::twocol)

    specific humidity above cloud given fixed RH=0.2
    and saturation calculated at Tft and fixed 1500 m 
"""
function qjump(u, p, LWP, fttype::twocol)
    zi, sM, qM, SST, CF = u;
    SST_trop = trop_sst(u, p, LWP);
    Tft = temp_ft(SST_trop, zft, p);
    qft = p.RHtropft * q_sat(zft, Tft);
    qj = qft - qM;
    return qj
end

"""
    sjump(u, p, LWP, p.fttype::twocol)
"""
function sjump(u, p, LWP, fttype::twocol)
    zi, sM, qM, SST, CF = u;
    SST_trop = trop_sst(u, p, LWP);
    Tft = temp_ft(SST_trop, zft, p);
    sft = Cp*Tft + g*zft;
    sj = sft - sM;
    return sj
end

"""
    S_zi(u, p, ent, LWP)

    moist enthalpy flux into the mixed-layer from above at z=zi
    S_zi = -we * (sft - sM)
"""
function S_zi(u, p, ent, LWP)
    sj = sjump(u, p, LWP, p.fttype);
    return -ent * sj
end

"""
    Q_zi(u, p, ent, LWP)

    moisture flux into the mixed-layer from above at z=zi
    Q_zi = -we * (qft - qM)
"""
function Q_zi(u, p, ent, LWP)
    qj = qjump(u, p, LWP, p.fttype);
    return -ent * qj
end