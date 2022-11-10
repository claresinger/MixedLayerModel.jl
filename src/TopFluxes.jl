export fixEIS, fixedFT, co2EIS
export sjump, qjump, S_zi, Q_zi

###########
# create structure for fttype
###########
abstract type ft_type end
struct fixedFT <: ft_type end
struct fixEIS <: ft_type end
struct co2EIS <: ft_type end

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
    sjump(u, p, LWP, p.fttype::fixEIS)
    defines s+(z) in free troposphere given EIS and dTdz
"""
function sjump(u, p, LWP, fttype::fixEIS)
    zi, sM, qM, SST, CF = u;
    Tft = SST + p.EIS0 + p.dTdz*zi;
    sft = Cp*Tft + g*zi;
    sj = sft - sM;
    return sj
end

"""
    sjump(u, p, LWP, p.fttype::co2EIS)
    defines s+(z) in free troposphere given EIS and dTdz
"""
function sjump(u, p, LWP, fttype::co2EIS)
    zi, sM, qM, SST, CF = u;
    EIS = p.EIS0 + (p.ECS/log(2))*log(p.CO2 / 400) - p.Eexport*(p.CFmax - CF);
    Tft = SST + EIS + p.dTdz*zi;
    sft = Cp*Tft + g*zi;
    sj = sft - sM;
    return sj
end

"""
    qjump(u, p, LWP, fttype::Union{twocol, fixEIS, fixedFT, co2EIS})

    specific humidity above cloud given fixed RHft
    and saturation calculated at Tft
"""
function qjump(u, p, LWP, fttype::Union{fixedFT, fixEIS, co2EIS})
    zi, sM, qM, SST, CF = u;
    sft = sjump(u, p, LWP, p.fttype) + sM;
    Tft = (sft - g*zi)/Cp;
    qft = p.RHft * q_sat(zi, Tft);
    qj = qft - qM;
    return qj
end

"""
    S_zi(u, p, ent, LWP)

    energy flux into the mixed-layer from above at z=zi
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