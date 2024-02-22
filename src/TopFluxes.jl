export fixEIS, fixedFT, co2EIS
export sjump, qjump, S_zi, Q_zi

###########
# create structure for fttype
###########
abstract type ft_type end
struct fixedFT <: ft_type end
struct fixEIS <: ft_type end
struct co2EIS <: ft_type end

# adjustment factor based on degree of decoupling (or CF)
# if CF = p.CFmax then λt = 1, if CF = p.CFmin then λt = 1 - p.λtop
function λt(u,p)
    zi, sM, qM, SST, CF = u;
    return 1 - p.λtop * ((p.CFmax - CF) / (p.CFmax - p.CFmin));
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
    sjump(u, p, LWP, p.fttype::fixEIS)
    defines s+(z) in free troposphere given fixed EIS
"""
function sjump(u, p, LWP, fttype::fixEIS)
    zi, sM, qM, SST, CF = u;
    Tft = temp(zi, sM, qM) + p.EIS0;
    # Tft = SST + p.EIS0 + p.dTdz*zi;
    sft = Cp*Tft + g*zi;
    sj = sft - sM;
    return sj
end

"""
    sjump(u, p, LWP, p.fttype::co2EIS)
    defines s+(z) in free troposphere given EIS = f(CO2)
"""
function sjump(u, p, LWP, fttype::co2EIS)
    zi, sM, qM, SST, CF = u;
    EIS = p.EIS0 + (p.ECS/log(2))*log(p.CO2 / 400) - p.Eexport*(p.CFmax - CF);
    Tft = temp(zi, sM, qM) + EIS;
    # Tft = SST + EIS + p.dTdz*zi;
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
    return -ent * λt(u,p) * sj
end

"""
    Q_zi(u, p, ent, LWP)

    moisture flux into the mixed-layer from above at z=zi
    Q_zi = -we * (qft - qM)
"""
function Q_zi(u, p, ent, LWP)
    qj = qjump(u, p, LWP, p.fttype);
    return -ent * λt(u,p) * qj
end