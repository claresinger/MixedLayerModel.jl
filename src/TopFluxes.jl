export co2dep, fixEIS, fixedFT, co2EIS, twocol
export sjump, qjump, S_zi, Q_zi

###########
# create structure for fttype
###########
abstract type ft_type end
struct co2dep <: ft_type end
struct fixedFT <: ft_type end
struct fixEIS <: ft_type end
struct co2EIS <: ft_type end
struct twocol <: ft_type end

"""
    qjump(u, p, LWP, p.fttype::co2dep)
    defines the inversion jump for qt 
        via linear regression to LES results
"""
function qjump(u, p, LWP, fttype::co2dep)
    zi, sM, qM, SST, CF = u;
    qj = -2.19e-6 * p.CO2 - 4.04e-3; # kg/kg
    qj /= min(1, CF*2.5);
    if qj + qM < 2e-3
        println("~~~~~~~~~~~~SMALL QFT~~~~~~~~~~~~~~~~")
    end
    qj = max(qj, 2e-3 - qM);

    # [-0.00468956  0.00172823  0.01020646]
    # qj = -0.0047 - log(p.CO2/400)*0.0017 - (p.CFmax-CF)*0.01
    return qj
end

"""
    sjump(u, p, LWP, p.fttype::co2dep)
    defines the inversion jump for s
        via linear regression to LES results
"""
function sjump(u, p, LWP, fttype::co2dep)
    zi, sM, qM, SST, CF = u;
    hj = -6.07 * p.CO2 + 157.0; # m^2/s^2 = J/kg
    hj /= min(1, CF*1.5);
    qj = qjump(u, p, LWP, p.fttype);
    sj = hj - L0*qj;

    # 10051.39210848   404.81418053  3024.85246107
    # sj = 10050 - log(p.CO2/400)*405 - (p.CFmax-CF)*3025
    return sj
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
    EIS = p.EIS0 + p.ECS*log(p.CO2 / 400) - p.Eexport*(p.CFmax - CF);
    Tft = SST + p.dTdz*zi + EIS;
    sft = Cp*Tft + g*zi;
    sj = sft - sM;
    return sj
end

"""
    sjump(u, p, LWP, p.fttype::twocol)
"""
function sjump(u, p, LWP, fttype::twocol)
    zi, sM, qM, SST, CF = u;
    SST_trop = trop_sst(u, p, LWP);
    Tft = temp_ft(SST_trop, zi, p);
    sft = Cp*Tft + g*zi;
    sj = sft - sM;
    return sj
end

"""
    qjump(u, p, LWP, fttype::Union{twocol, fixEIS, fixedFT, co2EIS})

    specific humidity above cloud given fixed RHft
    and saturation calculated at Tft
"""
function qjump(u, p, LWP, fttype::Union{twocol, fixEIS, fixedFT, co2EIS})
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