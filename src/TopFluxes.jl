export sstdep, co2dep, fixEIS, EISco2, fixedFT, twocol
export sjump, qjump, S_zi, Q_zi

###########
# create structure for fttype
###########
abstract type ft_type end
# struct sstdep <: ft_type end
struct co2dep <: ft_type end
struct fixEIS <: ft_type end
# struct EISco2 <: ft_type end
# struct fixedFT <: ft_type end
struct twocol <: ft_type end

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
    # qj = max(qj, 2e-3 - qM);
    # println("qj: ", qj*1e3, " ", (qj+qM)*1e3, " ", qM*1e3)

    # [-0.00468956  0.00172823  0.01020646]
    # qj = -0.0047 - log(p.CO2/400)*0.0017 - (CFmax-CF)*0.01
    # qj = -0.0047 - log(p.CO2/400)*0.0017 - (CFmax-CF)*0.02
    
    # if qj + qM < 2e-3
    #     println("~~~~~~~~~~~~SMALL QFT~~~~~~~~~~~~~~~~")
    # end
    # qj = max(qj, 2e-3 - qM);

    # sft = sjump(u, p, LWP, p.fttype) + sM;
    # Tft = (sft - g*zi)/Cp;
    # qft = p.RHft * q_sat(zi, Tft);
    # qj = qft - qM;

    return qj
end

"""
    sjump(u, p, LWP, p.fttype::co2dep)
    defines the inversion jump for s
        via linear regression to LES results
"""
function sjump(u, p, LWP, fttype::co2dep)
    zi, sM, qM, SST, CF = u;
    # sj = (-6.07 + 2.5*2.19) * p.CO2 + (157.0 + 2.5e3*4.04); # m^2/s^2 = J/kg
    # sj /= min(1, CF*5)
    # sj = max(sj, Cp*SST + g*zi - sM)

    hj = -6.07 * p.CO2 + 157.0; # m^2/s^2 = J/kg
    hj /= min(1, CF*1.5);
    qj = qjump(u, p, LWP, p.fttype);
    sj = hj - L0*qj;
    # println("sj: ", sj/Cp, " ", hj/Cp, " ", qj*1e3)

    # 10051.39210848   404.81418053  3024.85246107
    # sj = 10050 - log(p.CO2/400)*405 - (CFmax-CF)*3025
    # sj = 10050 - log(p.CO2/400)*405 - (CFmax-CF)*3000
    return sj
end

# """
#     qjump(u, p, LWP, p.fttype::fixEIS)
#     defines qft from p.RHft and Tft(sft)
# """
# function qjump(u, p, LWP, fttype::fixEIS)
#     zi, sM, qM, SST, CF = u;
#     sft = sjump(u, p, LWP, p.fttype) + sM;
#     Tft = (sft - g*zi)/Cp;
#     qft = p.RHft * q_sat(zi, Tft);
#     qj = qft - qM;
#     return qj
# end

# """
#     sjump(u, p, LWP, p.fttype::fixEIS)
#     defines s+(z) in free troposphere given EIS and dTdz
# """
# function sjump(u, p, LWP, fttype::fixEIS)
#     zi, sM, qM, SST, CF = u;
#     Tft = p.SST0 + p.EIS + p.dTdz*zi;
#     sft = Cp*Tft + g*zi;
#     sj = sft - sM;
#     return sj
# end

# """
#     sjump(u, p, LWP, p.fttype::EISco2)
#     defines s+(z) in free troposphere 
# """
# function sjump(u, p, LWP, fttype::EISco2)
#     zi, sM, qM, SST, CF = u;
#     # EIS = (p.EIS - log(p.CO2/400)) * CF;
#     EIS = (p.EIS - (p.CO2/400)) * CF;
#     # EIS = (10.5 - 0.5*log(p.CO2/400)) - (0.5/CF);
#     Tft = p.SST0 + EIS + p.dTdz*zi;
#     sft = Cp*Tft + g*zi;
#     sj = sft - sM;
#     return sj
# end

# """
#     qjump(u, p, LWP, p.fttype::fixedFT)
#     defines qft from p.RHft and Tft(sft)
# """
# function qjump(u, p, LWP, fttype::fixedFT)
#     zi, sM, qM, SST, CF = u;
#     sft = sjump(u, p, LWP, p.fttype) + sM;
#     qft = p.RHft * q_sat(zi, (sft - g*zi)/Cp);
#     qj = qft - qM;
#     return qj
# end

# """
#     sjump(u, p, LWP, p.fttype::fixedFT)
#     defines s+(z) in free troposphere -- given Gamma_s and Gamma_q
# """
# function sjump(u, p, LWP, fttype::fixedFT)
#     zi, sM, qM, SST, CF = u;
#     sft = p.sft0 + p.Gamma_s * zi;
#     sj = sft - sM;
#     return sj
# end

# """
#     qjump(u, p, LWP, p.fttype::twocol)

#     specific humidity above cloud given fixed RH=0.2
#     and saturation calculated at Tft
# """
# function qjump(u, p, LWP, fttype::twocol)
#     zi, sM, qM, SST, CF = u;
#     sft = sjump(u, p, LWP, p.fttype) + sM;
#     Tft = (sft - g*zi)/Cp;
#     qft = p.RHft * q_sat(zi, Tft);
#     qj = qft - qM;
#     return qj
# end

# """
#     sjump(u, p, LWP, p.fttype::twocol)
# """
# function sjump(u, p, LWP, fttype::twocol)
#     zi, sM, qM, SST, CF = u;
#     SST_trop = trop_sst(u, p, LWP);
#     Tft = temp_ft(SST_trop, zi, p);
#     sft = Cp*Tft + g*zi;
#     sj = sft - sM;
#     return sj
# end

# """
#     qjump(u, p, LWP, type)

#     specific humidity above cloud given fixed RHft
#     and saturation calculated at Tft
# """
# function qjump(u, p, LWP, type)
#     zi, sM, qM, SST, CF = u;
#     sft = sjump(u, p, LWP, p.fttype) + sM;
#     Tft = (sft - g*zi)/Cp;
#     qft = p.RHft * q_sat(zi, Tft);
#     qj = qft - qM;
#     return qj
# end

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