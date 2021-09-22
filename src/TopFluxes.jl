export sstdep, co2dep, fixedFT
export hjump, qjump, H_zi, Q_zi

###########
# create structure for fttype
###########
abstract type ft_type end
struct sstdep <: ft_type end
struct co2dep <: ft_type end
struct fixedFT <: ft_type end

"""
    qjump(u, p, p.fttype::sstdep)
    defines the inversion jump for qt 
        via linear regression to LES results
"""
function qjump(u, p, fttype::sstdep)
    zi, hM, qM, SST = u;
    qj = p.qj_m * (SST-290) + p.qj_b; # kg/kg
    return qj
end

"""
    hjump(u, p, p.fttype::sstdep)
    defines the inversion jump for h 
        via linear regression to LES results
"""
function hjump(u, p, fttype::sstdep)
    zi, hM, qM, SST = u;
    hj = p.hj_m * (SST-290) + p.hj_b; # m^2/s^2 = J/kg
    return hj
end

"""
    qjump(u, p, p.fttype::co2dep)
    defines the inversion jump for qt 
        via linear regression to LES results
"""
function qjump(u, p, fttype::co2dep)
    qj = -2.19e-6 * p.CO2 - 4.04e-3; # kg/kg
    return qj
end

"""
    hjump(u, p, p.fttype::co2dep)
    defines the inversion jump for h 
        via linear regression to LES results
"""
function hjump(u, p, fttype::co2dep)
    hj = -6.07 * p.CO2 + 157.0; # m^2/s^2 = J/kg
    return hj
end

"""
    qjump(u, p, p.fttype::fixedFT)
    defines qt+(z) in free troposphere -- given Gamma_q
"""
function qjump(u, p, fttype::fixedFT)
    zi, hM, qM, SST = u;
    qft = p.qft0 .+ p.Gamma_q .* zi;
    qft = max.(qft, 2e-3);
    qj = qft - qM;
    return qj
end

"""
    hjump(u, p, p.fttype::fixedFT)
    defines h+(z) in free troposphere -- given Gamma_s and Gamma_q
"""
function hjump(u, p, fttype::fixedFT)
    zi, hM, qM, SST = u;
    sft = p.sft0 .+ p.Gamma_s .* zi;
    hft = sft .* Cp .+ L0 .* q_ft(u, p, p.fttype);
    hj = hft - hM;
    return hj
end

"""
    H_zi(u, p)

    enthalpy flux into the mixed-layer from above at z=zi
    H_zi = -we * (hft - hM)
"""
function H_zi(u, p)
    zi, hM, qM, SST = u;
    hj = hjump(u, p, p.fttype);
    Hzi = -we(u, p, p.etype) * hj;
    return Hzi
end

"""
    Q_zi(u, p)

    moisture flux into the mixed-layer from above at z=zi
    Q_zi = -we * (qft - qM)
"""
function Q_zi(u, p)
    zi, hM, qM, SST = u;
    qj = qjump(u, p, p.fttype);
    Qzi = - we(u, p, p.etype) * qj;
    return Qzi
end