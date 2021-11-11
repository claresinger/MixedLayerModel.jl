export sstdep, co2dep, fixedFT
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
function qjump(u, p, fttype::sstdep)
    zi, hM, qM, SST, CF = u;
    qj = p.qj_m * (SST-290) + p.qj_b; # kg/kg
    return qj
end

"""
    hjump(u, p, p.fttype::sstdep)
    defines the inversion jump for h 
        via linear regression to LES results
"""
function hjump(u, p, fttype::sstdep)
    zi, hM, qM, SST, CF = u;
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
    minimum value of qft of 2 g/kg
"""
function qjump(u, p, fttype::fixedFT)
    zi, hM, qM, SST, CF = u;
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
    zi, hM, qM, SST, CF = u;
    sft = p.sft0 .+ p.Gamma_s .* zi;
    qft = qjump(u, p, p.fttype) + qM;
    hft = Cp .* sft .+ L0 .* qft;
    hj = hft - hM;
    return hj
end

"""
    trop_sst(u,p)
"""
function trop_sst(u, p)
    zi, hM, qM, SST, CF = u;

    T = 0.8;
    LWP = incloud_LWP(u)*1e3; # kg/m^2 \to g/m^2
    αc = cloud_albedo(LWP, CF);
    α = T*CF*αc + (1-CF)*α_ocean;
    SW_net = (1-α)*S_subtrop;
    
    OLR = -491.0 - 2.57*SST;
    R_s = SW_net - OLR;
    ΔR_s = R_s - p.R_s_400;
    
    AreaFrac = 0.065;
    ΔR_t = -AreaFrac/(1-AreaFrac)*ΔR_s;
    
    sst_t = 301.0 - 2.4*ΔR_t;
    return sst_t
end

"""
    qjump(u, p, p.fttype::twocol)
"""
function qjump(u, p, fttype::twocol)
    zi, hM, qM, SST, CF = u;

    qj = qft - qM;
    return qj
end

"""
    hjump(u, p, p.fttype::twocol)
"""
function hjump(u, p, fttype::twocol)
    zi, hM, qM, SST, CF = u;
    
    qft = qjump(u, p, p.fttype) + qM;
    hft = Cp .* sft .+ L0 .* qft;
    hj = hft - hM;
    return hj
end

"""
    H_zi(u, p)

    enthalpy flux into the mixed-layer from above at z=zi
    H_zi = -we * (hft - hM)
"""
function H_zi(u, p)
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
    qj = qjump(u, p, p.fttype);
    Qzi = - we(u, p, p.etype) * qj;
    return Qzi
end
