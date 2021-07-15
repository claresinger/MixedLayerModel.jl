export h_ft, q_ft, hjump, qjump, H_zi, Q_zi

function hjump(p)
    hj = -6.07 * p.CO2 + 157.0; # m^2/s^2
    return hj
end

function qjump(p)
    qj = -2.19e-6 * p.CO2 - 4.04e-3; # kg/kg
    return qj
end

"""
    H_zi(u, p)

    enthalpy flux into the mixed-layer from above at z=zi
    H_zi = -we * (hft - hM)
"""
function H_zi(u, p)
    zi, hM, qM, SST = u;
    #hj = h_ft(zi, p) - hM;
    hj = hjump(p);
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
    #qj = q_ft(zi, p) - qM;
    qj = qjump(p);
    Qzi = - we(u, p, p.etype) * qj;
    return Qzi
end

"""
    q_ft(z, p)

    defines qt+(z) in free troposphere -- given Gamma_q
"""
function q_ft(z, p)
    qft = p.qft0 .+ p.Gamma_q .* z
    qft = max.(qft, 2e-3)
    return qft
end

"""
    h_ft(z, p)

    defines h+(z) in free troposphere -- given Gamma_s and Gamma_q
"""
function h_ft(z, p)
    sft = p.sft0 .+ p.Gamma_s .* z
    hft = sft .* Cp .+ L0 .* q_ft(z, p)
    return hft
end