export h_ft, q_ft, H_zi, Q_zi

"""
    H_zi(u, p)

    enthalpy flux into the mixed-layer from above at z=zi
    H_zi = -we * (hft - hM)
"""
function H_zi(u, p)
    zi, hM, qM, SST = u;
    hft = h_ft(zi, p);
    Hzi = -we(u, p, p.etype) * (hft - hM);
    return Hzi
end

"""
    Q_zi(u, p)

    moisture flux into the mixed-layer from above at z=zi
    Q_zi = -we * (qft - qM)
"""
function Q_zi(u, p)
    zi, hM, qM, SST = u;
    qft = q_ft(zi, p);
    Qzi = - we(u, p, p.etype) * (qft - qM);
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