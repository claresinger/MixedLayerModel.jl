module TopFluxes

include("Definitions.jl")
using ..Entrainment

export H_zi, Q_zi

# define H(zi) function
function H_zi(u, p)
    zi, hM, qM, SST = u;
    hft = h_ft(zi, p);
    Hzi = -we(u, p, p.etype) * (hft - hM);
    return Hzi
end

# define Q(zi) function
function Q_zi(u, p)
    zi, hM, qM, SST = u;
    qft = q_ft(zi, p);
    Qzi = - we(u, p, p.etype) * (qft - qM);
    return Qzi
end

end