export calc_S, cloud_fraction, cloud_fraction_S

"""
    calc_S(u, p)

    calculates the stability parameter, s
    also called the minimal decloupling parameter
    in Bretherton and Wyant (1997)
    S = (LHF / ΔR) * (zc / zi)
    where zc is the cloud thickness (zi-zb)
"""
function calc_S(u, p, zb, LWP)
    zi, hM, qM, SST, CF = u;
    LHF = calc_LHF(u, p);
    ΔR = calc_cloudtop_RAD(u, p, LWP, p.rtype);
    zc = zi - zb;
    S = (LHF/ΔR)*(zc/zi);
    return S
end

"""
    cloud_fraction(u, p, zb, LWP)

    calculates the cloud fraction 
    inspired by Chung et al 2012
    uses a logistic function to smoothly interpolate between 
    100% and 20% cloud fraction based on the value of the
    stability parameter
"""
function cloud_fraction(u, p, zb, LWP)
    S = calc_S(u, p, zb, LWP);
    CF = cloud_fraction_S(S);
    return CF
end

function cloud_fraction_S(S)
    m = 6; # tunable parameter for the slope of the CF nonlinearity
    S_crit = 0.55 + log(7)/m; # 50% value such that the 90% value is at 0.55
    CF = 1 - 0.8 / (1 + exp(-m*(S-S_crit)));
    return CF
end