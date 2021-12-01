export calc_S, cloud_fraction

"""
    calc_S(u, p)

    calculates the stability parameter, s
    also called the minimal decloupling parameter
    in Bretherton and Wyant (1997)
    S = (LHF / ΔR) * (zc / zi)
    where zc is the cloud thickness (zi-zb)
"""
function calc_S(u, p, zb)
    zi, hM, qM, SST, CF = u;
    LHF = calc_LHF(u, p);
    ΔR = calc_cloudtop_RAD(u, p, zb, p.rtype);
    zc = zi - zb;
    S = (LHF/ΔR)*(zc/zi);
    return S
end

"""
    cloud_fraction(u, p)

    calculates the cloud fraction 
    inspired by Chung et al 2012
    uses a logistic function to smoothly interpolate between 
    100% and 20% cloud fraction based on the value of the
    stability parameter
"""
function cloud_fraction(u, p, zb)
    S = calc_S(u, p, zb);
    m = 10; # tunable parameter for the slope of the CF nonlinearity
    S_crit = 0.7; # tunable parameter for the halfway point of CF decrease
    CF = 1 - 0.8 / (1 + exp(-m*(S-S_crit)));
    
    # S_min = 1/(2*σ); # minimum decoupling parameter, ≈0.55
    # if S < S_min
    #     CF = 1.0
    # end
    return CF
end