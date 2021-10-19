export calc_S, cloud_fraction

"""
    calc_S(u, p)

    calculates the stability parameter, s
    S = (LHF / ΔR) * (zc / zi)
    where zc is the cloud thickness (zi-zb)
"""
function calc_S(u, p)
    zi, hM, qM, SST, CF = u;
    LHF = calc_LHF(u, p);
    ΔR = calc_cloudtop_RAD(u, p, p.rtype);
    zb = calc_LCL(u);
    zc = zi - zb;
    S = (LHF/ΔR)*(zc/zi);
    return S
end

"""
    cloud_fraction(u, p)

    calculates the cloud fraction 
    inspired by Chung et al 2012
    uses a logistic function to smoothly interpolate between 
    100% and 10% cloud fraction based on the value of the
    stability parameter
"""
function cloud_fraction(u, p)
    S = calc_S(u, p);
    m = 10; # tunable parameter for the slope of the CF nonlinearity
    S_crit = 0.75; # tunable parameter for the halfway point of CF decrease
    CF = 1 - 0.9 / (1 + exp(-m*(S-S_crit)));
    # S_min = 0.5; # minimum decoupling parameter
    # if S < S_min
    #     CF = 1.0
    # end
    return CF
end