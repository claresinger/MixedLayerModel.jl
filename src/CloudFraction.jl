export calc_decoupling, cloud_fraction, cloud_fraction_param

"""
    calc_decoupling(u, p, zb, LWP)

    calculates the decoupling parameter
    De = (σ LHF / ΔR)
    slightly different from Bretherton and Wyant (1997)
"""
function calc_decoupling(u, p, zb, LWP)
    LHF = calc_LHF(u, p);
    ΔR = calc_cloudtop_RAD(u, p, LWP, p.rtype);
    return σ * LHF / ΔR
end

"""
    cloud_fraction_param(De, p)

    calculates the cloud fraction as a function of decoupling parameter
    inspired by Chung et al 2012
    uses a logistic function to smoothly interpolate between 
    CFmax and CFmin cloud fraction based on the value of the
    stability parameter
"""
function cloud_fraction_param(De, p)
    return p.CFmax - (p.CFmax - p.CFmin) / (1 + (1/9)*exp(-p.Dslope*(De-p.Dcrit)))
end

"""
    cloud_fraction(u, p, zb, LWP)

    calculates the cloud fraction 
    as function of current state, u
"""
function cloud_fraction(u, p, zb, LWP)
    return cloud_fraction_param(calc_decoupling(u, p, zb, LWP), p);
end