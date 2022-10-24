export calc_decoupling, cloud_fraction, cloud_fraction_param

"""
    calc_decoupling(u, p)

    calculates the stability parameter, s
    also called the minimal decloupling parameter
    in Bretherton and Wyant (1997)
    De = (LHF / ΔR) * (zc / zi)
    where zc is the cloud thickness (zi-zb)
"""
function calc_decoupling(u, p, zb, LWP)
    zi, sM, qM, SST, CF = u;
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
    decoup = calc_decoupling(u, p, zb, LWP);
    CF = cloud_fraction_param(decoup, p);
    return CF
end

function cloud_fraction_param(decoup, p)
    m = p.decoup_slope; # tunable parameter for the slope of the CF nonlinearity
    dcrit = 1;
    CF = p.CFmax - (p.CFmax - p.CFmin) / (1 + (1/9)*exp(-m*(decoup-dcrit)));
    return CF
end