export calc_decoupling, cloud_fraction, cloud_fraction_param

"""
    calc_decoupling(u, p)

    calculates the decoupling parameter
    De = (σ LHF / ΔR)
    slightly different from Bretherton and Wyant (1997)
"""
function calc_decoupling(u, p, zb, LWP)
    zi, sM, qM, SST, CF = u;
    LHF = calc_LHF(u, p);
    ΔR = calc_cloudtop_RAD(u, p, LWP, p.rtype);
    De = σ * LHF / ΔR

    # χ = randn(1)[1] + 1
    # χ = 2^randn(1)[1]
    # χ = 2^min(max(randn(1)[1], -3), 3)
    # χ = 2^min(max(randn(1)[1], -2), 2)
    # if CF < (p.CFmax + p.CFmin)/2
    #     De = De*χ
    # end

    # χ = 2^min(max(randn(1)[1], -3), 0)
    
    # χ = 2^min(randn(1)[1]*1.5, 0)
    # De = De*χ
    
    # if De > p.Dcrit
    #     De = De*χ
    # end


    # if CF < p.CFmin*1.1
    #     De -= 0.1*abs(randn(1)[1])
    # end
    return De
end

"""
    cloud_fraction(u, p, zb, LWP)

    calculates the cloud fraction 
    inspired by Chung et al 2012
    uses a logistic function to smoothly interpolate between 
    CFmax and CFmin cloud fraction based on the value of the
    stability parameter
"""
function cloud_fraction(u, p, zb, LWP)
    # zi, sM, qM, SST, CF = u;
    # De = calc_decoupling(u, p, zb, LWP)
    # χ = randn(1)[1];
    # if De > p.Dcrit
    #     p.CFmin = max(min(p.CFmin+0.1*χ, p.CFmax-0.1), 0.2)
    #     # if χ > 1
    #     #     p.CFmin = min(p.CFmin+0.1*χ, 0.5)
    #     # elseif χ < -1
    #     #     p.CFmin = max(p.CFmin+0.1*χ, 0.2);
    #     # end
    # else
    #     p.CFmin = 0.2
    # end
    CF = cloud_fraction_param(calc_decoupling(u, p, zb, LWP), p);
    # if CF < p.CFmin*1.1
    #     CF *= (1 + abs(rand(1)[1]));
    # end
    return CF
end

function cloud_fraction_param(De, p)
    return p.CFmax - (p.CFmax - p.CFmin) / (1 + (1/9)*exp(-p.Dslope*(De-p.Dcrit)))
end