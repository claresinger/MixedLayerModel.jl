export calc_bflux
export calc_OHU
export calc_SWCRE

"""
    calc_bflux(u, p, zarr, etype::bflux)

    u is the state vector [zi, sM, qM, SST, CF]
    p is the parameter object
    zarr is an array of altitudes from 0 to some maxz > zi

    calculates the buoyancy flux for plotting
"""
function calc_bflux(u, p, zarr, etype::bflux)
    zi, sM, qM, SST, CF = u;
    zb = calc_LCL(u);
    LWP = incloud_LWP(u, zb);

    z1 = zarr[zarr .< zb];
    z2 = intersect(zarr[zarr .>= zb], zarr[zarr .< zi];)
    z3 = zarr[zarr .>= zi];
                                                                                        
    S0 = S_0(u, p, p.ftype);
    Q0 = Q_0(u, p, p.ftype);
    H0 = S0 + L0*Q0;

    ent = we(u, p, zb, LWP, p.etype);
    Szi = S_zi(u, p, ent, LWP);
    Qzi = Q_zi(u, p, ent, LWP);
    Hzi = Szi + L0*Qzi;

    wh(z) = (1 - z/zi) * H0 + (z/zi) * Hzi;
    wq(z) = (1 - z/zi) * Q0 + (z/zi) * Qzi;
    wsv_1(z) = wh(z) - (μ*L0) * wq(z);
    wsv_2(z) = β * wh(z) - (ϵ*L0) * wq(z);
    
    wsv_z = [wsv_1.(z1); wsv_2.(z2); zeros(length(z3))];
    Tsurf =  temp(0.0, sM, qM);
    bflux = wsv_z * (g / Cp / Tsurf);
    
    return bflux
end

"""
    calc_OHU(u, p, p.stype::fixSST)

    u is the state vector [zi, sM, qM, SST, CF]
    p is the parameter object

    calculates the ocean heat uptake as the residual 
    between the radiative fluxes and the LHF + SHF
"""
function calc_OHU(u, p, LWP, stype::fixSST)
    OHU = calc_surf_RAD(u, p, LWP) - calc_SHF(u,p) - calc_LHF(u,p);
    return OHU
end

"""
    calc_OHU(u, p, p.stype::varSST)

    u is the state vector [zi, sM, qM, SST, CF]
    p is the parameter object

    returns p.OHU
"""
function calc_OHU(u, p, LWP, stype::varSST)
    return p.OHU
end

function calc_SWCRE(u)
    zi, sM, qM, SST, CF = u;
    zb = calc_LCL(u);
    LWP = incloud_LWP(u, zb);
    αc = cloud_albedo(LWP);
    SWCRE = S_subtr * CF * (α_ocean - αc);
    return SWCRE
end
