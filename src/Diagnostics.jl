export calc_bflux

"""
    calc_bflux(u, p, zarr, etype::bflux)

    u is the state vector [zi, hM, qM, SST]
    p is the parameter object
    zarr is an array of altitudes from 0 to some maxz > zi

    calculates the buoyancy flux for plotting
"""
function calc_bflux(u, p, zarr, etype::bflux)
    zi,hM,qM,SST = u;
    zb = calc_LCL(zi,hM,qM);

    z1 = zarr[zarr .< zb];
    z2 = intersect(zarr[zarr .>= zb], zarr[zarr .< zi];)
    z3 = zarr[zarr .>= zi];
                                                                                        
    H0 = H_0(u, p, p.ftype);
    Q0 = Q_0(u, p, p.ftype);

    Hzi = H_zi(u, p);
    Qzi = Q_zi(u, p);

    wh(z) = (1 .- z./zi) .* H0 .+ (z./zi) .* Hzi;
    wq(z) = (1 .- z./zi) .* Q0 .+ (z./zi) .* Qzi;
    wsv_1(z) = wh(z) .- (μ*L0) .* wq(z);
    wsv_2(z) = β .* wh(z) .- (ϵ*L0) .* wq(z);
    
    wsv_z = [wsv_1(z1); wsv_2(z2); zeros(length(z3))];
    Tsurf =  temp(0.0, hM, qM);
    bflux = wsv_z * (g / Cp / Tsurf);
    
    return bflux
end