using MixedLayerModel
using MixedLayerModel: g, Cp, T0, L0, Rv, Rd
using FileIO
using Plots
using Statistics

function linear_flux_profile(u, p)
    zi, sM, qM, SST, CF = u;
    zb = calc_LCL(u);
    LWP = incloud_LWP(u, zb);
                                                                                        
    S0 = S_0(u, p, p.ftype);
    Q0 = Q_0(u, p, p.ftype);

    ent = we(u, p, zb, LWP, p.etype);
    Szi = S_zi(u, p, ent, LWP);
    ΔR = calc_cloudtop_RAD(u, p, LWP, p.rtype)/ρref(SST);
    Qzi = Q_zi(u, p, ent, LWP);
    # println(S0, "\t", Szi, "\t", ΔR, "\t", Szi + ΔR)
    # println(Q0, "\t", Qzi)

    ϵ = Rd/Rv
    ϵ_ = (Rv-Rd)/Rd
    ΔT = 0.1
    T(z) = temp(z, sM, qM)
    γ(z) = q_sat(z, T(z)+ΔT) - q_sat(z, T(z)) / ΔT
    βsm(z) = (1 + q_sat(z, T(z))/ϵ - qM + (T0/ϵ)*γ(z)) / (1 + L0*γ(z)/Cp)
    α(z) = (1 + ϵ_*qM)*(z < zb) + βsm(z)*(zb <= z < zi)
    β(z) = ϵ_*(z < zb) + βsm(z)*(L0/Cp/T0)*(zb <= z < zi)

    S(z) = (1 - z/zi)*S0 + (z/zi)*(Szi + ΔR)
    Q(z) = (1 - z/zi)*Q0 + (z/zi)*Qzi

    B(z) = α(z)*(g/(Cp*T0))*S(z) + β(z)*g*Q(z)
    
    zarr = 0:1:zi+100
    normz = zarr ./ zi

    # println(minimum(α.(zarr)), "\t", mean(α.(zarr)), "\t", maximum(α.(zarr)))
    # println(minimum(β.(zarr)), "\t", mean(β.(zarr)), "\t", maximum(β.(zarr)))

    return normz, S.(zarr), Q.(zarr), B.(zarr)
end

function const_flux_profile(u, p)
    zi, sM, qM, SST, CF = u;
    zb = calc_LCL(u);
    LWP = incloud_LWP(u, zb);
                                                                                        
    S0 = S_0(u, p, p.ftype);
    Q0 = Q_0(u, p, p.ftype);

    ent = we(u, p, zb, LWP, p.etype);
    Szi = S_zi(u, p, ent, LWP);
    ΔR = calc_cloudtop_RAD(u, p, LWP, p.rtype)/ρref(SST);
    Qzi = Q_zi(u, p, ent, LWP);
    # println(S0, "\t", Szi, "\t", ΔR, "\t", Szi + ΔR)
    # println(Q0, "\t", Qzi)

    ϵ = Rd/Rv
    ϵ_ = (Rv-Rd)/Rd
    ΔT = 0.1
    T(z) = temp(z, sM, qM)
    γ(z) = q_sat(z, T(z)+ΔT) - q_sat(z, T(z)) / ΔT
    βsm(z) = (1 + q_sat(z, T(z))/ϵ - qM + (T0/ϵ)*γ(z)) / (1 + L0*γ(z)/Cp)
    α(z) = (1 + ϵ_*qM)*(z < zb) + βsm(z)*(zb <= z < zi)
    β(z) = ϵ_*(z < zb) + βsm(z)*(L0/Cp/T0)*(zb <= z < zi)

    S(z) = S0*(z < zb) + (Szi + ΔR)*(zb <= z < zi)
    Q(z) = Q0*(z < zb) + Qzi*(zb <= z < zi)

    B(z) = α(z)*(g/(Cp*T0))*S(z) + β(z)*g*Q(z)
    
    zarr = 0:1:zi+100
    normz = zarr ./ zi

    # println(minimum(α.(zarr)), "\t", mean(α.(zarr)), "\t", maximum(α.(zarr)))
    # println(minimum(β.(zarr)), "\t", mean(β.(zarr)), "\t", maximum(β.(zarr)))

    return normz, S.(zarr), Q.(zarr), B.(zarr)
end

path = "experiments/output/20230207/"
b = plot(0,0, xlabel="<w'b'> [m²/s³]", ylabel="z [m]", dpi=300)
sl = plot(0,0, xlabel="g/CₚT <w'sₗ'> [m/s]", ylabel="z [m]", dpi=300)
q = plot(0,0, xlabel="g <w'qₜ'> [m/s]", ylabel="z [m]", dpi=300)
for (i,co2) in enumerate([200, 300, 400, 600, 800, 1000, 1200, 1300, 1400, 1600])
    file = path*"co2_upstep_"*string(co2)*".jld2"
    dat = load(file)
    uf = dat["uf"];
    par = dat["p"];
    zarr, sarr, qarr, barr = const_flux_profile(uf, par)
    plot!(b, barr * 1e5, zarr, label=string(co2), palette=palette([:blue, :red], 11), legend=:bottomleft)
    plot!(sl, sarr*g/(Cp*T0) * 1e5, zarr, label=false, palette=palette([:blue, :red], 11))
    plot!(q, qarr*g * 1e5, zarr, label=false, palette=palette([:blue, :red], 11))

    # file = path*"co2_downstep_"*string(co2)*".jld2"
    # dat = load(file)
    # uf = dat["uf"];
    # par = dat["p"];
    # zarr, sarr, qarr, barr = flux_profile(uf, par)
    # plot!(b, barr * 1e5, zarr, ls=:dot, lw=2, label=false, palette=palette([:blue, :red], 11))
    # plot!(sl, sarr*g/(Cp*T0) * 1e5, zarr, ls=:dot, lw=2, label=false, palette=palette([:blue, :red], 11))
    # plot!(q, qarr*g * 1e5, zarr, ls=:dot, lw=2, label=false, palette=palette([:blue, :red], 11))
end
p = plot(b, sl, q, layout=(1,3))
# savefig(p, "flux_profiles.png")
