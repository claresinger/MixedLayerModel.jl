exp_path = "critical_co2/";

# create parameters
using MixedLayerModel

function co2_loop(par)
     # 400 ppm
     par.stype = fixSST();
     par.CO2 = 400;
     u0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
     uf = sol.u[end];
     zb = calc_LCL(uf);
     LWP = incloud_LWP(uf, zb);
     OHU_400 = calc_OHU(uf,par,LWP,par.stype);
 
     # upsteps
     par.stype = varSST();
     for co2i in 600:200:2400
        par.CO2 = co2i;
        par.OHU = OHU_400;
        
        u0, sol = run_mlm_from_init(uf, par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax), quiet=true);
        uf = sol.u[end];
        println(co2i, "\t", uf[5])
        if uf[5] < 0.3
            return co2i
        end
    end
    return NaN
end

par = upCO2();
par.etype = enBal();
par.fttype = co2EIS();
par.rtype = varRad();
par.stype = fixSST();
dt, tmax = 10*24.0, 100.0;

# adjust tunable parameters
par.Cd = 8e-4;
par.α_vent = 1e-3;
par.EIS0 = 8.0;
par.ECS = 3.0;
par.Eexport = 15.0;
par.SW_b = 160;

x_list = 0.5e-3:0.1e-3:1.5e-3
co2c_list = zeros(length(x_list))
for (i,x) in enumerate(x_list)
    par.α_vent = x
    co2c_list[i] = co2_loop(par)
    println("co2 crit = ", co2c_list[i])
    println()
end
par.α_vent = 1e-3

plot(figsize=(6,4), dpi=300)
plot!(x_list, co2c_list, marker=:circle, 
    xlabel="α_vent", ylabel="Critical CO2 [ppm]", label=false)
savefig("experiments/figures/co2_crit_α.png")


x_list = 7e-4:0.5e-4:10e-4
co2c_list = zeros(length(x_list))
for (i,x) in enumerate(x_list)
    par.Cd = x
    co2c_list[i] = co2_loop(par)
    println("co2 crit = ", co2c_list[i])
    println()
end
par.Cd = 8e-4

plot(figsize=(6,4), dpi=300)
plot!(x_list, co2c_list, marker=:circle, 
    xlabel="Cd", ylabel="Critical CO2 [ppm]", label=false)
savefig("experiments/figures/co2_crit_Cd.png")

x_list = 100:10:200
co2c_list = zeros(length(x_list))
for (i,x) in enumerate(x_list)
    par.SW_b = x
    co2c_list[i] = co2_loop(par)
    println("co2 crit = ", co2c_list[i])
    println()
end
par.SW_b = 160

plot(figsize=(6,4), dpi=300)
plot!(x_list, co2c_list, marker=:circle, 
    xlabel="SW b coeff", ylabel="Critical CO2 [ppm]", label=false)
savefig("experiments/figures/co2_crit_SWb.png")
