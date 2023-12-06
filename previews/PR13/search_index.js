var documenterSearchIndex = {"docs":
[{"location":"exp/#Running-an-experiment","page":"CO_2 Perturbation Experiment","title":"Running an experiment","text":"","category":"section"},{"location":"exp/","page":"CO_2 Perturbation Experiment","title":"CO_2 Perturbation Experiment","text":"First use the MLM to calculate the steady-state clouds given an atmosphere with 400 ppm CO_2 and fix the SST at 290 K.","category":"page"},{"location":"exp/","page":"CO_2 Perturbation Experiment","title":"CO_2 Perturbation Experiment","text":"using MixedLayerModel\nusing OrdinaryDiffEq\nusing Plots\ninclude(\"../../experiments/mlm_solve_funcs.jl\")\n\n# run simulation\npar = upCO2();\npar.etype = enBal();\npar.fttype = co2dep();\npar.rtype = varRad();\npar.stype = fixSST();\ndt = 2.0;\ntmax = 40.0;\nENV[\"GKSwstype\"]=\"nul\"\nu0, sol = run_mlm(par, dt=3600.0*dt, tspan=(0.0,3600.0*24.0*tmax));\n\n# plot #hide\nt = sol.t / 3600.0 / 24.0 #hide\nzi = getindex.(sol.u,1) #hide\nzb = calc_LCL.(sol.u) #hide\nsM = getindex.(sol.u,2) * 1e-3 #hide\nqtM = getindex.(sol.u,3) * 1e3 #hide\nSST = getindex.(sol.u,4) #hide\nCF = getindex.(sol.u,5) #hide\nLWP = incloud_LWP.(sol.u, zb) #hide\nplot(size=(600,500), layout=(5,1)) #hide\nplot!(t, zi, marker=\"o-\", label=\"\", subplot=1, ylabel=\"zi, zb [m]\") #hide\nplot!(t, zb, marker=\"o-\", subplot=1, label=\"\") #hide\nplot!(t, sM, marker=\"o-\", label=\"\", subplot=2, ylabel=\"sM [kJ/kg]\") #hide\nplot!(t, qtM, marker=\"o-\", label=\"\", subplot=3, ylabel=\"qtM [g/kg]\") #hide\nplot!(t, SST, marker=\"o-\", label=\"\", subplot=4, ylabel=\"SST [K]\") #hide\nplot!(t, CF*1e2, marker=\"o-\", label=\"\", subplot=5, ylabel=\"CF [%]\") #hide","category":"page"},{"location":"ode_solver/#Solving-the-coupled-MLM-ODEs","page":"ODE Solver","title":"Solving the coupled MLM ODEs","text":"","category":"section"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"We want to solve the MLM in two scenarios: 1) from a specified initial steady-state condition or 2) from an arbitrary initial guess.","category":"page"},{"location":"ode_solver/#From-initial-condition","page":"ODE Solver","title":"From initial-condition","text":"","category":"section"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"In this first case, we are taking the steady-state solution from a prior simulation of some climate state (usually present-day CO_2 levels of 400 ppm) and perturbing the CO_2 and letting the system evolve to reach a new steady state. We will need to do using DifferentialEquations to define the ODEProblem.","category":"page"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"This is one way we could code this with an explicit timestep of 5 hours, running for 10 days. We are using the ","category":"page"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"function run_mlm_from_init(u0, params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))\n    prob = ODEProblem(mlm, u0, tspan, params);\n    @time begin\n        sol = solve(prob, Euler(), dt=dt, progress=true, progress_steps=50);\n    end\n    return u0, sol\nend","category":"page"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"Alternatively, we could specify this as a SteadyStateProblem like this:","category":"page"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"function run_mlm_ss_from_init(u0, params; dt=3600.0*5.0, tspan=3600.0*24.0*10.0)\n    prob = SteadyStateProblem(mlm, u0, params);\n    tol = 1e-9;\n    @time begin\n        sol = solve(prob, DynamicSS(Euler(); abstol=0.0, reltol=tol, tspan=tspan), dt=dt, progress=true, progress_steps=50);\n    end\n    return u0, sol\nend","category":"page"},{"location":"ode_solver/#No-initial-condition","page":"ODE Solver","title":"No initial condition","text":"","category":"section"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"On the other hand, we sometimes may want to solve the MLM without having an initial condition in mind and we will need to make a good guess to start. We could do that as follows.","category":"page"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"function run_mlm(params; dt=3600.0*5.0, tspan=(0.0,3600.0*24.0*10.0))\n    zi0 = 1200.0\n    qtM0 = params.RHsurf * q_sat(0.0, params.SST0);\n    sM0 = MixedLayerModel.Cp * params.SST0 + MixedLayerModel.L0 * qtM0;\n    u0 = [zi0, sM0, qtM0, params.SST0];\n    prob = ODEProblem(mlm, u0, tspan, params);\n    @time begin\n        sol = solve(prob, Euler(), dt=dt, progress=true, progress_steps=50);\n    end\n    return u0, sol\nend","category":"page"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"And like above, we could also solve this directly for the steady-state using SteadyStateProblem rather than ODEProblem.","category":"page"},{"location":"ode_solver/","page":"ODE Solver","title":"ODE Solver","text":"function run_mlm_ss(params; dt=3600.0*5.0, tspan=3600.0*24.0*10.0)    \n    zi0 = 1200.0;\n    qtM0 = params.RHsurf * q_sat(0.0, params.SST0);\n    sM0 = MixedLayerModel.Cp * params.SST0 + MixedLayerModel.L0 * qtM0;\n    u0 = [zi0, sM0, qtM0, params.SST0];\n\n    prob = SteadyStateProblem(mlm, u0, params);\n    tol = 1e-9;\n    @time begin\n        sol = solve(prob, DynamicSS(Euler(); abstol=0.0, reltol=tol, tspan=tspan), dt=dt, progress=true, progress_steps=50);\n    end\n    return u0, sol\nend","category":"page"},{"location":"theory/#Mixed-Layer-Theory","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"","category":"section"},{"location":"theory/#Bulk-boundary-layers","page":"Mixed Layer Theory","title":"Bulk boundary layers","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The mixed-layer model (MLM) is based on the assumption that the boundary layer between the surface and inversion height z_i is well-mixed. With these assumptions the governing equations simplify to a set of coupled ODEs for the inversion height and two thermodynamic variables for the energy and the water in the system.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The two thermodynamic quantities (psi) we use are h = C_p T + gz + L_v q_v, the moist static energy, and q_t = q_v + q_l, the total water specific humidity. Fig. 1 shows a sketch of the profiles of moist static energy, total water specific humidity, relative humidity, and liquid water specific humidity from the surface into the free-troposphere. The cloud is indicated by the grey shading between altitudes z_b (diagnosed as the lifting condensation level) and z_i where relative humidity is equal to 100%.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"include(\"MakeDiagram.jl\") #hide","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"(Image: ) Fig. 1. Sketch of the vertical profiles of psi_1 and psi_2. The grey shading denotes where the cloud layer is predicted to form.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The temporal evolution of these is governed by a balance between turbulent fluxes at the surface and across the inversion, and a diabatic source term Delta F_psi.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"z_i fracdpsi_Mdt = V (psi_0 - psi_M ) + w_e (psi_+ - psi_M) - Delta F_psi","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The prognostic equation for inversion height is found by vertically integrating the continuity equation. ","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"We extend the traditional MLM from e.g. Bretherton and Wyant (1997) to couple it to a slab ocean to ensure a closed surface energy budget and a simple radiation scheme that includes the effects of CO_2 on cloud-top radiative cooling. We also add a prognostic equation for cloud fraction.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The sea surface temperature (SST) is found by a enforcing a closed surface energy budget. In this model we neglect precipitation (Delta F_q_t = 0).","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"beginaligned \n    fracdz_idt = w_e - Dz_i  \n    z_i fracdh_Mdt = V (h_0 - h_M) + w_e (h_+ - h_M) - Delta R  rho  \n    z_i fracdq_tMdt = V (q_t0 - q_tM) + w_e (q_t+ - q_tM)  \n    C fracdSSTdt = (1-alpha) fracS_04 - LW_net - rho V (h_0 - h_M) - OHU \n    fracdCFdt = fracCF - CFtau_CF \nendaligned","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"To close these equations we must specify the cloud-top entrainment velocity (w_e) and the cloud-top radiative cooling (Delta R) as well as the functional form for the cloud fraction.","category":"page"},{"location":"theory/#Cloud-top-entrainment","page":"Mixed Layer Theory","title":"Cloud-top entrainment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"There are two choices for the entrainment parameterization in this mixed-layer model.","category":"page"},{"location":"theory/#Energy-balance-entrainment","page":"Mixed Layer Theory","title":"Energy balance entrainment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The first assumes that the entrainment velocity is such as to satisfy a steady-state energy balance.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"w_e = fracDelta R  rho_refDelta_i s_vl","category":"page"},{"location":"theory/#Buoyancy-flux-entrainment","page":"Mixed Layer Theory","title":"Buoyancy flux entrainment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"This alternative described by Bretherton and Wyant (1997) as the ''minimal'' model calculates the entrainment velocity as being proportional to the average sub-cloud buoyancy flux and inversely proportional to the buoyancy jump across the inversion.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"w_e = frac25 A overlinelangle w s_v rangleDelta_i s_v","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"where s_v is the virtual dry static energy. ","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"This is solved via the method described in Appendix A of Bretherton and Wyant (1997) by splitting the buoyancy flux into two parts and solving each integral analytically and then inverting. ","category":"page"},{"location":"theory/#Radiation","page":"Mixed Layer Theory","title":"Radiation","text":"","category":"section"},{"location":"theory/#Surface-energy-balance","page":"Mixed Layer Theory","title":"Surface energy balance","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The surface energy budget equation can be written as,","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"fracd SSTdt = SW^down - SW^up + LW^down - LW^up- LHF - SHF - OHU","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"We write the shortwave terms as SW_net = (1 - alpha_cloud)(1 - alpha_ocean) fracS_04","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"We approximate the net longwave radiation as constant LW_net = -30 W/m^2.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"Or we write the net longwave radiation as the difference of two blackbody terms","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"LW_net = sigma (SST - t)^4 - sigma SST^4","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"where t is proportional to the mixed-layer specific humidity because as the air becomes  more moist, the effective emission level gets closer to the surface. We use a  simple fit to LES data and write, t = 500 cdot qM where qM is in kg/kg. ","category":"page"},{"location":"theory/#Cloud-shortwave-albedo","page":"Mixed Layer Theory","title":"Cloud shortwave albedo","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"alpha_cloud = 1 - fracL_12L_12 + LWP","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"where L_12 is the liquid water path (LWP) value such that alpha_cloud = 05. This is based on Stephens (1978b) equations 1 and 7.  We can write L_12 = frac2 mu r_e3 beta where mu = costheta is the cosine of the solar zenith angle,  r_e is the droplet effective radius, and beta is the backscatter coefficient. We take beta = 007 from Table 2, theta = 60^circ, and r_e = 10 mum, which yields a value of L_12 approx 71 g/m^2.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"Alternatively, the cloud albedo can be parameterized empirically based on the LES results from Schneider et al. (2019) as,","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"alpha_cloud = a left( 1 - fracL_12L_12 + LWP right)","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"where a = 0795 and L_12 = 19136 g/m^2.","category":"page"},{"location":"theory/#Cloud-top-longwave-cooling","page":"Mixed Layer Theory","title":"Cloud-top longwave cooling","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The amount of longwave cooling at the cloud top is dependent on the infrared energy radiating up from the cloud and the infrared energy radiating back down from higher in the atmosphere.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"Delta R = epsilon_cloud sigma T(z_i)^4 - sigma T_eff^4","category":"page"},{"location":"theory/#Effective-emissions-temperature-of-downwelling-longwave-radiation-to-cloud-top","page":"Mixed Layer Theory","title":"Effective emissions temperature of downwelling longwave radiation to cloud-top","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"T_eff = a_0 + a_1 ln left( fracCO_2400 right) = 2635 + 108 ln left( fracCO_2400 right)","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"This is an empirical fit to the LES results from Schneider et al. (2019). ","category":"page"},{"location":"theory/#Cloud-longwave-emissivity","page":"Mixed Layer Theory","title":"Cloud longwave emissivity","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"epsilon_cloud = 1 - exp(-LWPL_tau)","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"where LWP_tau = 7 g/m^2 is the optical thickness of the cloud. This is based on Stephens (1978b) equations 15 and 16, taking an intermediate value of the parameter a_0 = 015.","category":"page"},{"location":"theory/#Cloud-Fraction","page":"Mixed Layer Theory","title":"Cloud Fraction","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"We parameterize the cloud fraction as a function of the stability parameter (aka decoupling parameter), S = left( fracLHFDelta R right) left( fracz_i - z_bz_i right), inspired by Chung and Teixeira (2012).","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"Specifically, we use a smooth function  CF = 1 - frac081 + exp(-m(S-S_crit))  where m=10 is a tunable parameter that sets the strength of the nonlinear feedback and S_crit=07 is the value of the stability parameter that corresponds to CF=06, the halfway point of the transition. The theoretical limit for the stability threshold, where decoupling occurs is at S approx 055.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"include(\"CloudFrac_vs_S.jl\") #hide","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"(Image: )","category":"page"},{"location":"theory/#Thermodynamics","page":"Mixed Layer Theory","title":"Thermodynamics","text":"","category":"section"},{"location":"theory/#Saturation-adjustment","page":"Mixed Layer Theory","title":"Saturation adjustment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The liquid water specific humidity is determined according to a standard saturation adjustment procedure. Simply, the mass of condensed water is the excess total water specific humidity exceeding the saturation specific humidity. ","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"In this mixed-layer model, all thermodynamic equations are written in terms of moist static energy h and total water specific humidity q_t. The temperature then is calculated implicitly by requiring that ","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"h = C_p T + gz + L_0 q_v.","category":"page"},{"location":"theory/#Lifting-condensation-level,-cloud-base","page":"Mixed Layer Theory","title":"Lifting condensation level, cloud base","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The cloud base is calculated as the lifting condensation level (LCL), which is a thermodynamic property depending only on the mixed-layer properties. The LCL is defined as the altitude z_b such that,","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"q_tM - q_sat(z T(z h_M q_tM)) = 0","category":"page"},{"location":"theory/#Liquid-water-path","page":"Mixed Layer Theory","title":"Liquid water path","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The liquid water path is the mass of condensed liquid water along a vertical path through the cloud, i.e.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"LWP = int_z_b^z_i rho q_l(z) dz","category":"page"},{"location":"library/#Application-Programming-Interface-(APIs)","page":"APIs","title":"Application Programming Interface (APIs)","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Documenting the public user interface","category":"page"},{"location":"library/#Thermodynamics","page":"APIs","title":"Thermodynamics","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Thermodynamics.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.calc_LCL-Tuple{Any}","page":"APIs","title":"MixedLayerModel.calc_LCL","text":"calc_LCL(u)\n\ncalculate the lifiting condensation level\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_qft0-NTuple{4,Any}","page":"APIs","title":"MixedLayerModel.calc_qft0","text":"calc_qft0(RHft, Gamma_q, sft0, Gamma_s)\n\ncalculate the initial free-tropospheric humidity given\nft RH, ft humidity lapse rate, initial ft dry static energy, \nand ft dry static energy lapse rate\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.incloud_LWP-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.incloud_LWP","text":"incloud_LWP(u)\n\ncalulcate the in-cloud liquid water path\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.moist_adiabat-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.moist_adiabat","text":"moist_adiabat(Tsurf, zft, p)\ncalculate moist adiabat given a surface temperature (Tsurf),\nup to an altitude zft, with the parameters p\n- first calculates the zLCL\n- then calculates the moist adiabatic profile with dz=10m up to zft\nreturns (T,z) profile\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.pres-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.pres","text":"pres(z, T)\n\nassumes hydrostatic balance\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.q_l-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.q_l","text":"q_l(z, T, qt)\n\nas difference between total specific humidity and saturation specific humidity\nif undersaturated, then ql=0\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.q_sat-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.q_sat","text":"q_sat(z, T)\n\nuses Clasius-Clapeyron relation with assumed constant \nlatent heat of vaporization term L0=2.5e6\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.rho-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.rho","text":"rho(z, T)    \n\ncalculate density given altitude and temperature\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.temp-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.temp","text":"temp(z, s, qt)\n\nuses saturation adjustment on the enthalpy\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.temp_ft-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.temp_ft","text":"calculate actual moist adiabat by integrating\ngo up dry adiabat to LCL and then saturated adiabat\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.Γs-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.Γs","text":"Γs - saturated adiabatic lapse rate\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.ρref-Tuple{Any}","page":"APIs","title":"MixedLayerModel.ρref","text":"ρref(T)    \n\ncalculate reference density given temperature\nand reference pressure pref\n\n\n\n\n\n","category":"method"},{"location":"library/#Surface-fluxes","page":"APIs","title":"Surface fluxes","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"SurfaceFluxes.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.Q_0-Tuple{Any,Any,fixFlux}","page":"APIs","title":"MixedLayerModel.Q_0","text":"define the surface moisture flux, Q_surf\ngiven prescribed latent heat flux\n\nQ_surf = LHF / (Lv * ρref)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.Q_0-Tuple{Any,Any,varFlux}","page":"APIs","title":"MixedLayerModel.Q_0","text":"define surface moisture flux, Q_surf\nusing bulk aerodynamic formula\n\nQ_surf = C * V * (q0 - q)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.S_0-Tuple{Any,Any,fixFlux}","page":"APIs","title":"MixedLayerModel.S_0","text":"define the surface enthalpy flux, H_surf\ngiven prescribed sensible and latent heat fluxes\n\nH_surf = (SHF + LHF) / ρref\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.S_0-Tuple{Any,Any,varFlux}","page":"APIs","title":"MixedLayerModel.S_0","text":"define surface liquid static energy flux, S_surf\nusing bulk aerodynamic formula\n\nS_surf = C * V * (s0 - s)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_LHF-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.calc_LHF","text":"calculate the latent heat flux\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_SHF-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.calc_SHF","text":"calculate the sensible heat flux\n\n\n\n\n\n","category":"method"},{"location":"library/#Top-fluxes","page":"APIs","title":"Top fluxes","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"TopFluxes.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.Q_zi-NTuple{4,Any}","page":"APIs","title":"MixedLayerModel.Q_zi","text":"Q_zi(u, p, ent, LWP)\n\nmoisture flux into the mixed-layer from above at z=zi\nQ_zi = -we * (qft - qM)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.S_zi-NTuple{4,Any}","page":"APIs","title":"MixedLayerModel.S_zi","text":"S_zi(u, p, ent, LWP)\n\nmoist enthalpy flux into the mixed-layer from above at z=zi\nS_zi = -we * (sft - sM)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.qjump-Tuple{Any,Any,Any,co2dep}","page":"APIs","title":"MixedLayerModel.qjump","text":"qjump(u, p, LWP, p.fttype::co2dep)\ndefines the inversion jump for qt \n    via linear regression to LES results\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.qjump-Tuple{Any,Any,Any,fixedFT}","page":"APIs","title":"MixedLayerModel.qjump","text":"qjump(u, p, LWP, p.fttype::fixedFT)\ndefines qt+(z) in free troposphere -- given Gamma_q\nminimum value of qft of 2 g/kg\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.qjump-Tuple{Any,Any,Any,twocol}","page":"APIs","title":"MixedLayerModel.qjump","text":"qjump(u, p, LWP, p.fttype::twocol)\n\nspecific humidity above cloud given fixed RH=0.2\nand saturation calculated at Tft and fixed 1500 m\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.sjump-Tuple{Any,Any,Any,co2dep}","page":"APIs","title":"MixedLayerModel.sjump","text":"sjump(u, p, LWP, p.fttype::co2dep)\ndefines the inversion jump for s\n    via linear regression to LES results\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.sjump-Tuple{Any,Any,Any,fixedFT}","page":"APIs","title":"MixedLayerModel.sjump","text":"sjump(u, p, LWP, p.fttype::fixedFT)\ndefines s+(z) in free troposphere -- given Gamma_s and Gamma_q\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.sjump-Tuple{Any,Any,Any,twocol}","page":"APIs","title":"MixedLayerModel.sjump","text":"sjump(u, p, LWP, p.fttype::twocol)\n\n\n\n\n\n","category":"method"},{"location":"library/#Radiation","page":"APIs","title":"Radiation","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Radiation.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.calc_cloudtop_RAD-Tuple{Any,Any,Any,fixRad}","page":"APIs","title":"MixedLayerModel.calc_cloudtop_RAD","text":"returns the prescribed cloud-top radiative cooling ΔR\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_cloudtop_RAD-Tuple{Any,Any,Any,varRad}","page":"APIs","title":"MixedLayerModel.calc_cloudtop_RAD","text":"calculate the net ΔR at cloud-top based on CO2\nbalance between upwelling and downwelling longwave\ndownwelling longwave is based on an effective temperature\nwhich is empirically fit to LES\n\ngives ΔR ≈ 80 W/m2 for 400 ppm CO2\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_surf_RAD-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.calc_surf_RAD","text":"calculate net SW and LW radiation at the surface\n\n\n\n\n\n","category":"method"},{"location":"library/#Entrainment","page":"APIs","title":"Entrainment","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Entrainment.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.sv_jump-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.sv_jump","text":"sv_jump(u, p, LWP)\n\nΔsv = Δs + Cp[(Rv/Rd - 1)(Tft*qvft - T*qv) + T*ql]\n\njump in virtual liquid static energy across inversion\nproportional to buoyancy jump\nused in energy balance entrainment\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,Any,Any,bflux}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, zb, LWP, etype::bflux)\n\nentrainment velocity based on buoyancy flux\nwithout radiation\n\nintegral is calculated analytically\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,Any,Any,enBal}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, zb, LWP, etype::enBal)\n\nentrainment velocity obtained via energy balance requirement\nw = ΔR / (Δsv * ρref)\nΔsv = the jump in virtual liquid static energy\nsv = Cp*Tv + g*z - Lv*ql = s + Cp(Tv - T)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,Any,Any,fixed}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, zb, LWP, etype::fixed)\n\nfixed entrainment velocity of 7 mm/s\n\n\n\n\n\n","category":"method"},{"location":"library/#MLM-ODE","page":"APIs","title":"MLM ODE","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"MLMode.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.mlm-NTuple{4,Any}","page":"APIs","title":"MixedLayerModel.mlm","text":"mlm(du, u, p, t)    \n\ndefine the coupled ODE\n  dzi/dt = we - D*zi\n  dsM/dt = -dE/dz = 1/zi * (Szi - S0 + dR/rho)\n  dqM/dt = -dW/dz = 1/zi * (Qzi - Q0)\n  dSST/dt = 1/c * (SWnet - LWnet - SHF - LHF - OHU)\n\n\n\n\n\n","category":"method"},{"location":"library/#Diagnostics","page":"APIs","title":"Diagnostics","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Diagnostics.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.calc_OHU-Tuple{Any,Any,Any,fixSST}","page":"APIs","title":"MixedLayerModel.calc_OHU","text":"calc_OHU(u, p, p.stype::fixSST)\n\nu is the state vector [zi, sM, qM, SST, CF]\np is the parameter object\n\ncalculates the ocean heat uptake as the residual \nbetween the radiative fluxes and the LHF + SHF\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_OHU-Tuple{Any,Any,Any,varSST}","page":"APIs","title":"MixedLayerModel.calc_OHU","text":"calc_OHU(u, p, p.stype::varSST)\n\nu is the state vector [zi, sM, qM, SST, CF]\np is the parameter object\n\nreturns p.OHU\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_bflux-Tuple{Any,Any,Any,bflux}","page":"APIs","title":"MixedLayerModel.calc_bflux","text":"calc_bflux(u, p, zarr, etype::bflux)\n\nu is the state vector [zi, sM, qM, SST, CF]\np is the parameter object\nzarr is an array of altitudes from 0 to some maxz > zi\n\ncalculates the buoyancy flux for plotting\n\n\n\n\n\n","category":"method"},{"location":"#MixedLayerModel.jl","page":"Home","title":"MixedLayerModel.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MixedLayerModel","category":"page"},{"location":"","page":"Home","title":"Home","text":"MixedLayerModel.jl is a library for making generic mixed layer models for the atmosphere and subtropical stratocumulus clouds. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Based on Bretherton and Wyant (1997).","category":"page"},{"location":"","page":"Home","title":"Home","text":"Visit on github: MixedLayerModel.jl","category":"page"}]
}
