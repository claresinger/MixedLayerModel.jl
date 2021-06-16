var documenterSearchIndex = {"docs":
[{"location":"exp/#Running-an-experiment","page":"Running an experiment","title":"Running an experiment","text":"","category":"section"},{"location":"exp/","page":"Running an experiment","title":"Running an experiment","text":"First use the MLM to calculate the steady-state clouds given an atmosphere with 400 ppm CO_2. Save this result to a file. In this example we use the energy balance entrainment parameterization and have interactive radiation, but fixed sea-surface temperatures.","category":"page"},{"location":"exp/","page":"Running an experiment","title":"Running an experiment","text":"using MixedLayerModel\nusing FileIO\n\n# run simulation with 400 ppm CO2\npar = basic_params();\npar.etype = enBal();\npar.rtype = varRad();\nu0, sol = run_mlm(par);\n\n# get output\ncode = sol.retcode;\nprintln(code);\nuf = sol.u;\ndu = zeros(4);\nmlm(du, uf, par, 0.0);\nzi,hM,qM,SST = uf;\nzb = calc_LCL(zi,hM,qM);\n\n# save output\noutput = Dict(\"code\" => code, \"p\"=>par, \"u0\" => u0, \"uf\" => uf, \"du/u\" => du./uf, \n\"we\" => we(uf,par,par.etype), \"zb\" => zb, \"zc\" => zi-zb,\n\"RHsurf\" => RH(0.0, hM, qM), \"LHF\" => calc_LHF(uf,par), \"SHF\" => calc_SHF(uf,par),\n\"ΔR\" => calc_cloudtop_RAD(uf,par,par.rtype), \"OHU\" => calc_OHU(uf,par,par.stype))\npath = \"experiments/output/exp_name/\";\nsave(path*\"co2_400.jld2\", output)","category":"page"},{"location":"exp/","page":"Running an experiment","title":"Running an experiment","text":"Now we can load in the results from the 400 ppm simulation and use the calculated steady-state to initialize a new simulation with 800 ppm CO_2. In this example we still have fixed SSTs. ","category":"page"},{"location":"exp/","page":"Running an experiment","title":"Running an experiment","text":"using MixedLayerModel\nusing FileIO\n\n# load initial condition from file\npath = \"experiments/output/exp_name/\";\noutput = load(path*\"co2_400.jld2\");\nu0 = output[\"uf\"];\n\n# run experiment now with 800 ppm CO2\npar = basic_params();\nnewCO2 = 800.0;\npar.CO2 = newCO2;\npar.etype = enBal();\npar.rtype = varRad();\nu0, sol = run_mlm_from_init(u0, par);\n\n# get output\ncode = sol.retcode;\nprintln(code);\nuf = sol.u;\ndu = zeros(4);\nmlm(du, uf, par, 0.0);\nzi,hM,qM,SST = uf;\nzb = calc_LCL(zi,hM,qM);\n\n# save output\noutput = Dict(\"code\" => code, \"p\" => par, \"u0\" => u0, \"uf\" => uf, \"du/u\" => du./uf, \n\"we\" => we(uf,par,par.etype), \"zb\" => zb, \"zc\" => zi-zb,\n\"RHsurf\" => RH(0.0, hM, qM), \"LHF\" => calc_LHF(uf,par), \"SHF\" => calc_SHF(uf,par),\n\"ΔR\" => calc_cloudtop_RAD(uf,par,par.rtype), \"OHU\" => calc_OHU(uf,par,par.stype))\nsave(path*\"co2_upstep_fixSST_\"*string(Int(newCO2))*\".jld2\", output)","category":"page"},{"location":"exp/","page":"Running an experiment","title":"Running an experiment","text":"Instead if we want to allow SSTs to change interactively, then we can set the ocean heat uptake (OHU) to the value defined by the 400 ppm simulation. Also be sure to switch to varSST() mode.","category":"page"},{"location":"exp/","page":"Running an experiment","title":"Running an experiment","text":"using MixedLayerModel\nusing FileIO\n\n# load initial condition from file\npath = \"experiments/output/exp_name/\";\noutput = load(path*\"co2_400.jld2\");\nu0 = output[\"uf\"];\nOHU = output[\"OHU\"];\n\n# run experiment now with 800 ppm CO2\npar = basic_params();\nnewCO2 = 800.0;\npar.CO2 = newCO2;\npar.OHU = OHU;\npar.etype = enBal();\npar.rtype = varRad();\npar.stype = varSST();\nu0, sol = run_mlm_from_init(u0, par);\n\n# get output\ncode = sol.retcode;\nprintln(code);\nuf = sol.u;\ndu = zeros(4);\nmlm(du, uf, par, 0.0);\nzi,hM,qM,SST = uf;\nzb = calc_LCL(zi,hM,qM);\n\n# save output\noutput = Dict(\"code\" => code, \"p\" => par, \"u0\" => u0, \"uf\" => uf, \"du/u\" => du./uf, \n\"we\" => we(uf,par,par.etype), \"zb\" => zb, \"zc\" => zi-zb,\n\"RHsurf\" => RH(0.0, hM, qM), \"LHF\" => calc_LHF(uf,par), \"SHF\" => calc_SHF(uf,par),\n\"ΔR\" => calc_cloudtop_RAD(uf,par,par.rtype), \"OHU\" => calc_OHU(uf,par,par.stype))\nsave(path*\"co2_upstep_\"*string(Int(newCO2))*\".jld2\", output)","category":"page"},{"location":"results/#MLM-behavior-with-fixed-SST","page":"Results","title":"MLM behavior with fixed SST","text":"","category":"section"},{"location":"results/","page":"Results","title":"Results","text":"(Image: ) Fig. 1. xxx","category":"page"},{"location":"results/","page":"Results","title":"Results","text":"(Image: ) Fig. 2. xxx","category":"page"},{"location":"theory/#Mixed-Layer-Theory","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"","category":"section"},{"location":"theory/#Bulk-boundary-layers","page":"Mixed Layer Theory","title":"Bulk boundary layers","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The mixed-layer model (MLM) is based on the assumption that the boundary layer between the surface and inversion height z_i is well-mixed. With these assumptions the governing equations simplify to a set of coupled ODEs for the inversion height and two thermodynamic variables for the energy and the water in the system.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The two thermodynamic quantities (psi) we use are h = C_p T + gz + L_v q_v, the moist static energy, and q_t = q_v + q_l, the total water specific humidity. Fig. 1 shows a sketch of the profiles of moist static energy, total water specific humidity, relative humidity, and liquid water specific humidity from the surface into the free-troposphere. The cloud is indicated by the grey shading between altitudes z_b (diagnosed as the lifting condensation level) and z_i where relative humidity is equal to 100%.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"(Image: ) Fig. 1. Sketch of the vertical profiles of psi_1 and psi_2. The grey shading denotes where the cloud layer is predicted to form.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The temporal evolution of these is governed by a balance between turbulent fluxes at the surface and across the inversion, and a diabatic source term Delta F_psi.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"z_i fracdpsi_Mdt = V (psi_0 - psi_M ) + w_e (psi_+ - psi_M) - Delta F_psi","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The prognostic equation for inversion height is found by vertically integrating the continuity equation. The sea surface temperature (SST) is found by a enforcing a closed surface energy budget. In this model we neglect precipitation (Delta F_q_t = 0).","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"beginaligned \n    fracdz_idt = w_e - Dz_i  \n    z_i fracdh_Mdt = V (h_0 - h_M) + w_e (h_+ - h_M) - Delta R  rho  \n    z_i fracdq_tMdt = V (q_t0 - q_tM) + w_e (q_t+ - q_tM)  \n    C fracdSSTdt = (1-alpha) fracS_04 - LW_net - rho V (h_0 - h_M) - OHU \nendaligned","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"To close these equations we must specify the cloud-top entrainment velocity (w_e) and the cloud-top radiative cooling (Delta R).","category":"page"},{"location":"theory/#Cloud-top-entrainment","page":"Mixed Layer Theory","title":"Cloud-top entrainment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"There are two choices for the entrainment parameterization in this mixed-layer model.","category":"page"},{"location":"theory/#Energy-balance-entrainment","page":"Mixed Layer Theory","title":"Energy balance entrainment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The first assumes that the entrainment velocity is such as to satisfy a steady-state energy balance.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"w_e = fracDelta R  rho_refDelta_i s_vl","category":"page"},{"location":"theory/#Buoyancy-flux-entrainment","page":"Mixed Layer Theory","title":"Buoyancy flux entrainment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"This alternative described by Bretherton and Wyant (1997) as the ''minimal'' model calculates the entrainment velocity as being proportional to the average sub-cloud buoyancy flux and inversely proportional to the buoyancy jump across the inversion.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"w_e = frac25 A overlinelangle w s_v rangleDelta_i s_v","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"where s_v is the virtual dry static energy. ","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"This is solved via the method described in Appendix A of Bretherton and Wyant (1997) by splitting the buoyancy flux into two parts and solving each integral analytically and then inverting. ","category":"page"},{"location":"theory/#Radiation","page":"Mixed Layer Theory","title":"Radiation","text":"","category":"section"},{"location":"theory/#Surface-energy-balance","page":"Mixed Layer Theory","title":"Surface energy balance","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The surface energy budget equation can be written as,","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"fracd SSTdt = SW^down + LW^down - SW^up - LW^up- LHF - SHF - OHU","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"We call the first four terms RAD and can write them as","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"beginaligned \n    RAD = SW^down + LW^down - SW^up - LW^up  \n    = (1 - alpha_cloud)(1 - alpha_ocean) fracS_04 + epsilon_cloud sigma T(z_b)^4 - sigma SST^4\nendaligned","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"We then need to define the cloud albedo and cloud emissivity.","category":"page"},{"location":"theory/#Cloud-shortwave-albedo","page":"Mixed Layer Theory","title":"Cloud shortwave albedo","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"alpha_cloud = 1 - fracL_12L_12 + LWP","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"where L_12 is the liquid water path (LWP) value such that alpha_cloud = 05.","category":"page"},{"location":"theory/#Cloud-longwave-emissivity","page":"Mixed Layer Theory","title":"Cloud longwave emissivity","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"epsilon_cloud = 1 - exp(-LWPL_tau)","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"where LWP_tau = 7 g/m^2 is the optical thickness of the cloud.","category":"page"},{"location":"theory/#Cloud-top-longwave-cooling","page":"Mixed Layer Theory","title":"Cloud-top longwave cooling","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The amount of longwave cooling at the cloud top is dependent on the infrared energy radiating up from the cloud and the infrared energy radiating back down from higher in the atmosphere.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"Delta R = epsilon_cloud sigma T(z_i)^4 - sigma T_eff^4","category":"page"},{"location":"theory/#Effective-emissions-temperature-of-downwelling-longwave-radiation-to-cloud-top","page":"Mixed Layer Theory","title":"Effective emissions temperature of downwelling longwave radiation to cloud-top","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"T_eff = a_0 + a_1 ln left( fracCO_2400 right) = 2635 + 108 ln left( fracCO_2400 right)","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"This is an empirical fit to the LES results from Schneider et al. (2019). ","category":"page"},{"location":"theory/#Thermodynamics","page":"Mixed Layer Theory","title":"Thermodynamics","text":"","category":"section"},{"location":"theory/#Saturation-adjustment","page":"Mixed Layer Theory","title":"Saturation adjustment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The liquid water specific humidity is determined according to a standard saturation adjustment procedure. Simply, the mass of condensed water is the excess total water specific humidity exceeding the saturation specific humidity. ","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"In this mixed-layer model, all thermodynamic equations are written in terms of moist static energy h and total water specific humidity q_t. The temperature then is calculated implicitly by requiring that ","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"h = C_p T + gz + L_0 q_v.","category":"page"},{"location":"theory/#Lifting-condensation-level,-cloud-base","page":"Mixed Layer Theory","title":"Lifting condensation level, cloud base","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The cloud base is calculated as the lifting condensation level (LCL), which is a thermodynamic property depending only on the mixed-layer properties. The LCL is defined as the altitude z_b such that,","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"q_tM - q_sat(z T(z h_M q_tM)) = 0","category":"page"},{"location":"theory/#Liquid-water-path","page":"Mixed Layer Theory","title":"Liquid water path","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The liquid water path is the mass of condensed liquid water along a vertical path through the cloud, i.e.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"LWP = int_z_b^z_i rho q_l(z) dz","category":"page"},{"location":"library/#Application-Programming-Interface-(APIs)","page":"APIs","title":"Application Programming Interface (APIs)","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Documenting the public user interface","category":"page"},{"location":"library/#Thermodynamics","page":"APIs","title":"Thermodynamics","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Thermodynamics.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.RH-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.RH","text":"RH(z, h, qt)\n\nrelative humidity is the ratio of \ntotal specific humidity to saturation\nmax value is 1\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_LCL-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.calc_LCL","text":"calc_LCL(zi, hM, qtM)\n\ncalculate the lifiting condensation level\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_LWP-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.calc_LWP","text":"calc_LWP(zi, hM, qtM)\n\ncalulcate the liquid water path\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_qft0-NTuple{4,Any}","page":"APIs","title":"MixedLayerModel.calc_qft0","text":"calc_qft0(RHft, Gamma_q, sft0, Gamma_s)\n\ncalculate the initial free-tropospheric humidity given\nft RH, ft humidity lapse rate, initial ft dry static energy, \nand ft dry static energy lapse rate\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.pres-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.pres","text":"pres(z, T)\n\nassumes hydrostatic balance\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.q_l-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.q_l","text":"q_l(z, T, qt)\n\nas difference between total specific humidity and saturation specific humidity\nif undersaturated, then ql=0\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.q_sat-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.q_sat","text":"q_sat(z, T)\n\nuses Clasius-Clapeyron relation with assumed constant \nlatent heat of vaporization term L0=2.5e6\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.q_v-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.q_v","text":"q_v(z, T, qt)\n\nas the minimum between the total specific humidty and saturation specific humidity\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.rho-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.rho","text":"rho(z, T)    \n\ncalculate density given altitude and temperature\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.temp-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.temp","text":"temp(z, h, qt)\n\nuses saturation adjustment on the enthalpy\nif no zero is found, then set temp = 0°C\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.theta-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.theta","text":"theta(z,h,qt)\n\ncalculate the potential temperature\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.ρref-Tuple{Any}","page":"APIs","title":"MixedLayerModel.ρref","text":"ρref(T)    \n\ncalculate reference density given temperature\nand reference pressure pref\n\n\n\n\n\n","category":"method"},{"location":"library/#Surface-fluxes","page":"APIs","title":"Surface fluxes","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"SurfaceFluxes.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.H_0-Tuple{Any,Any,fixFlux}","page":"APIs","title":"MixedLayerModel.H_0","text":"define the surface enthalpy flux, H_surf\ngiven prescribed sensible and latent heat fluxes\n\nH_surf = (SHF + LHF) / ρref\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.H_0-Tuple{Any,Any,varFlux}","page":"APIs","title":"MixedLayerModel.H_0","text":"define surface enthalpy flux, H_surf\nusing bulk aerodynamic formula\n\nH_surf = C * V * (h0 - h+)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.Q_0-Tuple{Any,Any,fixFlux}","page":"APIs","title":"MixedLayerModel.Q_0","text":"define the surface moisture flux, Q_surf\ngiven prescribed latent heat flux\n\nQ_surf = LHF / (Lv * ρref)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.Q_0-Tuple{Any,Any,varFlux}","page":"APIs","title":"MixedLayerModel.Q_0","text":"define surface moisture flux, Q_surf\nusing bulk aerodynamic formula\n\nQ_surf = C * V * (q0 - q+)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_LHF-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.calc_LHF","text":"calculate the latent heat flux\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_SHF-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.calc_SHF","text":"calculate the sensible heat flux\n\n\n\n\n\n","category":"method"},{"location":"library/#Top-fluxes","page":"APIs","title":"Top fluxes","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"TopFluxes.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.H_zi-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.H_zi","text":"H_zi(u, p)\n\nenthalpy flux into the mixed-layer from above at z=zi\nH_zi = -we * (hft - hM)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.Q_zi-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.Q_zi","text":"Q_zi(u, p)\n\nmoisture flux into the mixed-layer from above at z=zi\nQ_zi = -we * (qft - qM)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.h_ft-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.h_ft","text":"h_ft(z, p)\n\ndefines h+(z) in free troposphere -- given Gamma_s and Gamma_q\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.q_ft-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.q_ft","text":"q_ft(z, p)\n\ndefines qt+(z) in free troposphere -- given Gamma_q\n\n\n\n\n\n","category":"method"},{"location":"library/#Radiation","page":"APIs","title":"Radiation","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Radiation.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.calc_cloudtop_RAD-Tuple{Any,Any,fixRad}","page":"APIs","title":"MixedLayerModel.calc_cloudtop_RAD","text":"returns the prescribed cloud-top radiative cooling ΔR\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_cloudtop_RAD-Tuple{Any,Any,varRad}","page":"APIs","title":"MixedLayerModel.calc_cloudtop_RAD","text":"calculate the net ΔR at cloud-top based on CO2\nbalance between upwelling and downwelling longwave\ndownwelling longwave is based on an effective temperature\nwhich is empirically fit to LES\n\ngives ΔR ≈ 80 W/m2 for 400 ppm CO2\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_surf_RAD-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.calc_surf_RAD","text":"calculate net SW and LW radiation at the surface\n\n\n\n\n\n","category":"method"},{"location":"library/#Entrainment","page":"APIs","title":"Entrainment","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Entrainment.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,Sally}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, etype::Sally)\n\nentrainment velocity obtained via energy balance requirement\nw = a * ΔR / (Δs_vli * ρref)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,bflux}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, etype::bflux)\n\nentrainment velocity based on buoyancy flux\nwithout radiation\n\nintegral is calculated analytically\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,enBal}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, etype::enBal)\n\nentrainment velocity obtained via energy balance requirement\nw = ΔR / (Δs_vli * ρref)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,fixed}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, etype::fixed)\n\nfixed entrainment velocity of 7 mm/s\n\n\n\n\n\n","category":"method"},{"location":"library/#MLM-ODE","page":"APIs","title":"MLM ODE","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"MLMode.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.mlm-NTuple{4,Any}","page":"APIs","title":"MixedLayerModel.mlm","text":"mlm(du, u, p, t)    \n\ndefine the coupled ODE\n  dzi/dt = D*zi - we\n  dhM/dt = -dE/dz = 1/zi * (Hzi - H0 + dR/rho)\n  dqM/dt = -dW/dz = 1/zi * (Qzi - Q0)\n  dSST/dt = 1/c * (SWnet - LWnet - SHF - LHF - OHU)\n\n\n\n\n\n","category":"method"},{"location":"library/#MLM-solve","page":"APIs","title":"MLM solve","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"MLMsolve.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.run_mlm-Tuple{Any}","page":"APIs","title":"MixedLayerModel.run_mlm","text":"run_mlm(params, filename=\"default.txt\")\n\nrun MLM simulation with given parameters\nand save output to file\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.run_mlm_from_init","page":"APIs","title":"MixedLayerModel.run_mlm_from_init","text":"run_mlm_from_init(u0, params, filename=\"default.txt\")\n\nrun MLM simulation with given parameters\nfrom initial state u0\nand save output to file\n\n\n\n\n\n","category":"function"},{"location":"library/#Diagnostics","page":"APIs","title":"Diagnostics","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Diagnostics.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.calc_OHU-Tuple{Any,Any,fixSST}","page":"APIs","title":"MixedLayerModel.calc_OHU","text":"calc_OHU(u, p, p.stype::fixSST)\n\nu is the state vector [zi, hM, qM, SST]\np is the parameter object\n\ncalculates the ocean heat uptake as the residual \nbetween the radiative fluxes and the LHF + SHF\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_OHU-Tuple{Any,Any,varSST}","page":"APIs","title":"MixedLayerModel.calc_OHU","text":"calc_OHU(u, p, p.stype::varSST)\n\nu is the state vector [zi, hM, qM, SST]\np is the parameter object\n\nreturns p.OHU\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_bflux-Tuple{Any,Any,Any,bflux}","page":"APIs","title":"MixedLayerModel.calc_bflux","text":"calc_bflux(u, p, zarr, etype::bflux)\n\nu is the state vector [zi, hM, qM, SST]\np is the parameter object\nzarr is an array of altitudes from 0 to some maxz > zi\n\ncalculates the buoyancy flux for plotting\n\n\n\n\n\n","category":"method"},{"location":"#MixedLayerModel.jl","page":"Home","title":"MixedLayerModel.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MixedLayerModel","category":"page"},{"location":"","page":"Home","title":"Home","text":"MixedLayerModel.jl is a library for making generic mixed layer models for the atmosphere and subtropical stratocumulus clouds. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Based on Bretherton and Wyant (1997).","category":"page"},{"location":"","page":"Home","title":"Home","text":"Visit on github: MixedLayerModel.jl","category":"page"}]
}
