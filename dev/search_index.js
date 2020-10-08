var documenterSearchIndex = {"docs":
[{"location":"radiation/#Radiation","page":"Radiation","title":"Radiation","text":"","category":"section"},{"location":"radiation/#Surface-energy-balance","page":"Radiation","title":"Surface energy balance","text":"","category":"section"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"The surface energy budget equation can be written as,","category":"page"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"fracd SSTdt = SW^down + LW^down - SW^up - LW^up- LHF - SHF - OHU","category":"page"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"We call the first four terms RAD and can write them as","category":"page"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"RAD = SW^down + LW^down - SW^up - LW^up = (1 - alpha_cloud)(1 - alpha_ocean) fracS_04 + epsilon_cloud sigma T(z_b)^4 - sigma SST^4","category":"page"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"We then need to define the cloud albedo and cloud emissivity.","category":"page"},{"location":"radiation/#Cloud-shortwave-albedo","page":"Radiation","title":"Cloud shortwave albedo","text":"","category":"section"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"alpha_cloud = 1 - fracLWP_12LWP_12 + LWP","category":"page"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"where LWP_12 is the LWP value such that alpha_cloud = 05.","category":"page"},{"location":"radiation/#Cloud-longwave-emissivity","page":"Radiation","title":"Cloud longwave emissivity","text":"","category":"section"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"epsilon_cloud = 1 - exp(-LWPLWP_0)","category":"page"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"where LWP_0 = 7 g/m^2 is the ``lifetime.''","category":"page"},{"location":"radiation/#Cloud-top-longwave-cooling","page":"Radiation","title":"Cloud-top longwave cooling","text":"","category":"section"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"The amount of longwave cooling at the cloud top is dependent on the infrared energy radiating up from the cloud and the infrared energy radiating back down from higher in the atmosphere.","category":"page"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"Delta R = epsilon_cloud sigma T(z_i)^4 - epsilon_a sigma T(z_a)^4","category":"page"},{"location":"radiation/#Atmospheric-longwave-emissivity","page":"Radiation","title":"Atmospheric longwave emissivity","text":"","category":"section"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"epsilon = 08 fraclog(CO_2)log(400)","category":"page"},{"location":"radiation/","page":"Radiation","title":"Radiation","text":"This is an empirical fit. ","category":"page"},{"location":"theory/#Mixed-Layer-Theory","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"","category":"section"},{"location":"theory/#Bulk-boundary-layers","page":"Mixed Layer Theory","title":"Bulk boundary layers","text":"","category":"section"},{"location":"theory/#Cloud-top-entrainment","page":"Mixed Layer Theory","title":"Cloud-top entrainment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"There are two choices for the entrainment parameterization in this mixed-layer model.","category":"page"},{"location":"theory/#Energy-balance-entrainment","page":"Mixed Layer Theory","title":"Energy balance entrainment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"The first assumes that the entrainment velocity is such as to satisfy a steady-state energy balance.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"w_e = fracDelta R  rho_refDelta_i s_vl","category":"page"},{"location":"theory/#Buoyancy-flux-entrainment","page":"Mixed Layer Theory","title":"Buoyancy flux entrainment","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"This alternative described by Bretherton and Wyant (1997) as the ''minimal'' model calculates the entrainment velocity as being proportional to the average sub-cloud buoyancy flux and inversely proportional to the buoyancy jump across the inversion.","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"w_e = frac25 A overlinelangle w s_v rangleDelta_i s_v","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"where s_v is the virtual dry static energy. ","category":"page"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"This is solved via the method described in Appendix A of Bretherton and Wyant (1997) by splitting the buoyancy flux into two parts and solving each integral analytically and then inverting. ","category":"page"},{"location":"theory/#Assumptions","page":"Mixed Layer Theory","title":"Assumptions","text":"","category":"section"},{"location":"theory/","page":"Mixed Layer Theory","title":"Mixed Layer Theory","text":"Hydrostatic balance\nNo precipitation","category":"page"},{"location":"library/#Application-Programming-Interface-(APIs)","page":"APIs","title":"Application Programming Interface (APIs)","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Documenting the public user interface","category":"page"},{"location":"library/#Thermodynamics","page":"APIs","title":"Thermodynamics","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Thermodynamics.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.RH-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.RH","text":"RH(z, h, qt)\n\nrelative humidity is the ratio of \ntotal specific humidity to saturation\nmax value is 1\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_LCL-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.calc_LCL","text":"calc_LCL(zi, hM, qtM)\n\ncalculate the lifiting condensation level\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_LWP-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.calc_LWP","text":"calc_LWP(zi, hM, qtM)\n\ncalulcate the liquid water path\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_qft0-NTuple{4,Any}","page":"APIs","title":"MixedLayerModel.calc_qft0","text":"calc_qft0(RHft, Gamma_q, sft0, Gamma_s)\n\ncalculate the initial free-tropospheric humidity given\nft RH, ft humidity lapse rate, initial ft dry static energy, \nand ft dry static energy lapse rate\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.pres-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.pres","text":"pres(z, T)\n\nassumes hydrostatic balance\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.q_l-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.q_l","text":"q_l(z, T, qt)\n\nas difference between total specific humidity and saturation specific humidity\nif undersaturated, then ql=0\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.q_sat-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.q_sat","text":"q_sat(z, T)\n\nuses Clasius-Clapeyron relation with assumed constant \nlatent heat of vaporization term L0=2.5e6\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.q_v-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.q_v","text":"q_v(z, T, qt)\n\nas the minimum between the total specific humidty and saturation specific humidity\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.rho-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.rho","text":"rho(z, T)    \n\ncalculate density given altitude and temperature\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.temp-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.temp","text":"temp(z, h, qt)\n\nuses saturation adjustment on the enthalpy\nif no zero is found, then set temp = 0°C\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.theta-Tuple{Any,Any,Any}","page":"APIs","title":"MixedLayerModel.theta","text":"theta(z,h,qt)\n\ncalculate the potential temperature\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.ρref-Tuple{Any}","page":"APIs","title":"MixedLayerModel.ρref","text":"ρref(T)    \n\ncalculate reference density given temperature\nand reference pressure pref\n\n\n\n\n\n","category":"method"},{"location":"library/#Surface-fluxes","page":"APIs","title":"Surface fluxes","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"SurfaceFluxes.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.H_0-Tuple{Any,Any,fixFlux}","page":"APIs","title":"MixedLayerModel.H_0","text":"define the surface enthalpy flux, H_surf\ngiven prescribed sensible and latent heat fluxes\n\nH_surf = (SHF + LHF) / ρref\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.H_0-Tuple{Any,Any,varFlux}","page":"APIs","title":"MixedLayerModel.H_0","text":"define surface enthalpy flux, H_surf\nusing bulk aerodynamic formula\n\nH_surf = C * V * (h0 - h+)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.Q_0-Tuple{Any,Any,fixFlux}","page":"APIs","title":"MixedLayerModel.Q_0","text":"define the surface moisture flux, Q_surf\ngiven prescribed latent heat flux\n\nQ_surf = LHF / (Lv * ρref)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.Q_0-Tuple{Any,Any,varFlux}","page":"APIs","title":"MixedLayerModel.Q_0","text":"define surface moisture flux, Q_surf\nusing bulk aerodynamic formula\n\nQ_surf = C * V * (q0 - q+)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_LHF-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.calc_LHF","text":"calculate the latent heat flux\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_SHF-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.calc_SHF","text":"calculate the sensible heat flux\n\n\n\n\n\n","category":"method"},{"location":"library/#Top-fluxes","page":"APIs","title":"Top fluxes","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"TopFluxes.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.H_zi-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.H_zi","text":"H_zi(u, p)\n\nenthalpy flux into the mixed-layer from above at z=zi\nH_zi = -we * (hft - hM)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.Q_zi-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.Q_zi","text":"Q_zi(u, p)\n\nmoisture flux into the mixed-layer from above at z=zi\nQ_zi = -we * (qft - qM)\n\n\n\n\n\n","category":"method"},{"location":"library/#Radiation","page":"APIs","title":"Radiation","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Radiation.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.calc_cloudtop_RAD-Tuple{Any,Any,fixRad}","page":"APIs","title":"MixedLayerModel.calc_cloudtop_RAD","text":"calculate the net OLR at cloud-top based on CO2\nthis is the 3-layer atmosphere model\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_cloudtop_RAD-Tuple{Any,Any,varRad}","page":"APIs","title":"MixedLayerModel.calc_cloudtop_RAD","text":"calculate the net OLR at cloud-top based on CO2\nthis is the 3-layer atmosphere model\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.calc_surf_RAD-Tuple{Any,Any}","page":"APIs","title":"MixedLayerModel.calc_surf_RAD","text":"calculate net SW and LW radiation at the surface\n\n\n\n\n\n","category":"method"},{"location":"library/#Entrainment","page":"APIs","title":"Entrainment","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Entrainment.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,bflux}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, etype::bflux)\n\nentrainment velocity based on buoyancy flux\nwithout radiation\n\nintegral is calculated analytically\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,enBal}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, etype::enBal)\n\nentrainment velocity obtained via energy balance requirement\nw = ΔR / (Δs_vli * ρref)\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.we-Tuple{Any,Any,fixed}","page":"APIs","title":"MixedLayerModel.we","text":"we(u, p, etype::fixed)\n\nfixed entrainment velocity of 1.5 mm/s\n\n\n\n\n\n","category":"method"},{"location":"library/#MLM-ODE","page":"APIs","title":"MLM ODE","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"MLMode.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.mlm-NTuple{4,Any}","page":"APIs","title":"MixedLayerModel.mlm","text":"mlm(du, u, p, t)    \n\ndefine the coupled ODE\n\n\n\n\n\n","category":"method"},{"location":"library/#MLM-solve","page":"APIs","title":"MLM solve","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"MLMsolve.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.run-Tuple{Any}","page":"APIs","title":"MixedLayerModel.run","text":"run(params, print=false)\n\nrun MLM simulation with given parameters\n\n\n\n\n\n","category":"method"},{"location":"library/#MixedLayerModel.run_with_output","page":"APIs","title":"MixedLayerModel.run_with_output","text":"run_with_output(params)\n\nrun MLM simulation and print output to file\n\n\n\n\n\n","category":"function"},{"location":"library/#Diagnostics","page":"APIs","title":"Diagnostics","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [MixedLayerModel]\nPrivate = false\nPages   = [\"Diagnostics.jl\"]","category":"page"},{"location":"library/#MixedLayerModel.calc_bflux-Tuple{Any,Any,Any,bflux}","page":"APIs","title":"MixedLayerModel.calc_bflux","text":"calc_bflux(u, p, zarr, etype::bflux)\n\nu is the state vector [zi, hM, qM, SST]\np is the parameter object\nzarr is an array of altitudes from 0 to some maxz > zi\n\ncalculates the buoyancy flux for plotting\n\n\n\n\n\n","category":"method"},{"location":"thermo/#Thermodynamics","page":"Thermodynamics","title":"Thermodynamics","text":"","category":"section"},{"location":"thermo/#Saturation-adjustment","page":"Thermodynamics","title":"Saturation adjustment","text":"","category":"section"},{"location":"thermo/","page":"Thermodynamics","title":"Thermodynamics","text":"The liquid water specific humidity is determined according to a standard saturation adjustment procedure. Simply, the mass of condensed water is the excess total water specific humidity exceeding the saturation specific humidity. ","category":"page"},{"location":"thermo/","page":"Thermodynamics","title":"Thermodynamics","text":"In this mixed-layer model, all thermodynamic equations are written in terms of moist static energy h and total water specific humidity q_t. The temperature then is calculated implicitly by requiring that ","category":"page"},{"location":"thermo/","page":"Thermodynamics","title":"Thermodynamics","text":"h = C_p T + gz + L_0 q_v.","category":"page"},{"location":"thermo/#Lifting-condensation-level,-cloud-base","page":"Thermodynamics","title":"Lifting condensation level, cloud base","text":"","category":"section"},{"location":"thermo/","page":"Thermodynamics","title":"Thermodynamics","text":"The cloud base is calculated as the lifting condensation level (LCL), which is a thermodynamic property depending only on the mixed-layer properties. The LCL is defined as the altitude z such that,","category":"page"},{"location":"thermo/","page":"Thermodynamics","title":"Thermodynamics","text":"q_tM - q_sat(z T(z h_M q_tM)) = 0","category":"page"},{"location":"thermo/#Liquid-water-path","page":"Thermodynamics","title":"Liquid water path","text":"","category":"section"},{"location":"thermo/","page":"Thermodynamics","title":"Thermodynamics","text":"The liquid water path is the mass of condensed liquid water along a vertical path through the cloud, i.e.","category":"page"},{"location":"thermo/","page":"Thermodynamics","title":"Thermodynamics","text":"LWP = int_0^z_i rho q_l(z) dz","category":"page"},{"location":"#MixedLayerModel.jl","page":"Home","title":"MixedLayerModel.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MixedLayerModel","category":"page"},{"location":"","page":"Home","title":"Home","text":"MixedLayerModel.jl is a library for making generic mixed layer models for the atmosphere and subtropical stratocumulus clouds. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Based on Bretherton and Wyant (1997).","category":"page"}]
}
