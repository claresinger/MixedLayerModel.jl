###
# Definitions and constants for atmospheric mixed-layer model
###

# physical constants
g = 9.8             # gravity (m/s^2)
T0 = 273.16         # absolute zero (K)
Rd = 287.0          # gas constant dry air (J/K/kg)
Rv = 461.0          # gas constant water vapor (J/K/kg)
e0 = 610.78         # reference water vapor partial pressure (Pa = J/m^3)
L0 = 2.5e6          # latent heat of vaporization (J/kg)
Cp = 1004.0         # heat capacity at constant pressure (J/K/kg)
σ_SB = 5.67e-8      # stefan-boltzmann constant (W/m^2/K^4)
ρw = 1e3            # density of water (kg/m^3)
Cw = 4.19e3         # (heat capacity of water (J/kg/K) 
Γd = g / Cp         # dry adiabatic lapse rate (K/m)

# definitions
pref = 100000.0     # reference pressure (Pa = J/m^3)
psurf = 101780.0    # surface pressure (Pa = J/m^3)
S_subtr = 471.0     # solar flux in subtropics (W/m^2)
α_ocean = 0.1       # surface albedo of ocean water

# TODO these params change between obs and LES (1 -> 0.2 vs 0.8 -> 0.1)
CFmax = 1           # maximum cloud fraction
CFmin = 0.2         # minimum cloud fraction

# thermodynamic constants
δ = (Rv-Rd)/Rd     
ϵ = Cp*T0/L0        
μ = 1 - δ*ϵ   
β = 0.5         
σ = β*μ - ϵ