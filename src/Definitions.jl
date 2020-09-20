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
pref = 100000.0     # reference pressure (Pa = J/m^3)
psurf = 101780.0    # surface pressure (Pa = J/m^3)
σ_SB = 5.67e-8      # stefan-boltzmann constant (W/m^2/K^4)
S_subtr = 471.0     # solar flux (W/m^2)
ρw = 1e3            # density of water (kg/m^3)
Cw = 4.19e3         # (heat capacity of water (J/kg/K) 

# thermodynamic constants
δ = (Rv-Rd)/Rd     
ϵ = Cp*T0/L0        
μ = 1 - δ*ϵ   
β = 0.5         
σ = β*μ - ϵ

# define qt+(z) in free troposphere -- given Gamma_q
function q_ft(z, p)
    qft = p.qft0 .+ p.Gamma_q .* z
    qft = max.(qft, 2e-3)
    return qft
end

# define h+(z) in free troposphere -- given Gamma_s and Gamma_q
function h_ft(z, p)
    sft = p.sft0 .+ p.Gamma_s .* z
    hft = sft .* Cp .+ L0 .* q_ft(z, p)
    return hft
end