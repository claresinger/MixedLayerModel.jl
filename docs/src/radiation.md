# Radiation

## Surface energy balance
The surface energy budget equation can be written as,

``\frac{d SST}{dt} = SW^{down} + LW^{down} - SW^{up} - LW^{up}- LHF - SHF - OHU.``

We call the first four terms ``RAD`` and can write them as

```math
\begin{aligned} 
    RAD &= SW^{down} + LW^{down} - SW^{up} - LW^{up} \\ 
    &= (1 - \alpha_{cloud})(1 - \alpha_{ocean}) \frac{S_0}{4} + \epsilon_{cloud} \sigma T(z_b)^4 - \sigma SST^4
\end{aligned}
```

We then need to define the cloud albedo and cloud emissivity.

### Cloud shortwave albedo
``\alpha_{cloud} = 1 - \frac{L_{1/2}}{L_{1/2} + LWP}``

where ``L_{1/2}`` is the LWP value such that ``\alpha_{cloud} = 0.5``.

### Cloud longwave emissivity 
``\epsilon_{cloud} = 1 - \exp(-LWP/L_\tau)``

where ``LWP_\tau = 7`` g/m``^2`` is the optical thickness of the cloud.

## Cloud-top longwave cooling 
The amount of longwave cooling at the cloud top is dependent on the infrared energy radiating up from the cloud and the infrared energy radiating back down from higher in the atmosphere.

``\Delta R = \epsilon_{cloud} \sigma T(z_i)^4 - \sigma T_{eff}^4``

### Effective emissions temperature of downwelling longwave radiation to cloud-top
``T_{eff} = a_0 + a_1 \ln \left( \frac{CO_2}{400} \right) = 263.5 + 10.8 \ln \left( \frac{CO_2}{400} \right)``

This is an empirical fit to the LES results from Schneider et al. (2019). 