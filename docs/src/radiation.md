# Radiation

## Surface energy balance
The surface energy budget equation can be written as,

``\frac{d SST}{dt} = SW^{down} + LW^{down} - SW^{up} - LW^{up}- LHF - SHF - OHU.``

We call the first four terms ``RAD`` and can write them as

``RAD = SW^{down} + LW^{down} - SW^{up} - LW^{up} = (1 - \alpha_{cloud})(1 - \alpha_{ocean}) \frac{S_0}{4} + \epsilon_{cloud} \sigma T(z_b)^4 - \sigma SST^4``

We then need to define the cloud albedo and cloud emissivity.

### Cloud shortwave albedo
``\alpha_{cloud} = 1 - \frac{LWP_{1/2}}{LWP_{1/2} + LWP}``

where ``LWP_{1/2}`` is the LWP value such that ``\alpha_{cloud} = 0.5``.

### Cloud longwave emissivity 
``\epsilon_{cloud} = 1 - \exp(-LWP/LWP_0)``

where ``LWP_0 = 7`` g/m``^2`` is the ``lifetime.''

## Cloud-top longwave cooling 
The amount of longwave cooling at the cloud top is dependent on the infrared energy radiating up from the cloud and the infrared energy radiating back down from higher in the atmosphere.

``\Delta R = \epsilon_{cloud} \sigma T(z_i)^4 - \epsilon_{a} \sigma T(z_a)^4``

### Atmospheric longwave emissivity
``\epsilon = 0.8 \frac{\log(CO_2)}{\log(400)}``

This is an empirical fit. 