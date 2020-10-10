# Thermodynamics

## Saturation adjustment
The liquid water specific humidity is determined according to a standard saturation adjustment procedure. Simply, the mass of condensed water is the excess total water specific humidity exceeding the saturation specific humidity. 

In this mixed-layer model, all thermodynamic equations are written in terms of moist static energy ``h`` and total water specific humidity ``q_t``. The temperature then is calculated implicitly by requiring that 

``h = C_p T + gz + L_0 q_v``.

## Lifting condensation level, cloud base
The cloud base is calculated as the lifting condensation level (LCL), which is a thermodynamic property depending only on the mixed-layer properties. The LCL is defined as the altitude ``z_b`` such that,

``q_{tM} - q_{sat}(z, T(z, h_M, q_{tM})) = 0.``

## Liquid water path
The liquid water path is the mass of condensed liquid water along a vertical path through the cloud, i.e.

``LWP = \int_{z_b}^{z_i} \rho q_l(z) dz``
