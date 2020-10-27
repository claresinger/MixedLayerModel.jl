# Mixed Layer Theory

## Bulk boundary layers

The mixed-layer model (MLM) is based on the assumption that the boundary layer between the surface and inversion height ``z_i`` is well-mixed. 

The two thermodynamic quantities (``\psi``) we use are ``h = C_p T + gz + L_v q_v``, the moist static energy, and ``q_t = q_v + q_l``, the total water specific humidity. The temporal evolution of these is governed by a balance between turbulent fluxes at the surface and across the inversion, and a diabatic source term ``\Delta F_\psi``.

``z_i \frac{d\psi_{M}}{dt} = V (\psi_0 - \psi_{M} ) + w_e (\psi_+ - \psi_{M}) - \Delta F_\psi``

The prognostic equation for inversion height is found by vertically integrating the continuity equation. The sea surface temperature (SST) is found by a enforcing a closed surface energy budget.

```math
\begin{align} 
    \frac{dz_i}{dt} &= w_e - Dz_i \\ 
    z_i \frac{dh_M}{dt} &= V (h_0 - h_{M}) + w_e (h_+ - h_{M}) - \Delta R / \rho \\ 
    z_i \frac{dq_{tM}}{dt} &= V (q_{t0} - q_{tM}) + w_e (q_{t+} - q_{tM}) \\ 
    C \frac{dSST}{dt} &= (1-\alpha) \frac{S_0}{4} - LW_{net} - \rho V (h_0 - h_M) - OHU 
\end{align}
```

In this model we neglect precipitation (``\Delta F_{q_t} = 0``).

To close these equations we must specify the [cloud-top entrainment velocity](entrainment.md) (``w_e``) and the [cloud-top radiative cooling](radiation.md) (``\Delta R``).