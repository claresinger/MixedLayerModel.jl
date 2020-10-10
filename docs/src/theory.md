# Mixed Layer Theory

## Bulk boundary layers

## Cloud-top entrainment
There are two choices for the entrainment parameterization in this mixed-layer model.

***Energy balance entrainment***  

The first assumes that the entrainment velocity is such as to satisfy a steady-state energy balance.

``w_e = \frac{\Delta R / \rho_{ref}}{\Delta_i s_{vl}}``

***Buoyancy flux entrainment***
This alternative described by Bretherton and Wyant (1997) as the ''minimal'' model calculates the entrainment velocity as being proportional to the average sub-cloud buoyancy flux and inversely proportional to the buoyancy jump across the inversion.

``w_e = \frac{2.5 A \overline{\langle w' s_v' \rangle}}{\Delta_i s_v}``

where ``s_v`` is the virtual dry static energy. 

This is solved via the method described in Appendix A of Bretherton and Wyant (1997) by splitting the buoyancy flux into two parts and solving each integral analytically and then inverting. 

## Assumptions
1. Hydrostatic balance
2. No precipitation
