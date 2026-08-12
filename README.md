Matlab code to analyze D-TOPS data for a thin film on a substrate.  
Geometry is air (or vacuum)/metal-film/thin-film-under-study/substrate.
The only anisotropy that is currently treated in the through-thickness and in-plane anisotropy of the thermal conductivity of the thin film. 
CTE and elastic constants are assumed to isotropic for all layers.
Includes mirage effect in the air layer, and thermal expansion of the other three layers.

Analysis of the relative sizes of signals for an example polymer film.

Changes in the signal for various mechanisms in a thin film D-TOPS measurement.
Use 200 nm, 40 ppm/K, 0.17 W/(m K), 1.4 J/cm^3-K  meant to model A1 sample from Jingyi coated with NbV 54 nm.


Concentrate on in-phase signal at low frequency 1 kHz.  3 mW pump 10x objective. Signal in microvolts

Full model in air	920
No air	670
No air, polymer CTE=0	120
No air, polymer CTE=0, metal CTE=0	26
No air, polymer CTE=0, Si elastic constants x10	100


Metal CTE is 1/5 of polymer so expect that metal contributes 54/100/5=0.11 of what the polymer contributes.   Polymer is 550 so metal should contribute 60.  Metal seems to contribute 74. Poisson ratio isn't calculated correctly so this is probably OK. Only 1% off in terms of total signal

Asked AI to calculate the effective elastic constants of 110 texture Nb-V from the elastic constants of Nb and V.  AI's answer is 
C11 = 215.11 GPa
C12 = 143.41
C13 = 132.23
C33 = 226.28
C44 = 40.06

The values that are used now are 111 texture but the end result is almost the same, within 10%


Bottom line:
thermal expansion of Si contributes about 3%
Deformation of Si due to metal film stress is about 2%
Direct thermal expansion of metal is about 7%
Air mirage effect is 27%
Polymer is 60%
