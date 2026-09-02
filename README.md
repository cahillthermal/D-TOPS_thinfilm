# D-TOPS Thin Film Thermal & Thermoelastic Modeling Engine

A Python modeling and parameter-fitting package for Dual Thermo-Optic and Photo-acoustic / Photothermal Beam Deflection Spectroscopy (**D-TOPS**) of thin films on substrates.

This software simulates thermal transport, thermoelastic surface deformation, and surrounding gas photothermal refractive index variation (mirage effect) caused by modulated laser heating. It provides optimization routines to extract out-of-plane thermal conductivity ($\sigma_{2,z}$) and linear coefficient of thermal expansion ($\alpha_{T,2}$) from experimental lock-in amplifier measurements.

---

## Table of Contents

- [Overview & Physics Model](#overview--physics-model)
- [System Geometry](#system-geometry)
- [Signal Mechanism Analysis](#signal-mechanism-analysis)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Usage](#usage)
  - [Command-Line Interface (CLI)](#command-line-interface-cli)
  - [Python API Usage](#python-api-usage)
- [Data File Format](#data-file-format)
- [Mathematical Formulation](#mathematical-formulation)

---

## Overview & Physics Model

D-TOPS measurements utilize a modulated pump laser to heat a sample surface and a probe laser beam to detect small angular deflections ($\theta$) resulting from:
1. **Surface Thermoelastic Deformation**: Thermal expansion of the thin film, metal transducer, and underlying substrate causing local surface distortion.
2. **Photothermal Mirage Effect**: Temperature gradients in the adjacent air/gas layer inducing a spatially non-uniform refractive index gradient ($dn/dT$).

The mathematical framework solves coupled Hankel-transformed heat conduction and elastic wave equations across a multi-layer stack, accounting for through-thickness and in-plane thermal conductivity anisotropy ($\sigma_{2,z}, \sigma_{2,r}$) as well as anisotropic/textured elastic constants in transducer layers.

---

## System Geometry

The modeled domain consists of a 4-layer structure:

| Layer Index | Layer Description | Example Material | Key Physical Properties |
|-------------|-------------------|------------------|-------------------------|
| **Layer 4** | Ambient Environment ($z > 0$) | Air / Vacuum | Thermal conductivity ($\sigma_4$), Heat capacity ($C_4$), Refractive index derivative ($dn/dT$) |
| **Layer 1** | Metal Transducer Absorber ($0 > z > -L_1$) | Nb-V, Al, Ni ($54\text{ nm}$) | Thermal conductivity ($\sigma_1$), CTE ($\alpha_{T,1}$), Textured Elastic Constants ($C_{ij}$) |
| **Layer 2** | Thin Film Under Study ($-L_1 > z > -(L_1 + L_2)$) | Polymer, Oxide ($200\text{ nm}$) | Anisotropic thermal conductivity ($\sigma_{2,z}, \sigma_{2,r}$), CTE ($\alpha_{T,2}$), Elastic modulus ($E_2, \nu_2$) |
| **Layer 3** | Substrate ($z < -(L_1 + L_2)$) | Silicon (Si) | Thermal conductivity ($\sigma_3$), CTE ($\alpha_{T,3}$), Elastic Constants ($C_{11,3}, C_{12,3}$) |

---

## Signal Mechanism Analysis

For a representative $200\text{ nm}$ polymer film coated with $54\text{ nm}$ Nb-V metal transducer on a Silicon substrate subjected to a $3\text{ mW}$ pump beam with a $10\times$ objective at $1\text{ kHz}$ modulation frequency:

### Signal Breakdown (In-Phase Signal at 1 kHz)

| Simulation Condition | In-Phase Signal ($\mu\text{V}$) | Relative Contribution |
|----------------------|--------------------------------|-----------------------|
| **Full Model in Air** | **920** | **100%** |
| **Vacuum / No Air Mirage Effect** | 670 | ~73% (Pure Elastic Deformation) |
| **No Air, Polymer CTE = 0** | 120 | ~13% |
| **No Air, Polymer CTE = 0, Metal CTE = 0** | 26 | ~3% (Si Substrate Contribution) |
| **No Air, Polymer CTE = 0, Si Elastic Constants $\times 10$** | 100 | ~11% |

### Physical Contributions Summary

- **Polymer Film Thermal Expansion**: ~60% of total signal ($550\mu\text{V}$)
- **Air Mirage Effect ($dn/dT$)**: ~27% of total signal ($250\mu\text{V}$)
- **Metal Transducer Thermal Expansion**: ~7% of total signal ($60\mu\text{V}$)
- **Silicon Substrate Thermal Expansion**: ~3% of total signal
- **Silicon Substrate Stress Deformation**: ~2% of total signal

### Transducer Texture Elastic Constants (110 vs 111 Nb-V)

Effective elastic constants for $\text{Nb-V}$ ($110$ texture derived from elastic constants of Nb and V):
- $C_{11} = 215.11\text{ GPa}$
- $C_{12} = 143.41\text{ GPa}$
- $C_{13} = 132.23\text{ GPa}$
- $C_{33} = 226.28\text{ GPa}$
- $C_{44} = 40.06\text{ GPa}$

*(Note: Default simulations use $111$ texture approximations, which yield end results within 10% of the $110$ textured tensor).*

---

## Repository Structure

```
D-TOPS_thinfilm/
├── README.md                           # Documentation and physics summary
├── thin_film_on_iso_substrate.py       # Main Python runner, simulation, CLI & plotting
├── fit_thin_film.py                    # Non-linear least squares parameter optimizer
├── calculate_theta.py                  # Forward physical model engine
├── SensThinfilm.py                     # Sensitivity calculation module
├── GetData_out_in_ratio_f_VSUM.py      # Experimental data loader (.txt)
├── datacorrection_complex_leaking.py   # Instrumental phase/crosstalk correction
├── simpson_inte.py                     # Simpson rule integration for Hankel transforms
├── a1.txt                              # Experimental dataset (Sample A1)
└── pmma_si_b.txt                       # Experimental dataset (PMMA on Si)
```

### Module Overview

- **[`thin_film_on_iso_substrate.py`](thin_film_on_iso_substrate.py)**: Entry point for running forward calculations or fitting routines, correcting data, and displaying interactive plots.
- **[`fit_thin_film.py`](fit_thin_film.py:9)**: Uses [`scipy.optimize.least_squares()`](fit_thin_film.py:109) to fit thin film thermal conductivity ($\sigma_{2,z}$) and CTE ($\alpha_{T,2}$) against lock-in voltage measurements over user-defined frequency bands.
- **[`calculate_theta.py`](calculate_theta.py:6)**: Computes complex deflection angle $\theta(\omega)$ and detection factor (`det_factor`) via Hankel reciprocal space spatial integration ($p$-space).
- **[`datacorrection_complex_leaking.py`](datacorrection_complex_leaking.py:4)**: Applies frequency-dependent polynomial amplitude and phase delay corrections for lock-in signal crosstalk.
- **[`GetData_out_in_ratio_f_VSUM.py`](GetData_out_in_ratio_f_VSUM.py:5)**: Parses 4-column experimental text files into $V_{\text{out}}, V_{\text{in}}, \text{Ratio}, V_{\text{SUM}}, f$.

---

## Installation

### Prerequisites

Ensure Python 3.8+ is installed along with the required scientific software packages:

```bash
pip install numpy scipy matplotlib
```

---

## Usage

### Command-Line Interface (CLI)

Run forward simulations or fit thin film parameters directly from terminal using [`thin_film_on_iso_substrate.py`](thin_film_on_iso_substrate.py:218):

#### 1. Forward Simulation Comparison
Compare baseline model against experimental dataset [`a1.txt`](a1.txt):

```bash
python thin_film_on_iso_substrate.py --data a1
```

#### 2. Parameter Fitting ($\sigma_{2,z}$ and $\alpha_{T,2}$)
Optimize out-of-plane thermal conductivity and CTE over frequency range $10\text{ kHz}$ to $200\text{ kHz}$:

```bash
python thin_film_on_iso_substrate.py --data a1 --fit True --fmin 10000 --fmax 200000
```

#### CLI Options

- `--data <name>`: Experimental data file name (default: `a1`, searches for `<name>.txt`).
- `--fit <bool>`: Enable parameter fitting (`True` or `False`).
- `--fmin <float>`: Minimum frequency limit for optimization in Hz (e.g., `1e4`).
- `--fmax <float>`: Maximum frequency limit for optimization in Hz (e.g., `2e5`).
- `--sigma_2_z <float>`: Initial or fixed thermal conductivity value in $\text{W/(m}\cdot\text{K)}$ (default: `0.17`).
- `--alphaT_2 <float>`: Initial or fixed coefficient of thermal expansion in $\text{K}^{-1}$ (default: `40e-6`).

---

### Python API Usage

You can also import functions into your own custom scripts:

```python
from thin_film_on_iso_substrate import run_simulation
from fit_thin_film import fit_thin_film

# 1. Run parameter optimization
fit_results = fit_thin_film(
    FileNames_data='a1',
    fit_freq_range=(1e4, 2e5),
    initial_sigma_2_z=0.17,
    initial_alphaT_2=40e-6
)

print(f"Fitted Thermal Conductivity: {fit_results['sigma_2_z']:.4f} W/(m K)")
print(f"Fitted CTE: {fit_results['alphaT_2']:.4e} 1/K")

# 2. Run simulation and generate plots
results = run_simulation(
    FileNames_data='a1',
    fit_parameters=True,
    fit_freq_range=(1e4, 2e5),
    plot_results=True
)
```

---

## Data File Format

Experimental data files (such as [`a1.txt`](a1.txt) and [`pmma_si_b.txt`](pmma_si_b.txt)) must be tab-delimited or space-delimited text files containing 4 numerical columns:

| Column 1 | Column 2 | Column 3 | Column 4 |
|----------|----------|----------|----------|
| $V_{\text{in}}$ (Volts) | $V_{\text{out}}$ (Volts) | Modulation Frequency $f$ (Hz) | Sum Signal $V_{\text{SUM}}$ (Volts) |

---

## Mathematical Formulation

The lock-in amplifier signal $V_{\text{signal}}$ is related to the calculated complex probe beam deflection angle $\theta$ by:

$$V_{\text{signal}} = \text{det\_factor} \times \theta \times \frac{V_{\text{SUM}}}{\sqrt{2}}$$

where the geometric detection conversion factor is defined as:

$$\text{det\_factor} = \sqrt{\frac{8}{\pi}} \left(\frac{f_{\text{length}}}{w_{1d}}\right)$$

Here $f_{\text{length}}$ is the objective focal length and $w_{1d}$ is the probe beam radius at the detector. Hankel transformation integrals over spatial frequency $p$ are numerically evaluated using Simpson quadrature routines.
