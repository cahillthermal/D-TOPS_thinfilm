import argparse
import matplotlib.pyplot as plt
import numpy as np

from calculate_theta import calculate_theta
from datacorrection_complex_leaking import datacorrection_complex_leaking
from fit_thin_film import fit_thin_film
from GetData_out_in_ratio_f_VSUM import GetData_out_in_ratio_f_VSUM


def run_simulation(
    FileNames_data='a1',
    fit_parameters=True,
    fit_freq_range=(1e4,2e5),
    sigma_2_z=0.17,
    alphaT_2=40e-6,
    plot_results=True,
    **model_kwargs,
):
    """
    Runs thin film simulation and fits sigma_2_z and alphaT_2 to experimental data.

    Parameters
    ----------
    FileNames_data : str
        Name of experimental data file.
    fit_parameters : bool
        If True, optimizes sigma_2_z and alphaT_2 to best fit experimental data.
    fit_freq_range : tuple of (float, float) or None
        (f_min, f_max) frequency range for parameter fitting (Hz).
    sigma_2_z : float
        Initial or fixed out-of-plane thermal conductivity of thin film (W/m-K).
    alphaT_2 : float
        Initial or fixed CTE of thin film (1/K).
    plot_results : bool
        If True, generates plots comparing simulation and experiment.
    **model_kwargs : dict
        Additional parameters passed to calculate_theta.

    Returns
    -------
    dict
        Dictionary containing frequencies, model outputs, experimental data, and fitted parameters.
    """
    fit_info = None
    if fit_parameters:
        fit_info = fit_thin_film(
            FileNames_data=FileNames_data,
            fit_freq_range=fit_freq_range,
            initial_sigma_2_z=sigma_2_z,
            initial_alphaT_2=alphaT_2,
            **model_kwargs,
        )
        sigma_2_z = fit_info['sigma_2_z']
        alphaT_2 = fit_info['alphaT_2']
        print("=" * 60)
        print(f"   D-TOPS PARAMETER FITTING RESULTS ({FileNames_data})")
        print("=" * 60)
        if fit_freq_range is not None:
            print(f"Fit Frequency Range : {fit_freq_range[0]:.2e} Hz to {fit_freq_range[1]:.2e} Hz")
        else:
            print("Fit Frequency Range : All experimental frequencies")
        print(f"Optimal sigma_2_z   : {sigma_2_z:.6e} W/m-K")
        print(f"Optimal alphaT_2    : {alphaT_2:.6e} K^-1")
        print(f"Optimization Status : {fit_info['message']}")
        print("=" * 60)

    # Simulation frequencies (logarithmic grid)
    f = np.logspace(np.log10(2e3), np.log10(200e3), 40)

    theta, theta_deformation, theta_dndT_air, det_factor = calculate_theta(
        f, sigma_2_z=sigma_2_z, alphaT_2=alphaT_2, **model_kwargs
    )

    # Read experimental data for comparison
    V_out_data, V_in_data, V_ratio_data, V_SUM_data, fd = (
        GetData_out_in_ratio_f_VSUM(FileNames_data)
    )
    V_SUM = np.mean(V_SUM_data) * 4.0

    # Leaking correction parameters
    Amplitude_corrected_3 = -6.08e-09
    Amplitude_corrected_2 = 1.18e-06
    Amplitude_corrected_1 = 1.50e-04
    Amplitude_corrected_0 = 1.0

    delay_2 = 7.33e-12
    delay_1 = -1.10e-05
    delay_0 = 5.08e-03

    complex_leaking = (
        Amplitude_corrected_0
        + Amplitude_corrected_1 * np.sqrt(fd)
        + Amplitude_corrected_2 * np.sqrt(fd) ** 2
        + Amplitude_corrected_3 * np.sqrt(fd) ** 3
    ) * np.exp(1j * (delay_0 + delay_1 * fd + delay_2 * fd**2))

    Vcorrected_in, Vcorrected_out, Vcorrected_ratio = (
        datacorrection_complex_leaking(
            V_out_data, V_in_data, complex_leaking
        )
    )

    V_signal_deformation = det_factor * theta_deformation * V_SUM / np.sqrt(2.0)
    V_signal_dndT_air = det_factor * theta_dndT_air * V_SUM / np.sqrt(2.0)
    V_signal = det_factor * theta * V_SUM / np.sqrt(2.0)

    if plot_results:
        # Plot 1: Probe beam deflection angle theta
        complex_theta_to_plot = theta

        fig1 = plt.figure(1, figsize=(10, 4.5))
        plt.subplot(1, 2, 1)
        plt.semilogx(
            f, 1e6 * np.real(complex_theta_to_plot), 'k-', linewidth=1.5, label='In-phase'
        )
        plt.semilogx(
            f, 1e6 * np.imag(complex_theta_to_plot), 'k--', linewidth=1.5, label='Out-of-phase'
        )
        plt.xlabel('f (Hz)')
        plt.ylabel('in, out-of-phase (μrad)')
        plt.grid(True)
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.loglog(
            f,
            -np.real(complex_theta_to_plot) / np.imag(complex_theta_to_plot),
            'k-',
            linewidth=1.5,
        )
        plt.xlabel('f (Hz)')
        plt.ylabel('ratio')
        plt.grid(True)
        plt.tight_layout()

        # Plot 2: Lock-in amplifier signal V_signal
        complex_V_signal_to_plot = V_signal

        fig2 = plt.figure(2, figsize=(10, 4.5))
        if fit_parameters:
            fig2.suptitle(
                f"Fitting Results ({FileNames_data}): $\\sigma_{{2,z}}$ = {sigma_2_z:.4e} W/m-K, "
                f"$\\alpha_{{T,2}}$ = {alphaT_2:.4e} K$^{{-1}}$\n"
                f"Fit Frequency Range: {f'{fit_freq_range[0]:.1e} to {fit_freq_range[1]:.1e} Hz' if fit_freq_range else 'All frequencies'}",
                fontsize=11,
                fontweight='bold',
            )

        plt.subplot(1, 2, 1)
        plt.semilogx(
            f, 1e6 * np.real(complex_V_signal_to_plot), 'k-', linewidth=1.5, label='Sim In-phase'
        )
        plt.semilogx(
            f, 1e6 * np.imag(complex_V_signal_to_plot), 'k--', linewidth=1.5, label='Sim Out-of-phase'
        )
        plt.semilogx(
            fd, 1e6 * Vcorrected_in, 'b-', linewidth=1.5, label='Exp In-phase'
        )
        plt.semilogx(
            fd, 1e6 * Vcorrected_out, 'b--', linewidth=1.5, label='Exp Out-of-phase'
        )
        plt.xlabel('f (Hz)')
        plt.ylabel('in, out-of-phase (μV)')
        plt.grid(True)
        plt.legend()

        plt.subplot(1, 2, 2)
        sim_label = 'Sim Ratio (Fitted)' if fit_parameters else 'Sim Ratio'
        plt.loglog(
            f,
            -np.real(complex_theta_to_plot) / np.imag(complex_theta_to_plot),
            'k-',
            linewidth=1.5,
            label=sim_label,
        )
        plt.loglog(fd, Vcorrected_ratio, 'b-', linewidth=1.5, label='Exp Ratio')

        if fit_parameters and fit_freq_range is not None:
            f_min, f_max = fit_freq_range
            if f_min is not None:
                plt.axvline(f_min, color='r', linestyle=':', label='Fit min f')
            if f_max is not None:
                plt.axvline(f_max, color='r', linestyle='--', label='Fit max f')

        plt.xlabel('f (Hz)')
        plt.ylabel('ratio')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

    return {
        'f': f,
        'theta': theta,
        'V_signal': V_signal,
        'fd': fd,
        'Vcorrected_in': Vcorrected_in,
        'Vcorrected_out': Vcorrected_out,
        'Vcorrected_ratio': Vcorrected_ratio,
        'sigma_2_z': sigma_2_z,
        'alphaT_2': alphaT_2,
        'fit_info': fit_info,
    }


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='D-TOPS Thin Film Model Simulation and Parameter Fitting'
    )
    parser.add_argument('--data', type=str, default='a1', help='Data file name (default: a1)')
    parser.add_argument('--fit', type=str2bool, nargs='?', const=True, default=None, help='Fit sigma_2_z and alphaT_2 (e.g. --fit True or --fit False)')
    parser.add_argument('--fmin', type=float, default=None, help='Minimum frequency for fitting (Hz)')
    parser.add_argument('--fmax', type=float, default=None, help='Maximum frequency for fitting (Hz)')
    parser.add_argument('--sigma_2_z', type=float, default=0.17, help='Initial/fixed thin film thermal conductivity (W/m-K)')
    parser.add_argument('--alphaT_2', type=float, default=40e-6, help='Initial/fixed thin film thermal expansion coefficient (1/K)')

    args = parser.parse_args()

    fit_range = None
    if args.fmin is not None or args.fmax is not None:
        fit_range = (args.fmin, args.fmax)

    run_kwargs = {
        'FileNames_data': args.data,
        'fit_freq_range': fit_range,
        'sigma_2_z': args.sigma_2_z,
        'alphaT_2': args.alphaT_2,
    }
    if args.fit is not None:
        run_kwargs['fit_parameters'] = args.fit

    results = run_simulation(**run_kwargs)
    plt.show()
