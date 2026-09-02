import numpy as np
import scipy.optimize

from calculate_theta import calculate_theta
from datacorrection_complex_leaking import datacorrection_complex_leaking
from GetData_out_in_ratio_f_VSUM import GetData_out_in_ratio_f_VSUM


def fit_thin_film(
    FileNames_data='a1',
    fit_freq_range=None,
    initial_sigma_2_z=0.17,
    initial_alphaT_2=40e-6,
    sigma_2_z_bounds=(1e-3, 100.0),
    alphaT_2_bounds=(1e-7, 1e-3),
    **model_kwargs,
):
    """
    Fits sigma_2_z and alphaT_2 to find best agreement between model and experimental data.

    Parameters
    ----------
    FileNames_data : str
        Name of experimental data file (without extension or with .txt).
    fit_freq_range : tuple of (float, float) or None
        (f_min, f_max) frequency range for fitting (Hz). If None, uses all experimental frequencies.
    initial_sigma_2_z : float
        Initial guess for thermal conductivity sigma_2_z (W/m-K).
    initial_alphaT_2 : float
        Initial guess for thermal expansion coefficient alphaT_2 (1/K).
    sigma_2_z_bounds : tuple
        (min, max) bounds for sigma_2_z.
    alphaT_2_bounds : tuple
        (min, max) bounds for alphaT_2.
    **model_kwargs : dict
        Additional parameter overrides passed to calculate_theta.

    Returns
    -------
    dict
        Dictionary with fitted parameters ('sigma_2_z', 'alphaT_2'), optimization details,
        and fitted frequency subset.
    """
    V_out_data, V_in_data, V_ratio_data, V_SUM_data, fd = (
        GetData_out_in_ratio_f_VSUM(FileNames_data)
    )
    V_SUM = np.mean(V_SUM_data) * 4.0

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

    fd = np.asarray(fd)
    mask = np.ones(len(fd), dtype=bool)

    if fit_freq_range is not None:
        if isinstance(fit_freq_range, (tuple, list, np.ndarray)) and len(fit_freq_range) == 2:
            f_min, f_max = fit_freq_range
            if f_min is not None:
                mask &= (fd >= f_min)
            if f_max is not None:
                mask &= (fd <= f_max)
        else:
            raise ValueError("fit_freq_range must be a tuple or list of format (f_min, f_max)")

    fd_fit = fd[mask]
    Vcorrected_in_fit = Vcorrected_in[mask]
    Vcorrected_out_fit = Vcorrected_out[mask]
    ratio_exp_fit = Vcorrected_ratio[mask]

    if len(fd_fit) == 0:
        raise ValueError(f"No experimental frequency data points found in range {fit_freq_range}")

    def residuals(params):
        s2z, aT2 = params
        theta, _, _, det_factor = calculate_theta(
            fd_fit, sigma_2_z=s2z, alphaT_2=aT2, **model_kwargs
        )
        V_signal = det_factor * theta * V_SUM / np.sqrt(2.0)
        V_sim_in = np.real(V_signal)
        V_sim_out = np.imag(V_signal)

        res_in = V_sim_in - Vcorrected_in_fit
        res_out = V_sim_out - Vcorrected_out_fit
        return np.concatenate([res_in, res_out])

    x0 = [initial_sigma_2_z, initial_alphaT_2]
    lower_bounds = [sigma_2_z_bounds[0], alphaT_2_bounds[0]]
    upper_bounds = [sigma_2_z_bounds[1], alphaT_2_bounds[1]]

    res = scipy.optimize.least_squares(
        residuals,
        x0,
        bounds=(lower_bounds, upper_bounds),
        x_scale=[initial_sigma_2_z, initial_alphaT_2],
    )

    sigma_2_z_fit = float(res.x[0])
    alphaT_2_fit = float(res.x[1])

    return {
        'sigma_2_z': sigma_2_z_fit,
        'alphaT_2': alphaT_2_fit,
        'success': res.success,
        'message': res.message,
        'res_obj': res,
        'fd_fit': fd_fit,
        'Vcorrected_in_fit': Vcorrected_in_fit,
        'Vcorrected_out_fit': Vcorrected_out_fit,
        'ratio_exp_fit': ratio_exp_fit,
        'fit_freq_range': fit_freq_range,
    }
