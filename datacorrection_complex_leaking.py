import numpy as np


def datacorrection_complex_leaking(Vout_data, Vin_data, complex_leaking):
    """
    Corrects measured lock-in data using complex leaking values.

    Parameters
    ----------
    Vout_data : np.ndarray
    Vin_data : np.ndarray
    complex_leaking : np.ndarray

    Returns
    -------
    Vcorrected_in, Vcorrected_out, Vcorrected_ratio : np.ndarray
    """
    Vcomplex_data = Vin_data + 1j * Vout_data
    Vcorrected_complex = Vcomplex_data / complex_leaking
    Vcorrected_in = np.real(Vcorrected_complex)
    Vcorrected_out = np.imag(Vcorrected_complex)
    Vcorrected_ratio = -Vcorrected_in / Vcorrected_out

    return Vcorrected_in, Vcorrected_out, Vcorrected_ratio
