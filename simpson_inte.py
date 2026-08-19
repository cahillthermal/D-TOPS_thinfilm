import numpy as np


def simpson_inte(array, pace):
    """
    1D Simpson's rule numerical integration matching MATLAB simpson_inte.

    Parameters
    ----------
    array : array-like
        Array of integrand values.
    pace : float
        Step size (d_p).

    Returns
    -------
    float or complex
        Integrated value.
    """
    array = np.asarray(array)
    steps = len(array)
    # MATLAB: array(3:2:steps-2) -> Python: array[2:steps-2:2]
    edge_sum = np.sum(array[2 : steps - 2 : 2])
    # MATLAB: array(2:2:steps-1) -> Python: array[1:steps-1:2]
    mid_sum = np.sum(array[1 : steps - 1 : 2])
    I = (1 / 6) * (2 * pace) * (array[0] + 2 * edge_sum + 4 * mid_sum + array[-1])
    return I
