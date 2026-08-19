import os
import numpy as np


def GetData_out_in_ratio_f_VSUM(file_name):
    """
    Reads experimental data from <file_name>.txt formatted with 4 columns:
    [Vin, Vout, f, V_SUM].

    Parameters
    ----------
    file_name : str
        Name of data file (without extension or with .txt).

    Returns
    -------
    Vout, Vin, Vratio, V_SUM, f : np.ndarray
    """
    if not file_name.endswith('.txt'):
        file_path = f"{file_name}.txt"
    else:
        file_path = file_name

    data = np.loadtxt(file_path)
    # data shape is (N, 4)
    Vin = data[:, 0]
    Vout = data[:, 1]
    f = data[:, 2]
    V_SUM = data[:, 3]
    Vratio = -Vin / Vout

    return Vout, Vin, Vratio, V_SUM, f
