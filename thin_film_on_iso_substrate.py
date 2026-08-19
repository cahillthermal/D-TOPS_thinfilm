import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg
import scipy.special

from datacorrection_complex_leaking import datacorrection_complex_leaking
from GetData_out_in_ratio_f_VSUM import GetData_out_in_ratio_f_VSUM
from simpson_inte import simpson_inte


def run_simulation(FileNames_data='a1'):
    # structure from top to bottom:
    # air/Al coating/thin film/substrate
    # index: 4: air (or vacuum); 1: metal film; 2: thin film; 3: substrate

    f = np.logspace(np.log10(2e3), np.log10(200e3), 40)  # modulation frequency (Hz)

    # lens
    lens_magnification = 10.0
    lens_transmittance = 0.87
    focal_length = (5.0 / lens_magnification) * 40e-3

    # optical beam spot size and beam offset
    r_rms = (5.0 / lens_magnification) * 12.8e-6
    xoffset = (5.0 / lens_magnification) * 13.5e-6
    C_probe = 0.90
    w_1_d = 0.90e-3

    # pump beam power
    incident_pump = 3e-3

    # absorbance of pump power
    n_metal = 2.9
    k_metal = 8.2
    sample_reflectance = (
        abs(n_metal - 1 + 1j * k_metal) ** 2
        / abs(n_metal + 1 + 1j * k_metal) ** 2
    )
    sample_absorbance = 1-sample_reflectance
    sample_absorbance = 0.4  # use 0.4 for NbV
    A_pump = (
        incident_pump
        * lens_transmittance
        * sample_absorbance
        * (4.0 / np.pi)
    )

    flag_vacuum = 0  # 1 for in air; 0 for in vacuum

    # 1: metal film coating
    sigma_1 = 18.0
    L_1 = 54e-9

    # 2: thin film
    L_2 = 200e-9
    sigma_2_r = 0.17
    sigma_2_z = 0.17
    capac_2 = 1.4e6
    rho_2 = 1.0e3
    E_2 = 10e9
    niu_2 = 0.30
    alphaT_2 = 40e-6

    # 3: substrate (Si)
    sigma_3 = 142.0
    capac_3 = 1.64e6
    rho_3 = 2.33e3
    alphaT_3 = 2.6e-6
    C11_0_3 = 167.4e9
    C12_0_3 = 65.2e9
    C44_0_3 = (C11_0_3 - C12_0_3) / 2.0

    # properties of air at room temperature
    sigma_4 = 0.026
    capac_4 = 1192.0
    Dif_4 = sigma_4 / capac_4
    dndT_4 = -9e-7

    # properties of Nb-V
    C11_0_1 = 242e9
    C12_0_1 = 129e9
    C44_0_1 = 28e9
    rho_1 = 6.1e3
    alphaT_1 = 0e-6
    capac_1 = 2.65e6

    # Preprocessing of parameters
    Dif_1 = sigma_1 / capac_1

    # effective elastic constants for Al with (111) texture
    C11_1 = (C11_0_1 + C12_0_1 + 2 * C44_0_1) / 2.0
    C33_1 = (C11_0_1 + 2 * C12_0_1 + 4 * C44_0_1) / 3.0
    C44_1 = (C11_0_1 - C12_0_1 + C44_0_1) / 3.0
    C12_1 = (C11_0_1 + 5 * C12_0_1 - 2 * C44_0_1) / 6.0
    C13_1 = (C11_0_1 + 2 * C12_0_1 - 2 * C44_0_1) / 3.0
    C46_1 = 0.0
    C22_1 = C11_1
    C23_1 = C13_1
    C55_1 = C44_1
    C66_1 = (C11_1 - C12_1) / 2.0

    beta_1 = (C11_0_1 + 2 * C12_0_1) * alphaT_1
    betax_1 = beta_1
    betay_1 = beta_1
    betaz_1 = beta_1

    C22C11_1 = C22_1 / C11_1
    C33C11_1 = C33_1 / C11_1
    C12C11_1 = C12_1 / C11_1
    C13C11_1 = C13_1 / C11_1
    C23C11_1 = C23_1 / C11_1
    C44C11_1 = C44_1 / C11_1
    C55C11_1 = C55_1 / C11_1
    C66C11_1 = C66_1 / C11_1
    C46C11_1 = C46_1 / C11_1
    sqrtC11rho_1 = np.sqrt((1.0 + (1e-6) * 1j) * C11_1 / rho_1)
    betaxC11_1 = betax_1 / C11_1
    betayC11_1 = betay_1 / C11_1
    betazC11_1 = betaz_1 / C11_1

    # thin film
    Dif_2_z = sigma_2_z / capac_2
    eta_2 = sigma_2_r / sigma_2_z

    C11_0_2 = E_2 * (1 - niu_2) / (1 + niu_2) / (1 - 2 * niu_2)
    C12_0_2 = E_2 * niu_2 / (1 + niu_2) / (1 - 2 * niu_2)
    C44_0_2 = E_2 / 2.0 / (1 + niu_2)
    C12C11_2 = C12_0_2 / C11_0_2
    C44C11_2 = C44_0_2 / C11_0_2
    sqrtC11rho_2 = np.sqrt((1.0 + (1e-4) * 1j) * C11_0_2 / rho_2)
    beta_2 = (C11_0_2 + 2 * C12_0_2) * alphaT_2
    betaC11_2 = beta_2 / C11_0_2

    # substrate
    Dif_3 = sigma_3 / capac_3
    C12C11_3 = C12_0_3 / C11_0_3
    C44C11_3 = C44_0_3 / C11_0_3
    sqrtC11rho_3 = np.sqrt((1.0 + (1e-4) * 1j) * C11_0_3 / rho_3)
    beta_3 = (C11_0_3 + 2 * C12_0_3) * alphaT_3
    betaC11_3 = beta_3 / C11_0_3

    # discretization of reciprocal space p
    n_p = 63
    up_p = 8.0 / r_rms
    d_p = up_p / n_p
    pp = np.arange(1, n_p + 1) * d_p

    n_f = len(f)
    theta_deformation = np.empty(n_f, dtype=complex)
    theta_dndT_air = np.empty(n_f, dtype=complex)

    det_factor = np.sqrt(8.0 / np.pi) * (focal_length / w_1_d)

    # Precompute Simpson integration weights
    simp_weights = np.zeros(n_p)
    simp_weights[0] = 1.0
    simp_weights[-1] = 1.0
    simp_weights[2:n_p - 2:2] = 2.0
    simp_weights[1:n_p - 1:2] = 4.0
    simp_weights *= (d_p / 3.0)

    # Precompute p-dependent arrays
    bes = -scipy.special.jv(1, pp * xoffset)
    exp_term = np.exp(-((r_rms * pp) ** 2) / 8.0) * (pp**2)
    bes_exp = bes * exp_term
    fs_1_p = A_pump * np.exp(-((r_rms * pp) ** 2) / 8.0)

    A_1_all = []
    A_2_all = []
    A_3_all = []
    N_1_all = []
    N_2_all = []
    N_3_all = []
    B_1_base_all = []
    B_2_base_all = []
    B_3_base_all = []

    for i_p in range(n_p):
        p = pp[i_p]
        k = p / np.sqrt(2.0)
        xi = p / np.sqrt(2.0)

        A_1 = np.zeros((6, 6), dtype=complex)
        A_1[0, 3] = 1.0; A_1[1, 4] = 1.0; A_1[2, 5] = 1.0
        A_1[3, 0] = C55C11_1; A_1[4, 1] = C44C11_1; A_1[5, 2] = C33C11_1
        A_1[0, 0] = -C46C11_1 * 1j * k
        A_1[1, 1] = C46C11_1 * 1j * k
        A_1[0, 1] = C46C11_1 * 1j * xi
        A_1[1, 0] = C46C11_1 * 1j * xi
        A_1[0, 2] = C13C11_1 * 1j * k
        A_1[1, 2] = C23C11_1 * 1j * xi
        D_1 = np.zeros(6, dtype=complex)
        D_1[0] = betaxC11_1 * 1j * k
        D_1[1] = betayC11_1 * 1j * xi
        D_1[5] = betazC11_1
        N_1 = np.linalg.solve(A_1, D_1)

        A_2 = np.zeros((6, 6), dtype=complex)
        A_2[0, 3] = 1.0; A_2[1, 4] = 1.0; A_2[2, 5] = 1.0
        A_2[3, 0] = C44C11_2; A_2[4, 1] = C44C11_2; A_2[5, 2] = 1.0
        A_2[0, 2] = C12C11_2 * 1j * k
        A_2[1, 2] = C12C11_2 * 1j * xi
        D_2 = np.zeros(6, dtype=complex)
        D_2[0] = betaC11_2 * 1j * k
        D_2[1] = betaC11_2 * 1j * xi
        D_2[5] = betaC11_2
        N_2 = np.linalg.solve(A_2, D_2)

        A_3 = np.zeros((6, 6), dtype=complex)
        A_3[0, 3] = 1.0; A_3[1, 4] = 1.0; A_3[2, 5] = 1.0
        A_3[3, 0] = C44C11_3; A_3[4, 1] = C44C11_3; A_3[5, 2] = 1.0
        A_3[0, 2] = C12C11_3 * 1j * k
        A_3[1, 2] = C12C11_3 * 1j * xi
        D_3 = np.zeros(6, dtype=complex)
        D_3[0] = betaC11_3 * 1j * k
        D_3[1] = betaC11_3 * 1j * xi
        D_3[5] = betaC11_3
        N_3 = np.linalg.solve(A_3, D_3)

        B_1_base = np.zeros((6, 6), dtype=complex)
        B_1_base[0, 0] = k**2 + C66C11_1 * xi**2
        B_1_base[1, 1] = C22C11_1 * xi**2 + C66C11_1 * k**2
        B_1_base[0, 1] = (C12C11_1 + C66C11_1) * k * xi
        B_1_base[1, 0] = (C12C11_1 + C66C11_1) * k * xi
        B_1_base[0, 2] = C46C11_1 * (xi**2 - k**2)
        B_1_base[1, 2] = 2.0 * C46C11_1 * k * xi
        B_1_base[2, 3] = -1j * k
        B_1_base[2, 4] = -1j * xi
        B_1_base[3, 0] = C46C11_1 * 1j * k
        B_1_base[3, 1] = -C46C11_1 * 1j * xi
        B_1_base[3, 2] = -C55C11_1 * 1j * k
        B_1_base[4, 0] = -C46C11_1 * 1j * xi
        B_1_base[4, 1] = -C46C11_1 * 1j * k
        B_1_base[4, 2] = -C44C11_1 * 1j * xi
        B_1_base[5, 0] = -C13C11_1 * 1j * k
        B_1_base[5, 1] = -C23C11_1 * 1j * xi
        B_1_base[3, 3] = 1.0; B_1_base[4, 4] = 1.0; B_1_base[5, 5] = 1.0

        B_2_base = np.zeros((6, 6), dtype=complex)
        B_2_base[0, 0] = k**2 + C44C11_2 * xi**2
        B_2_base[1, 1] = xi**2 + C44C11_2 * k**2
        B_2_base[0, 1] = (C12C11_2 + C44C11_2) * k * xi
        B_2_base[1, 0] = (C12C11_2 + C44C11_2) * k * xi
        B_2_base[2, 3] = -1j * k
        B_2_base[2, 4] = -1j * xi
        B_2_base[3, 2] = -C44C11_2 * 1j * k
        B_2_base[4, 2] = -C44C11_2 * 1j * xi
        B_2_base[5, 0] = -C12C11_2 * 1j * k
        B_2_base[5, 1] = -C12C11_2 * 1j * xi
        B_2_base[3, 3] = 1.0; B_2_base[4, 4] = 1.0; B_2_base[5, 5] = 1.0

        B_3_base = np.zeros((6, 6), dtype=complex)
        B_3_base[0, 0] = k**2 + C44C11_3 * xi**2
        B_3_base[1, 1] = xi**2 + C44C11_3 * k**2
        B_3_base[0, 1] = (C12C11_3 + C44C11_3) * k * xi
        B_3_base[1, 0] = (C12C11_3 + C44C11_3) * k * xi
        B_3_base[2, 3] = -1j * k
        B_3_base[2, 4] = -1j * xi
        B_3_base[3, 2] = -C44C11_3 * 1j * k
        B_3_base[4, 2] = -C44C11_3 * 1j * xi
        B_3_base[5, 0] = -C12C11_3 * 1j * k
        B_3_base[5, 1] = -C12C11_3 * 1j * xi
        B_3_base[3, 3] = 1.0; B_3_base[4, 4] = 1.0; B_3_base[5, 5] = 1.0

        A_1_all.append(A_1)
        A_2_all.append(A_2)
        A_3_all.append(A_3)
        N_1_all.append(N_1)
        N_2_all.append(N_2)
        N_3_all.append(N_3)
        B_1_base_all.append(B_1_base)
        B_2_base_all.append(B_2_base)
        B_3_base_all.append(B_3_base)

    I_p2 = np.empty(n_p, dtype=complex)
    I_air_p2 = np.empty(n_p, dtype=complex)

    c_p_pi = -C_probe / np.pi
    c11_ratio_21 = C11_0_2 / C11_1
    c11_ratio_32 = C11_0_3 / C11_0_2

    for i_fr in range(n_f):
        omega = 2 * np.pi * f[i_fr]
        qn2_1 = 1j * omega / Dif_1
        qn2_2 = 1j * omega / Dif_2_z
        qn2_3 = 1j * omega / Dif_3
        qn2_4 = 1j * omega / Dif_4

        w2_1 = -omega**2 / sqrtC11rho_1**2
        w2_2 = -omega**2 / sqrtC11rho_2**2
        w2_3 = -omega**2 / sqrtC11rho_3**2

        for i_p in range(n_p):
            p = pp[i_p]
            p_sq = p**2

            fs_1 = fs_1_p[i_p]
            zeta_1 = np.sqrt(qn2_1 + p_sq)
            zeta_2 = np.sqrt(qn2_2 + eta_2 * p_sq)
            zeta_3 = np.sqrt(qn2_3 + p_sq)
            zeta_4 = np.sqrt(qn2_4 + p_sq)

            zeta_1L_1 = zeta_1 * L_1
            zeta_2L_2 = zeta_2 * L_2

            sigma_1zeta_1 = sigma_1 * zeta_1
            sigma_2zeta_2 = sigma_2_z * zeta_2
            sigma_3zeta_3 = sigma_3 * zeta_3
            sigma_4zeta_4 = sigma_4 * zeta_4

            tanh1 = np.tanh(zeta_1L_1)
            tanh2 = np.tanh(zeta_2L_2)

            m1_01 = -tanh1 / sigma_1zeta_1
            m1_10 = -tanh1 * sigma_1zeta_1
            m2_01 = -tanh2 / sigma_2zeta_2
            m2_10 = -tanh2 * sigma_2zeta_2
            m3_01 = -1.0 / sigma_3zeta_3
            m3_10 = -sigma_3zeta_3

            m2m1_00 = 1.0 + m2_01 * m1_10
            m2m1_01 = m1_01 + m2_01
            m2m1_10 = m2_10 + m1_10
            m2m1_11 = m2_10 * m1_01 + 1.0

            m_10 = m3_10 * m2m1_00 + m2m1_10
            m_11 = m3_10 * m2m1_01 + m2m1_11

            G_d = -m_11 / m_10
            if flag_vacuum > 0:
                G_u = 1.0 / sigma_4zeta_4
                G_tot = 1.0 / (1.0 / G_u + 1.0 / G_d)
            else:
                G_tot = G_d

            thetas_1 = fs_1 * G_tot
            fs_1_d = thetas_1 / G_d
            cosh_z1 = np.cosh(zeta_1L_1)
            cosh_z2 = np.cosh(zeta_2L_2)

            thetas_2 = (thetas_1 + m1_01 * fs_1_d) * cosh_z1
            thetas_3 = (m2m1_00 * thetas_1 + m2m1_01 * fs_1_d) * cosh_z1 * cosh_z2

            exp_z1 = np.exp(zeta_1L_1)
            exp_neg_z1 = np.exp(-zeta_1L_1)
            sinh_z1_2 = exp_z1 - exp_neg_z1

            am_1 = (thetas_2 - thetas_1 * exp_neg_z1) / sinh_z1_2
            ap_1 = thetas_1 - am_1

            exp_z2 = np.exp(zeta_2L_2)
            exp_neg_z2 = np.exp(-zeta_2L_2)
            sinh_z2_2 = exp_z2 - exp_neg_z2

            am_2 = (thetas_3 - thetas_2 * exp_neg_z2) / sinh_z2_2
            ap_2 = thetas_2 - am_2

            am_3 = 0.0
            ap_3 = thetas_3

            B_1 = B_1_base_all[i_p].copy()
            B_1[0, 0] += w2_1; B_1[1, 1] += w2_1; B_1[2, 2] = w2_1

            B_2 = B_2_base_all[i_p].copy()
            B_2[0, 0] += w2_2; B_2[1, 1] += w2_2; B_2[2, 2] = w2_2

            B_3 = B_3_base_all[i_p].copy()
            B_3[0, 0] += w2_3; B_3[1, 1] += w2_3; B_3[2, 2] = w2_3

            A_1 = A_1_all[i_p]
            A_2 = A_2_all[i_p]
            A_3 = A_3_all[i_p]

            LAMBDA_1, Q_1 = scipy.linalg.eig(B_1, A_1, check_finite=False)
            w_raw_2, Q_raw_2 = scipy.linalg.eig(B_2, A_2, check_finite=False)
            w_raw_3, Q_raw_3 = scipy.linalg.eig(B_3, A_3, check_finite=False)

            sorted_idx_2 = np.argsort(w_raw_2.real >= 0, kind='stable')
            LAMBDA_2 = w_raw_2[sorted_idx_2]
            Q_2 = Q_raw_2[:, sorted_idx_2]

            sorted_idx_3 = np.argsort(w_raw_3.real >= 0, kind='stable')
            LAMBDA_3 = w_raw_3[sorted_idx_3]
            Q_3 = Q_raw_3[:, sorted_idx_3]

            U_1 = np.linalg.solve(Q_1, N_1_all[i_p])
            U_2 = np.linalg.solve(Q_2, N_2_all[i_p])
            U_3 = np.linalg.solve(Q_3, N_3_all[i_p])

            BCM = np.zeros((15, 15), dtype=complex)
            BCM[0:3, 0:6] = Q_1[3:6, :]
            exp_L1_1 = np.exp(LAMBDA_1 * L_1)
            BCM[3:9, 0:6] = Q_1[0:6, :] * exp_L1_1

            exp_L2_1 = np.exp(LAMBDA_2 * L_1)
            exp_L2_12 = np.exp(LAMBDA_2 * (L_1 + L_2))
            BCM[3:6, 6:12] = -Q_2[0:3, :] * exp_L2_1
            BCM[6:9, 6:12] = -c11_ratio_21 * Q_2[3:6, :] * exp_L2_1
            BCM[9:15, 6:12] = Q_2[0:6, :] * exp_L2_12

            exp_L3_12 = np.exp(LAMBDA_3[:3] * (L_1 + L_2))
            BCM[9:12, 12:15] = -Q_3[0:3, :3] * exp_L3_12
            BCM[12:15, 12:15] = -c11_ratio_32 * Q_3[3:6, :3] * exp_L3_12

            t1 = U_1 * (am_1 / (zeta_1 - LAMBDA_1) + ap_1 / (-zeta_1 - LAMBDA_1))
            t1_ext = U_1 * (am_1 * exp_z1 / (zeta_1 - LAMBDA_1) + ap_1 * exp_neg_z1 / (-zeta_1 - LAMBDA_1))

            t2 = U_2 * (am_2 / (zeta_2 - LAMBDA_2) + ap_2 / (-zeta_2 - LAMBDA_2))
            t2_ext = U_2 * (am_2 * exp_z2 / (zeta_2 - LAMBDA_2) + ap_2 * exp_neg_z2 / (-zeta_2 - LAMBDA_2))

            t3 = U_3[:3] * (am_3 / (zeta_3 - LAMBDA_3[:3]) + ap_3 / (-zeta_3 - LAMBDA_3[:3]))

            BCC = np.zeros(15, dtype=complex)
            BCC[0:3] = -Q_1[3:6, :] @ t1

            sum1_1_a = Q_1[0:3, :] @ t1_ext
            sum2_2_a = Q_2[0:3, :] @ t2
            BCC[3:6] = -sum1_1_a + sum2_2_a

            sum1_1_b = Q_1[3:6, :] @ t1_ext
            sum2_2_b = Q_2[3:6, :] @ t2
            BCC[6:9] = -sum1_1_b + c11_ratio_21 * sum2_2_b

            sum1_2_a = Q_2[0:3, :] @ t2_ext
            sum2_3_a = Q_3[0:3, :3] @ t3
            BCC[9:12] = -sum1_2_a + sum2_3_a

            sum1_2_b = Q_2[3:6, :] @ t2_ext
            sum2_3_b = Q_3[3:6, :3] @ t3
            BCC[12:15] = -sum1_2_b + c11_ratio_32 * sum2_3_b

            J = np.linalg.solve(BCM, BCC)

            w_s_H = np.dot(Q_1[2, :], J[:6])
            w_s_P = np.dot(Q_1[2, :], t1)

            Z_deformation = -(w_s_H + w_s_P)
            Z_dndT = -dndT_4 * thetas_1 / zeta_4

            be = bes_exp[i_p]
            I_p2[i_p] = Z_deformation * be
            I_air_p2[i_p] = Z_dndT * be

        theta_deformation[i_fr] = c_p_pi * np.dot(simp_weights, I_p2)
        theta_dndT_air[i_fr] = c_p_pi * np.dot(simp_weights, I_air_p2)

    if flag_vacuum > 0:
        theta = theta_deformation + theta_dndT_air
    else:
        theta = theta_deformation

    # read data for comparison
    V_out_data, V_in_data, V_ratio_data, V_SUM_data, fd = (
        GetData_out_in_ratio_f_VSUM(FileNames_data)
    )
    V_SUM = np.mean(V_SUM_data) * 4.0

    # leaking correction parameters
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

    # Plot 1: Probe beam deflection angle theta
    complex_theta_to_plot = theta

    plt.figure(1, figsize=(10, 4))
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

    plt.figure(2, figsize=(10, 4))
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
    plt.loglog(
        f,
        -np.real(complex_theta_to_plot) / np.imag(complex_theta_to_plot),
        'k-',
        linewidth=1.5,
        label='Sim Ratio',
    )
    plt.loglog(fd, Vcorrected_ratio, 'b-', linewidth=1.5, label='Exp Ratio')
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
    }


if __name__ == '__main__':
    results = run_simulation('a1')
    plt.show()
