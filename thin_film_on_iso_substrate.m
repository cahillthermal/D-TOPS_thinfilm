%clear variables;
close all;
%% READ ME!

% structure from top to bottom:
% air/Al coating/thin film/substrate

% the index of layer used in this code
% 4: air (or vacuum); 1: Al coating; 2: thin film ; 3: substrate

% model assumptions:
% 1) thin film has isotropic elastic and thermal expansion properties
% 2) thin film has in-plane thermal conductivity sigma_r and cross-plane thermal conductivity sigma_z
% 3) substrate has isotropic elastic and thermal expansion properties
% 4) substrate has isotropic isotropic thermal conductivity

% result of this code is probe beam deflection angle or signal due to:
% 1) thermal expansion induced surface deformation of the sample surface
% 2) in-plane gradient of refractive index in the air due to dn/dT

%% parameters:
FileNames_data = 'apr2226/alon_2x';   % data to compare with
ii=sqrt(-1);
f = logspace(log10(2e3),log10(200e3),40)'; % modulation frequency (Hz)

% lens
lens_magnification=10;                   % e.g., 5 for 5x lens; 10 for 10x lens
lens_transmittance = 0.80; % transmittance of objective lens for pump beam, 0.80 for 2x; 0.88 for 5x; 0.85 for 10x; 0.75 for 20x
focal_length = 5/lens_magnification*40e-3; % focal length of the objective lens

% optical beam spot size and beam offset
r_rms=5/lens_magnification*12.8e-6;     % root-mean-square focused pump and probe beam 1/e^2 radius (m)
xoffset=5/lens_magnification*13.5e-6;   % Beam offset (m); pump beam moved laterally by 9.2 um move of the actuator on the gimbal mount
C_probe = 0.90;                         % 0.90 approximately for TOPS 2.0.
w_1_d = 0.90e-3;                        % probe beam 1/e^2 radius at the detector, 0.90 for TOPS 2.0

% set pump beam power
incident_pump=4e-3;        % avarage power of square wave modulated pump at the bfp of objective (W)

% absorbance of pump power
% Index of refraction of metal coating at the wavelength of the pump beam
n_metal=2.9; % for Al at 780 nm (from Mathewson and Myers):
k_metal=8.2;
sample_reflectance=abs(n_metal-1+(1i)*k_metal)^2/abs(n_metal+1+(1i)*k_metal)^2; % reflectance of the sample surface
%sample_absorbance=1-sample_reflectance;  % absorbance of sample surface
sample_absorbance=0.4;  % 0.4 is for NbV; index of NbV is n=2.6; k=3.6
A_pump=incident_pump*lens_transmittance*sample_absorbance*(4.0/pi); % Amplitude of the fundamental Fourier component of the absorbed pump beam

flag_vacuum = 1; % 1 for in air; 0 for in vacuum

% 1: Al, NbV, Ni, or other metal film coating
sigma_1 = 20; % thermal conductivity of metal film; Al typically 150 W/(m K); NbV 20 W/(m K)
L_1 = 60e-9; % thickness of metal

% % 2: thin film example for elastically isotropic material with Poisson ratio of 0.3, ratio of C11, C12, C44 is 7:3:2
L_2 = 3e-6; % thickness
sigma_2_r = 0.3; % thermal conductivity
sigma_2_z = 0.2;
capac_2 = 2.36e6; % volumetric heat capacity
rho_2 = 1.0e3; % mass density doesn't play a role because sound velocity is high but is included in the calculation
E_2 = 10e9; % Young's modulus
niu_2 = 0.30; % Poisson's ratio
alphaT_2 = 100e-6; % linear coefficient of thermal expansion

% % 3: substrate (fused silica is used as example here)
% sigma_3 = 1.3;
% capac_3 = 1.63e6;
% rho_3 = 2.20e3;
% alphaT_3 = 0.55e-6; % linear coefficient of thermal expansion
% % the key point below is to set values for the three elastic constants
% % for amorphous materials, given Young's modulus and Poisson's ratio
% % the formula below can be used:
% E_3 = 72e9; % Young's modulus
% niu_3 = 0.17; % Poisson's ratio
% C11_0_3 = E_3*(1-niu_3)/(1+niu_3)/(1-2*niu_3);
% C12_0_3 = E_3*niu_3/(1+niu_3)/(1-2*niu_3);
% C44_0_3 = E_3/2/(1+niu_3);

% 3: substrate (Si approximated as being isotropic is used as example here)
sigma_3 = 140;  % thermal conductivity (W/m-K)
capac_3 = 1.64e6; % volumetric heat capacity (J/m^3-K)
rho_3 = 2.33e3; % mass density (kg/m^3)
alphaT_3 = 2.6e-6; % linear coefficient of thermal expansion
% % the key point below is to set values for the three elastic constants
% for crystal materials, generally the elastic property is not isotropic
% and thus this code is not rigorous
% still, it turns out to be fine to approximate a (100) substrate as being
% isotropic by setting C44_0_3 = (C11_0_3-C12_0_3)/2;
% below is as an example
% values of C11_0_3 and C12_0_3 are looked up from literature
% while that of C44_0_3 is calculated using this formula
C11_0_3 = 167.4e9; % elastic constant C11 (Pa)
C12_0_3 = 65.2e9; % elastic constant C12 (Pa)
C44_0_3 = (C11_0_3-C12_0_3)/2; % elastic constant C44 (Pa);

%% parameters by default

% properties of air at room temperature
sigma_4 = 0.026;  % thermal conductivity (W/m-K)
capac_4 = 1192; % volumetric heat capacity (J/m^3-K)
Dif_4 = sigma_4/capac_4;
dndT_4 = -9e-7; % thermo-optic coefficient dn/dT (K^(-1))

% properties of Al
%C11_0_1 = 107.4e9; % elastic constant C11 (Pa)
%C12_0_1 = 60.5e9; % elastic constant C12 (Pa)
%C44_0_1 = 28.3e9; % elastic constant C44 (Pa)
%rho_1 = 2.70e3; % mass density (kg/m^3)
%alphaT_1 = 23e-6; % linear coefficient of thermal expansion (K^(-1))
%capac_1 = 2.42e6; % volumetric heat capacity (J/m^3-K)

% properties of Nb-V
C11_0_1 = 242e9; % elastic constant C11 (Pa)
C12_0_1 = 129e9; % elastic constant C12 (Pa)
C44_0_1 = 28e9; % elastic constant C44 (Pa)
rho_1 = 6.1e3; % mass denstiy (kg/m^3)
alphaT_1 = 7.9e-6; % linear coefficient of thermal expansion (K^(-1)), Nb 7.3e-6/V 8.4e-6
capac_1 = 2.65e6; % volumetric heat capacity (J/m^3-K)
%coef_intf_rfl = 0;
% coef_intf_rfl is defined by coef_intf_rfl = -lambda_1/(4*pi)*dphase(r)/dT
% where lambda_1 is wavelength of probe beam
% r is the reflection coefficient
% for optically opaque Al film (thickness above 30 nm) in air:
% coef_intf_rfl = 0 was found experimentally
%% preprocessing of parameters

% for Al coating:

Dif_1 = sigma_1/capac_1;

% effective elastic constants for Al with (111) texture:
C11_1 = (C11_0_1 + C12_0_1 + 2*C44_0_1)/2;
C33_1 = (C11_0_1 + 2*C12_0_1 + 4*C44_0_1)/3;
C44_1 = (C11_0_1 - C12_0_1 + C44_0_1)/3;
C12_1 = (C11_0_1 + 5*C12_0_1 - 2*C44_0_1)/6;
C13_1 = (C11_0_1 + 2*C12_0_1 - 2*C44_0_1)/3;
C46_1 = 0;
C22_1 = C11_1;
C23_1 = C13_1;
C55_1 = C44_1;
C66_1 = (C11_1 - C12_1)/2;

beta_1 = (C11_0_1 + 2*C12_0_1)*alphaT_1;
betax_1 = beta_1;
betay_1 = beta_1;
betaz_1 = beta_1;

C22C11_1 = C22_1/C11_1;
C33C11_1 = C33_1/C11_1;
C12C11_1 = C12_1/C11_1;
C13C11_1 = C13_1/C11_1;
C23C11_1 = C23_1/C11_1;
C44C11_1 = C44_1/C11_1;
C55C11_1 = C55_1/C11_1;
C66C11_1 = C66_1/C11_1;
C46C11_1 = C46_1/C11_1;
sqrtC11rho_1 = sqrt((1+(1e-6)*(1i))*C11_1/rho_1);
betaxC11_1 = betax_1/C11_1;
betayC11_1 = betay_1/C11_1;
betazC11_1 = betaz_1/C11_1;


% for thin film

Dif_2_z = sigma_2_z/capac_2;
eta_2 = sigma_2_r/sigma_2_z;

C11_0_2 = E_2*(1-niu_2)/(1+niu_2)/(1-2*niu_2);
C12_0_2 = E_2*niu_2/(1+niu_2)/(1-2*niu_2);
C44_0_2 = E_2/2/(1+niu_2);
C12C11_2 = C12_0_2/C11_0_2;
C44C11_2 = C44_0_2/C11_0_2;
sqrtC11rho_2 = sqrt((1+(1e-4)*(1i))*C11_0_2/rho_2);
beta_2 = (C11_0_2 + 2*C12_0_2)*alphaT_2;
betaC11_2 = beta_2/C11_0_2;

% for substrate:

Dif_3 = sigma_3/capac_3;

C12C11_3 = C12_0_3/C11_0_3;
C44C11_3 = C44_0_3/C11_0_3;
sqrtC11rho_3 = sqrt((1+(1e-4)*(1i))*C11_0_3/rho_3);
beta_3 = (C11_0_3 + 2*C12_0_3)*alphaT_3;
betaC11_3 = beta_3/C11_0_3;

%% discretization of reciprocal space p
% real space to reciprocal space is by Hankel transform: r -> p

n_p = 63;
up_p = 8/r_rms;
d_p = up_p/n_p;
pp = d_p:d_p:up_p;

%% defining the through-plane displacement after transform (outward: +)
Z_deformation_p_omega = zeros(n_p,length(f));
% displacement due to thermoelastic deformation of the sample surface

Z_dndT_p_omega = zeros(n_p,length(f));
% effective displacement due to dn/dT of air

%% defining corresponding PBD angle, theta
theta_deformation = ones(length(f),1);
theta_dndT_air = ones(length(f),1);
theta = ones(length(f),1); % theta is the sum of theta_deformation and theta_dndT_air

%% defining corresponding lock-in amplifier signal, V_signal
V_signal_deformation = ones(length(f),1);
V_signal_dndT_air = ones(length(f),1);
V_signal = ones(length(f),1);

% relation between PBD angle theta and V_signal
det_factor=(8/pi)^0.5*(focal_length/w_1_d);
% det_factor is the factor such that:
% V_signal*sqrt(2)/V_SUM = det_factor*theta

%% computation

% defining matrix A, B and vector D for layers 1,2,3

A_1 = zeros(6,6);
B_1 = zeros(6,6);
D_1 = zeros(6,1);
A_1(1,4) = 1;
A_1(2,5) = 1;
A_1(3,6) = 1;
A_1(4,1) = C55C11_1;
A_1(5,2) = C44C11_1;
A_1(6,3) = C33C11_1;
B_1(4,4) = 1;
B_1(5,5) = 1;
B_1(6,6) = 1;
D_1(6) = betazC11_1;

A_2 = zeros(6,6);
B_2 = zeros(6,6);
D_2 = zeros(6,1);
A_2(1,4) = 1;
A_2(2,5) = 1;
A_2(3,6) = 1;
A_2(4,1) = C44C11_2;
A_2(5,2) = C44C11_2;
A_2(6,3) = 1;
B_2(4,4) = 1;
B_2(5,5) = 1;
B_2(6,6) = 1;
D_2(6) = betaC11_2;

A_3 = zeros(6,6);
B_3 = zeros(6,6);
D_3 = zeros(6,1);
A_3(1,4) = 1;
A_3(2,5) = 1;
A_3(3,6) = 1;
A_3(4,1) = C44C11_3;
A_3(5,2) = C44C11_3;
A_3(6,3) = 1;
B_3(4,4) = 1;
B_3(5,5) = 1;
B_3(6,6) = 1;
D_3(6) = betaC11_3;

for i_fr = 1:length(f)
    omega = 2*pi*f(i_fr);
    qn2_1 = (1i)*omega/Dif_1;
    qn2_2 = (1i)*omega/Dif_2_z;
    qn2_3 = (1i)*omega/Dif_3;
    qn2_4 = (1i)*omega/Dif_4;
    I_p2 = zeros(n_p,1);
    I_intf_rfl_p2 = zeros(n_p,1);
    I_air_p2 = zeros(n_p,1);

    for i_p = 1:n_p
        p = pp(i_p);

        % computation of interface temperature
        fs_1 = A_pump*exp(-r_rms^2*p^2/8);
        zeta_1 = sqrt(qn2_1 + p^2);
        zeta_2 = sqrt(qn2_2 + eta_2*p^2);
        zeta_3 = sqrt(qn2_3 + p^2);
        zeta_4 = sqrt(qn2_4 + p^2);
        zeta_1L_1 = zeta_1*L_1;
        zeta_2L_2 = zeta_2*L_2;
        sigma_1zeta_1 = sigma_1*zeta_1;
        sigma_2zeta_2 = sigma_2_z*zeta_2;
        sigma_3zeta_3 = sigma_3*zeta_3;
        sigma_4zeta_4 = sigma_4*zeta_4;
        M1 = zeros(2,2);
        M2 = zeros(2,2);
        M3 = zeros(2,2);
        M3(1,1) = 1;
        M3(1,2) = -1/sigma_3zeta_3;
        M3(2,1) = -sigma_3zeta_3;
        M3(2,2) = 1;
        M2(1,1) = 1;
        M2(1,2) = -tanh(zeta_2L_2)/sigma_2zeta_2;
        M2(2,1) = -tanh(zeta_2L_2)*sigma_2zeta_2;
        M2(2,2) = 1;
        M1(1,1) = 1;
        M1(1,2) = -tanh(zeta_1L_1)/sigma_1zeta_1;
        M1(2,1) = -tanh(zeta_1L_1)*sigma_1zeta_1;
        M1(2,2) = 1;
        
        M2M1 = M2*M1;
        M = M3*M2M1;
        G_d = -M(2,2)/M(2,1);
        G_u = 1/sigma_4zeta_4;
        if flag_vacuum>0
            G_tot = 1/(1/G_u + 1/G_d);
        else
            G_tot = G_d;
        end
        thetas_1 = fs_1*G_tot; % at air/Al
        fs_1_d = thetas_1/G_d;
        thetas_2 = (M1(1,1)*thetas_1 + M1(1,2)*fs_1_d)*cosh(zeta_1L_1); % at Al/ thin film
        thetas_3 = (M2M1(1,1)*thetas_1 + M2M1(1,2)*fs_1_d)*cosh(zeta_1L_1)*cosh(zeta_2L_2); % at thin film/ substrate

        % computation of the amplitude of thermal wave:
        am_1 = (thetas_2 - thetas_1*exp(-zeta_1L_1))/(exp(zeta_1L_1) - exp(-zeta_1L_1));
        ap_1 = thetas_1 - am_1;
        am_2 = (thetas_3 - thetas_2*exp(-zeta_2L_2))/(exp(zeta_2L_2) - exp(-zeta_2L_2));
        ap_2 = thetas_2 - am_2;
        am_3 = 0;
        ap_3 = thetas_3 - am_3;

        k = p/sqrt(2);
        xi = p/sqrt(2);

        % general solution of stress-stain of Al coating:
        A_1(1,1) = -C46C11_1*(1i)*k;
        A_1(2,2) = C46C11_1*(1i)*k;
        A_1(1,2) = C46C11_1*(1i)*xi;
        A_1(2,1) = C46C11_1*(1i)*xi;
        A_1(1,3) = C13C11_1*(1i)*k;
        A_1(2,3) = C23C11_1*(1i)*xi;
        B_1(1,1) = k^2 + C66C11_1*xi^2 - omega^2/sqrtC11rho_1^2;
        B_1(2,2) = C22C11_1*xi^2 + C66C11_1*k^2 - omega^2/sqrtC11rho_1^2;
        B_1(1,2) = (C12C11_1 + C66C11_1)*k*xi;
        B_1(2,1) = (C12C11_1 + C66C11_1)*k*xi;
        B_1(3,3) = -omega^2/sqrtC11rho_1^2;
        B_1(1,3) = C46C11_1*(xi^2 - k^2);
        B_1(2,3) = 2*C46C11_1*k*xi;
        B_1(3,4) = -(1i)*k;
        B_1(3,5) = -(1i)*xi;
        B_1(4,1) = C46C11_1*(1i)*k;
        B_1(4,2) = -C46C11_1*(1i)*xi;
        B_1(4,3) = -C55C11_1*(1i)*k;
        B_1(5,1) = -C46C11_1*(1i)*xi;
        B_1(5,2) = -C46C11_1*(1i)*k;
        B_1(5,3) = -C44C11_1*(1i)*xi;
        B_1(6,1) = -C13C11_1*(1i)*k;
        B_1(6,2) = -C23C11_1*(1i)*xi;
        D_1(1) = betaxC11_1*(1i)*k;
        D_1(2) = betayC11_1*(1i)*xi;
%        these two lines in Jinchi's code are replaced
%        by the following line suggested by Matlab is faster and more
%        accurate
%        inv_A_1 = inv(A_1);  
%        N_1 = inv_A_1*D_1;
        N_1 = A_1\D_1;
        [Q_1, R_1] = eig(B_1,A_1);
        LAMBDA_1 = zeros(6,1);
        for i = 1:6
            LAMBDA_1(i) = R_1(i,i);
        end
        condi_1 = rcond(Q_1);
        if condi_1 < 2e-16
            pause;
        end
%        inv_Q_1 = inv(Q_1);
        U_1 = Q_1\N_1;

        % general solution of stress-stain of thin film:
        A_2(1,3) = C12C11_2*(1i)*k;
        A_2(2,3) = C12C11_2*(1i)*xi;
        B_2(1,1) = k^2 + C44C11_2*xi^2 - omega^2/sqrtC11rho_2^2;
        B_2(2,2) = xi^2 + C44C11_2*k^2 - omega^2/sqrtC11rho_2^2;
        B_2(1,2) = (C12C11_2 + C44C11_2)*k*xi;
        B_2(2,1) = (C12C11_2 + C44C11_2)*k*xi;
        B_2(3,3) = -omega^2/sqrtC11rho_2^2;
        B_2(3,4) = -(1i)*k;
        B_2(3,5) = -(1i)*xi;
        B_2(4,3) = -C44C11_2*(1i)*k;
        B_2(5,3) = -C44C11_2*(1i)*xi;
        B_2(6,1) = -C12C11_2*(1i)*k;
        B_2(6,2) = -C12C11_2*(1i)*xi;
        D_2(1) = betaC11_2*(1i)*k;
        D_2(2) = betaC11_2*(1i)*xi;
%        inv_A_2 = inv(A_2);
        N_2 = A_2\D_2;
        [Q_raw_2, R_raw_2] = eig(B_2,A_2);
        Q_2 = zeros(6,6);
        R_2 = zeros(6,6);
        count_2 = 0;
        for i = 1:6
            if real(R_raw_2(i,i)) < 0
                count_2 = count_2 + 1;
                R_2(count_2,count_2) = R_raw_2(i,i);
                Q_2(:,count_2) = Q_raw_2(:,i);
            else
                R_2(i-count_2+3,i-count_2+3) = R_raw_2(i,i);
                Q_2(:,i-count_2+3) = Q_raw_2(:,i);
            end
        end
        LAMBDA_2 = zeros(6,1);
        for i = 1:6
            LAMBDA_2(i) = R_2(i,i);
        end
        condi_2 = rcond(Q_2);
        if condi_2 < 2e-16
            pause;
        end
  %      inv_Q_2 = inv(Q_2);
        U_2 = Q_2\N_2;

        % general solution of stress-stain of substrate:
        A_3(1,3) = C12C11_3*(1i)*k;
        A_3(2,3) = C12C11_3*(1i)*xi;
        B_3(1,1) = k^2 + C44C11_3*xi^2 - omega^2/sqrtC11rho_3^2;
        B_3(2,2) = xi^2 + C44C11_3*k^2 - omega^2/sqrtC11rho_3^2;
        B_3(1,2) = (C12C11_3 + C44C11_3)*k*xi;
        B_3(2,1) = (C12C11_3 + C44C11_3)*k*xi;
        B_3(3,3) = -omega^2/sqrtC11rho_3^2;
        B_3(3,4) = -(1i)*k;
        B_3(3,5) = -(1i)*xi;
        B_3(4,3) = -C44C11_3*(1i)*k;
        B_3(5,3) = -C44C11_3*(1i)*xi;
        B_3(6,1) = -C12C11_3*(1i)*k;
        B_3(6,2) = -C12C11_3*(1i)*xi;
        D_3(1) = betaC11_3*(1i)*k;
        D_3(2) = betaC11_3*(1i)*xi;
%        inv_A_3 = inv(A_3);
        N_3 = A_3\D_3;
        [Q_raw_3, R_raw_3] = eig(B_3,A_3);
        Q_3 = zeros(6,6);
        R_3 = zeros(6,6);
        count_3 = 0;
        for i = 1:6
            if real(R_raw_3(i,i)) < 0
                count_3 = count_3 + 1;
                R_3(count_3,count_3) = R_raw_3(i,i);
                Q_3(:,count_3) = Q_raw_3(:,i);
            else
                R_3(i-count_3+3,i-count_3+3) = R_raw_3(i,i);
                Q_3(:,i-count_3+3) = Q_raw_3(:,i);
            end
        end
        LAMBDA_3 = zeros(6,1);
        for i = 1:6
            LAMBDA_3(i) = R_3(i,i);
        end
        condi_3 = rcond(Q_3);
        if condi_3 < 2e-16
            pause;
        end
 %       inv_Q_3 = inv(Q_3);
        U_3 = Q_3\N_3;

        % unique solution determined by boundary and continuous conditions:
        BCM = zeros(15,15);
        BCC = zeros(15,1);
        for i_cl = 1:6
            BCM(1:3,i_cl) = Q_1(4:6,i_cl);
            BCM(4:9,i_cl) = Q_1(1:6,i_cl)*exp(LAMBDA_1(i_cl)*L_1);
        end
        for i_cl = 7:12
            BCM(4:6,i_cl) = -Q_2(1:3,i_cl-6)*exp(LAMBDA_2(i_cl-6)*L_1);
            BCM(7:9,i_cl) = -C11_0_2/C11_1*Q_2(4:6,i_cl-6)*exp(LAMBDA_2(i_cl-6)*L_1);
            BCM(10:15,i_cl) = Q_2(1:6,i_cl-6)*exp(LAMBDA_2(i_cl-6)*(L_1+L_2));
        end
        for i_cl = 13:15
            BCM(10:12,i_cl) = -Q_3(1:3,i_cl-12)*exp(LAMBDA_3(i_cl-12)*(L_1+L_2));
            BCM(13:15,i_cl) = -C11_0_3/C11_0_2*Q_3(4:6,i_cl-12)*exp(LAMBDA_3(i_cl-12)*(L_1+L_2));
        end
        for i_rw = 1:3
            sum = 0;
            for jj = 1:6
                temp = Q_1(i_rw+3,jj)*U_1(jj)*(am_1/(zeta_1-LAMBDA_1(jj)) + ap_1/(-zeta_1-LAMBDA_1(jj)));
                sum = sum + temp;
            end
            BCC(i_rw) = -sum;
        end
        for i_rw = 4:6
            sum1 = 0;
            sum2 = 0;
            for jj = 1:6
                temp1 = Q_1(i_rw-3,jj)*U_1(jj)*(am_1*exp(zeta_1L_1)/(zeta_1-LAMBDA_1(jj)) + ap_1*exp(-zeta_1L_1)/(-zeta_1-LAMBDA_1(jj)));
                temp2 = Q_2(i_rw-3,jj)*U_2(jj)*(am_2/(zeta_2-LAMBDA_2(jj)) + ap_2/(-zeta_2-LAMBDA_2(jj)));
                sum1 = sum1 + temp1;
                sum2 = sum2 + temp2;
            end
            BCC(i_rw) = -sum1 + sum2;
        end
        for i_rw = 7:9
            sum1 = 0;
            sum2 = 0;
            for jj = 1:6
                temp1 = Q_1(i_rw-3,jj)*U_1(jj)*(am_1*exp(zeta_1L_1)/(zeta_1-LAMBDA_1(jj)) + ap_1*exp(-zeta_1L_1)/(-zeta_1-LAMBDA_1(jj)));
                temp2 = Q_2(i_rw-3,jj)*U_2(jj)*(am_2/(zeta_2-LAMBDA_2(jj)) + ap_2/(-zeta_2-LAMBDA_2(jj)));
                sum1 = sum1 + temp1;
                sum2 = sum2 + temp2;
            end
            BCC(i_rw) = -sum1 + (C11_0_2/C11_1)*sum2;
        end
        for i_rw = 10:12
            sum1 = 0;
            sum2 = 0;
            for jj = 1:6
                temp1 = Q_2(i_rw-9,jj)*U_2(jj)*(am_2*exp(zeta_2L_2)/(zeta_2-LAMBDA_2(jj)) + ap_2*exp(-zeta_2L_2)/(-zeta_2-LAMBDA_2(jj)));
                temp2 = Q_3(i_rw-9,jj)*U_3(jj)*(am_3/(zeta_3-LAMBDA_3(jj)) + ap_3/(-zeta_3-LAMBDA_3(jj)));
                sum1 = sum1 + temp1;
                sum2 = sum2 + temp2;
            end
            BCC(i_rw) = -sum1 + sum2;
        end
        for i_rw = 13:15
            sum1 = 0;
            sum2 = 0;
            for jj = 1:6
                temp1 = Q_2(i_rw-9,jj)*U_2(jj)*(am_2*exp(zeta_2L_2)/(zeta_2-LAMBDA_2(jj)) + ap_2*exp(-zeta_2L_2)/(-zeta_2-LAMBDA_2(jj)));
                temp2 = Q_3(i_rw-9,jj)*U_3(jj)*(am_3/(zeta_3-LAMBDA_3(jj)) + ap_3/(-zeta_3-LAMBDA_3(jj)));
                sum1 = sum1 + temp1;
                sum2 = sum2 + temp2;
            end
            BCC(i_rw) = -sum1 + (C11_0_3/C11_0_2)*sum2;
        end

        J = BCM\BCC;

        % surface thermoelastic displacement:
        w_s_H = Q_1(3,1)*J(1) + Q_1(3,2)*J(2) + Q_1(3,3)*J(3) + Q_1(3,4)*J(4) + Q_1(3,5)*J(5) + Q_1(3,6)*J(6);
        sssum = 0;
        for jj = 1:6
            temp = Q_1(3,jj)*U_1(jj)*(am_1/(zeta_1-LAMBDA_1(jj)) + ap_1/(-zeta_1-LAMBDA_1(jj)));
            sssum = sssum + temp;
        end
        w_s_P = sssum;
        Z_deformation_p_omega(i_p,i_fr) = -(w_s_H + w_s_P);
        Z_dndT_p_omega(i_p,i_fr) = -dndT_4*thetas_1/zeta_4;

        % integrand for calculating PBD
        I_p2(i_p) = Z_deformation_p_omega(i_p,i_fr)*(-besselj(1,p*xoffset))*exp(-r_rms^2*p^2/8)*p^2;
        I_air_p2(i_p) = Z_dndT_p_omega(i_p,i_fr)*(-besselj(1,p*xoffset))*exp(-r_rms^2*p^2/8)*p^2;
    end
    theta_deformation(i_fr,1) = -C_probe/pi*simpson_inte(I_p2,d_p);
    theta_dndT_air(i_fr,1) = -C_probe/pi*simpson_inte(I_air_p2,d_p);
    % minus sign for the experiment in which pump beam is below probe beam
end

% summing up all soures of PBD angle
if flag_vacuum>0
    theta(:) = theta_deformation(:)+theta_dndT_air(:);
else
    theta(:) = theta_deformation(:);
end

%% read data for comparison
% load the arrays for the out-of-phase signal, in-phase signal,
% ratio (-in-phase/out-of-phase), frequency, and the detector SUM voltage


[V_out_data,V_in_data,V_ratio_data,V_SUM_data,fd]=GetData_out_in_ratio_f_VSUM(FileNames_data);
V_SUM = mean(V_SUM_data)*4.0;

% calculate the "leaking" data for correction of the frequency response due
% to the imperfection of pump modulation and detector response; these values are for TOPS 2.0
Amplitude_corrected_3 = -6.08e-09; % 3rd order
Amplitude_corrected_2 = 1.18e-06; % 2nd order
Amplitude_corrected_1 = 1.50e-04; % 1st order
Amplitude_corrected_0 = 1 ; % constant

delay_2 = 7.33e-12; % 2nd order
delay_1 = -1.10e-05; % 1st order
delay_0 = 5.08e-03;  % const
complex_leaking = (Amplitude_corrected_0 + Amplitude_corrected_1 * sqrt(fd) + Amplitude_corrected_2 * sqrt(fd).^2+ Amplitude_corrected_3 * sqrt(fd).^3).* (exp(1i*(delay_0 + delay_1*fd + delay_2*fd.^2)));

% correct the measured data using the "leaking" data
[Vcorrected_in,Vcorrected_out,Vcorrected_ratio]=datacorrection_complex_leaking(V_out_data,V_in_data,complex_leaking);


%% conversion from PBD angle to signal of lock-in amplifier
V_signal_deformation(:) = det_factor*theta_deformation(:)*V_SUM/sqrt(2);
V_signal_dndT_air(:) = det_factor*theta_dndT_air(:)*V_SUM/sqrt(2);
V_signal(:) = det_factor*theta(:)*V_SUM/sqrt(2);



%% plot the probe beam deflection angle, theta
complex_theta_to_plot = theta;
% put whatever that needs plotting: theta, theta_deformation, or, theta_dndT_air

figure(1)
subplot(1,2,1)
semilogx(f,1e6*real(complex_theta_to_plot),'k-','linewidth',1.5); hold on
semilogx(f,1e6*imag(complex_theta_to_plot),'k--','linewidth',1.5); hold on
box on; axis tight;
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Arial');
xlabel('f (Hz)');
ylabel('in, out-of-phase (\murad)')

subplot(1,2,2)
loglog(f, -real(complex_theta_to_plot)./imag(complex_theta_to_plot),'k-','linewidth',1.5); hold on
box on; axis tight;
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Arial');
xlabel('f (Hz)');
ylabel('ratio')

%% plot the lock-in amplifier signal, V_signal
complex_V_signal_to_plot = V_signal;
% put whatever that needs plotting: V_signal, V_signal_deformation, or, V_signal_dndT_air

figure(2)
subplot(1,2,1)
semilogx(f,1e6*real(complex_V_signal_to_plot),'k-','linewidth',1.5); hold on
semilogx(f,1e6*imag(complex_V_signal_to_plot),'k--','linewidth',1.5); hold on
semilogx(fd,1e6*Vcorrected_in,'b-','linewidth',1.5); hold on
semilogx(fd,1e6*Vcorrected_out,'b-','linewidth',1.5); hold on
box on; axis tight;
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Arial');
xlabel('f (Hz)');
ylabel('in, out-of-phase (\muV)')
subplot(1,2,2)
loglog(f, -real(complex_theta_to_plot)./imag(complex_theta_to_plot),'k-','linewidth',1.5); hold on
loglog(fd, Vcorrected_ratio,'b-','linewidth',1.5); hold on
box on; axis tight;
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Arial');
xlabel('f (Hz)');
ylabel('ratio')



%% function needed in the model

function I = simpson_inte(array,pace)
steps = length(array);
edge_sum = sum(array(3:2:steps-2));
mid_sum = sum(array(2:2:steps-1));
I = (1/6)*(2*pace)*(array(1) + 2*edge_sum + 4*mid_sum + array(steps));
end