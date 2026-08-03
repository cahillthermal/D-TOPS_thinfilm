%%Spilt part for sensitivity calculation of thin film TOPS
function[theta]=SensThinfilm(sigma_1,capac_1,C11_0_1,C12_0_1,C44_0_1,alphaT_1,rho_1,sigma_2_z,capac_2,sigma_2_r,niu_2,E_2,rho_2,alphaT_2,sigma_3,capac_3,C12_0_3,C11_0_3,C44_0_3,alphaT_3,rho_3,r_rms,f,focal_length,w_1_d,det_factor,V_SUM,flag_vacuum,A_pump,L_1,L_2,xoffset,C_probe)


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
%air
sigma_4 = 0.028;  % thermal conductivity (W/m-K)
capac_4 = 1192; % volumetric heat capacity (J/m^3-K)
Dif_4 = sigma_4/capac_4;
dndT_4 = -9e-7; % thermo-optic coefficient dn/dT (K^(-1))
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
        inv_A_1 = inv(A_1);
        N_1 = inv_A_1*D_1;
        [Q_1, R_1] = eig(B_1,A_1);
        LAMBDA_1 = zeros(6,1);
        for i = 1:6
            LAMBDA_1(i) = R_1(i,i);
        end
        condi_1 = rcond(Q_1);
        if condi_1 < 2e-16
            pause;
        end
        inv_Q_1 = inv(Q_1);
        U_1 = inv_Q_1*N_1;

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
        inv_A_2 = inv(A_2);
        N_2 = inv_A_2*D_2;
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
        inv_Q_2 = inv(Q_2);
        U_2 = inv_Q_2*N_2;

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
        inv_A_3 = inv(A_3);
        N_3 = inv_A_3*D_3;
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
        inv_Q_3 = inv(Q_3);
        U_3 = inv_Q_3*N_3;

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

% conversion from PBD angle to signal of lock-in amplifier
V_signal_deformation(:) = det_factor*theta_deformation(:)*V_SUM/sqrt(2);
V_signal_dndT_air(:) = det_factor*theta_dndT_air(:)*V_SUM/sqrt(2);
V_signal(:) = det_factor*theta(:)*V_SUM/sqrt(2);

end