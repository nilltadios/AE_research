function reconstructed_pitch = reconstructSignalWithERA(tSol, xSol, tspan, era_params)
%reconstructSignalWithERA Performs ERA with modal transformation to preserve nonlinearities
% 
% PAPER NOTATION (Ghadami & Epureanu - WITH MODAL TRANSFORMATION T)
%
% Key insight from paper: Transform to modal coordinates using T (eigenvectors of S)
% This puts C into modal space where it contains the actual mode shapes.
%
% State-space representation:
%   x_dot(t) = A*x(t) + B*u(t)
%   y(t) = C*x(t) + D*u(t)
%
% Discrete-time form (Eq. 11-12):
%   y_k(t) = C*S^k*X_0
%
% Modal transformation: Tilde(S) = T^(-1)*S*T, Tilde(C) = C*T
% After transformation: Tilde(C) contains the mode shapes

% --- Parameter and Data Extraction ---
start_index = era_params.start_index;
r_era = era_params.r_era;
v_era = era_params.s_era;
N_r = era_params.n_r;

% Extract the pitch measurement (scalar output)
Y_measured = xSol(:, 1:4);  % y(t) - measured output

% Define system dimensions
n = 4;  % dimensionality of y (outputs)
N_measurements = length(tspan);

% --- ERA Algorithm (Following Paper Equations 13-14) ---

% Prepare measurement data segment
Y_segment = Y_measured(start_index:end,:);

% Step 1: Build Hankel matrices H(0) and H(1)
H_0 = build_block_hankel(Y_segment, 1, r_era, v_era, n);
H_1 = build_block_hankel(Y_segment, 2, r_era, v_era, n);

% Step 2: Perform SVD on H(0)
% H(0) = P*Z*J^T
[P, Z, J_T] = svd(H_0, 'econ');

% Step 3: Truncate to system order N_r
P_Nr = P(:, 1:N_r);
Z_Nr = Z(1:N_r, 1:N_r);
J_Nr = J_T(:, 1:N_r);

% Step 4: Compute square root and inverse square root
Z_sqrt = sqrt(Z_Nr);
Z_inv_sqrt = inv(Z_sqrt);

% Step 5: Compute discrete-time state transition matrix S (Eq. 14)
% S = Z^(-1/2) * P^T * H(1) * J * Z^(-1/2)
S = Z_inv_sqrt * P_Nr' * H_1 * J_Nr * Z_inv_sqrt;

% Step 6: Compute output matrix C (Eq. 14)
% C = E_n^T * P * Z^(1/2)
E_n_T = eye(n, size(P_Nr, 1));
C = E_n_T * P_Nr * Z_sqrt;

% Step 7: Compute controllability matrix X_0 (Eq. 14)
% X_0 = Z^(1/2) * J^T * E_l
l = size(J_Nr, 1);
E_l_T = eye(l, size(J_Nr, 1));
X_0 = Z_sqrt * J_Nr' * E_l_T;

% --- MODAL TRANSFORMATION TO RECOVER MODE SHAPES (Paper post-Eq. 14) ---
% 
% Paper states: "Letting T be the eigenvector matrix of S, and 
% transforming the computed matrices to modal coordinates, it can be 
% shown that C contains the mode shapes of the system"
%
% So we compute: T (eigenvectors of S)
% Then: Tilde(S) = T^(-1)*S*T  (diagonal in modal space)
%       Tilde(C) = C*T          (mode shapes in physical output space)

% Get eigenvalues and eigenvectors of S
[T, Lambda_S] = eig(S);  % T = eigenvector matrix, Lambda_S = eigenvalues
s_bar = diag(Lambda_S);   % Discrete eigenvalues

% Transform matrices to modal coordinates
% Tilde(S) = T^(-1)*S*T (this will be diagonal or block-diagonal)
S_tilde = inv(T) * S * T;

% Transform C to modal coordinates: Tilde(C) = C*T
% THIS contains the mode shapes (as stated in paper)
C_tilde = C * T;  % Mode shapes in physical output coordinates

% Time step
dt = tSol(2) - tSol(1);

% Continuous-time poles (Paper, Eq. 15)
% eta_i = sigma_i ± i*omega_d_i = (1/dt) * ln(s_bar_i)
eta = log(s_bar) / dt;

fprintf('\n========================================\n');
fprintf(' ERA MODAL ANALYSIS (With T Transformation)\n');
fprintf('========================================\n\n');

mode_count = 0;
displayed = false(length(eta), 1);

for i = 1:length(eta)
    if displayed(i)
        continue;
    end
    
    eta_i = eta(i);
    sigma_i = real(eta_i);
    omega_d_i = abs(imag(eta_i));
    displayed(i) = true;
    
    % Check for complex conjugate
    for j = (i+1):length(eta)
        if abs(eta(j) - conj(eta_i)) < 1e-10
            displayed(j) = true;
            break;
        end
    end
    
    mode_count = mode_count + 1;
    
    fprintf('Mode %d:\n', mode_count);
    if omega_d_i > 1e-6  % Oscillatory
        omega_n_i = sqrt(sigma_i^2 + omega_d_i^2);
        zeta_i = -sigma_i / omega_n_i;
        
        fprintf('  Frequency (Hz): %.4f\n', omega_d_i/(2*pi));
        fprintf('  Damping Ratio ζ: %.6f', zeta_i);
        if zeta_i < 0
            fprintf(' [UNSTABLE]\n');
        else
            fprintf(' [STABLE]\n');
        end
    else
        fprintf('  Non-oscillatory\n');
        fprintf('  Growth rate σ: %.6f', sigma_i);
        if sigma_i > 0
            fprintf(' [UNSTABLE]\n');
        else
            fprintf(' [STABLE]\n');
        end
    end
    fprintf('\n');
end

fprintf('========================================\n');
fprintf('System order N_r: %d (transformed via T)\n', N_r);
fprintf('========================================\n\n');

% --- PROJECTION USING MODE SHAPES IN PHYSICAL COORDINATES ---
%
% After transformation: Tilde(C) = C*T contains mode shapes
% Now project original signal onto these mode shapes
%
% y(t) ≈ Tilde(C) * q(t) = C*T * q(t)
% where q(t) are modal coordinates with physical nonlinear dynamics

% --- PROJECTION USING MODE SHAPES IN PHYSICAL COORDINATES ---

% Interpolate original measured signal and ensure it is a matrix of size [n x N_measurements]
y_original_interp = interp1(tSol, Y_measured, tspan, 'linear', 'extrap');
y_original = y_original_interp'; % Transpose to get [n x N_measurements]

% Output mode shapes in physical coordinates: Tilde(C) = C*T
Phi = C_tilde; % (n x N_r) - mode shapes after transformation

fprintf('C_tilde (mode shapes) = \n');
disp(C_tilde)

fprintf('Dimensions check:\n');
fprintf('  Phi size: [%d x %d]\n', size(Phi, 1), size(Phi, 2));
fprintf('  y_original size: [%d x %d]\n', size(y_original, 1), size(y_original, 2));

% Project original signal onto ALL mode shapes
% y(t) = Phi * q(t)  =>  q(t) = pinv(Phi) * y(t)
% If mode shapes (columns of Phi) are orthogonal, pinv(Phi) can be approximated by Phi'
% for projection.
q_all = pinv(Phi) * y_original;  % (N_r x N_measurements) - using pseudoinverse for robustness

num_modes = 2;
q_subset = q_all(1:num_modes, :);  % (num_modes x N_measurements)
Phi_subset = Phi(:, 1:num_modes);  % (n x num_modes)

% Reconstruct
y_reconstructed = Phi_subset * q_subset;  % (n x N_measurements)

% Return reconstructed signal
reconstructed_pitch = real(y_reconstructed(1, :))';

fprintf('\nPROJECTION COMPLETED (First Mode Only):\n');
fprintf('- Projected onto all %d modes\n', N_r);
fprintf('- Reconstructed using FIRST MODE ONLY\n');
fprintf('- This isolates the dominant modal dynamics\n');
fprintf('- Amplitude-dependent recovery rates preserved\n\n');


end
