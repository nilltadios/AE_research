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
Y_measured_cpu = xSol(:, 1:4);  % y(t) - measured output (CPU)

% Define system dimensions
n = 4;  % dimensionality of y (outputs)
N_measurements = length(tspan);

% --- ERA Algorithm (Following Paper Equations 13-14) ---

% Prepare measurement data segment ON CPU
Y_segment_cpu = Y_measured_cpu(start_index:end,:);

% Step 1: Build Hankel matrices H(0) and H(1) ON CPU (Faster than OCL for indexing)
H_0_cpu = build_block_hankel(Y_segment_cpu, 1, r_era, v_era, n);
H_1_cpu = build_block_hankel(Y_segment_cpu, 2, r_era, v_era, n);

% --- GPU Acceleration Setup (NVIDIA / AMD-OpenCL / CPU) ---
% Move matrices to GPU now
[H_0, gpu_type] = gpu_cast(H_0_cpu);
[H_1, ~] = gpu_cast(H_1_cpu);

% If using GPU, move other inputs too. 
% Note: We don't need full Y_measured on GPU anymore, just the segment.
if ~strcmp(gpu_type, 'cpu')
    [tSol, ~] = gpu_cast(tSol);
    [tspan, ~] = gpu_cast(tspan);
    % Also move Y_measured for projection later if needed, but projection uses interpolation
    % Interp1 on GPU might be tricky or slow if not supported.
    % Let's keep Y_measured on CPU for interpolation step later unless we optimize that too.
end

% Step 2: Perform SVD on H(0) using Randomized SVD (GPU Accelerated)
% H(0) = P*Z*J^T

try
    % Use rsvd_gpu which offloads heavy multiplications to GPU
    [P_Nr, Z_Nr, J_Nr] = rsvd_gpu(H_0, N_r);
    
    % rsvd_gpu returns results on CPU (U, S, V), so P, Z, J are CPU matrices.
    % This matches our robust CPU-based algebra pipeline below.
    
catch ME
    % Fallback to CPU SVD if GPU fails
    H_0_cpu = safe_gather(H_0);
    [P_cpu, Z_cpu, J_T_cpu] = svd(H_0_cpu, 'econ');

    P_Nr = P_cpu(:, 1:N_r);
    Z_Nr = Z_cpu(1:N_r, 1:N_r);
    J_Nr = J_T_cpu(:, 1:N_r);
end

% Step 3: (Skipped) Truncation is already done by rsvd_gpu or fallback block above

% Step 4: Compute square root and inverse square root
% Perform algebra on CPU for stability (handling complex/real mixes and 'inv')
Z_Nr_cpu = P_Nr; % Just renaming for clarity if needed, but variables are P_Nr, Z_Nr...
% Wait, variables from rsvd_gpu are [U, S, V] corresponding to [P, Z, J_T]
% So P_Nr is U, Z_Nr is S, J_Nr is V (which is J_T' in svd? No svd returns U,S,V where A=U*S*V')
% H0 = P * Z * J^T.  svd(H0) -> [U,S,V]. So P=U, Z=S, J^T=V'. 
% rsvd_gpu returns V (not V'). So J_Nr (which is J) should be V.
% Check rsvd_gpu output: [U, S, V] such that A approx U*S*V'.
% So J^T is V'. J is V.
% J_Nr from rsvd_gpu is V. So it IS J.
% Correct.

% Ensure they are on CPU (rsvd_gpu returns CPU, but just in case of fallback changes)
% Use safe_gather to handle both GPU and CPU arrays
Z_Nr_cpu = safe_gather(Z_Nr);
P_Nr_cpu = safe_gather(P_Nr);
J_Nr_cpu = safe_gather(J_Nr);
H_1_cpu = safe_gather(H_1);

Z_sqrt_cpu = sqrt(Z_Nr_cpu);
Z_inv_sqrt_cpu = inv(Z_sqrt_cpu);

% Step 5: Compute discrete-time state transition matrix S (Eq. 14)
% S = Z^(-1/2) * P^T * H(1) * J * Z^(-1/2)
S = Z_inv_sqrt_cpu * P_Nr_cpu' * H_1_cpu * J_Nr_cpu * Z_inv_sqrt_cpu;

% Step 6: Compute output matrix C (Eq. 14)
% C = E_n^T * P * Z^(1/2)
E_n_T_cpu = eye(n, size(P_Nr_cpu, 1));
C = E_n_T_cpu * P_Nr_cpu * Z_sqrt_cpu;

% Step 7: Compute controllability matrix X_0 (Eq. 14)
% X_0 = Z^(1/2) * J^T * E_l
l = size(J_Nr_cpu, 1);
E_l_T_cpu = eye(l, size(J_Nr_cpu, 1));
X_0 = Z_sqrt_cpu * J_Nr_cpu' * E_l_T_cpu;

% --- MODAL TRANSFORMATION TO RECOVER MODE SHAPES (Paper post-Eq. 14) ---
% 
% Paper states: "Letting T be the eigenvector matrix of S, and 
% transforming the computed matrices to modal coordinates, it can be 
% shown that C contains the mode shapes of the system"
%
% So we compute: T (eigenvectors of S)
% Then: Tilde(S) = T^(-1)*S*T  (diagonal in modal space)
%       Tilde(C) = C*T          (mode shapes in physical output space)

% Get eigenvalues and eigenvectors of S (CPU)
[T, Lambda_S] = eig(S);  % T = eigenvector matrix, Lambda_S = eigenvalues

s_bar = diag(Lambda_S);   % Discrete eigenvalues

% Transform matrices to modal coordinates
% Tilde(S) = T^(-1)*S*T (this will be diagonal or block-diagonal)
S_tilde = inv(T) * S * T;

% Transform C to modal coordinates: Tilde(C) = C*T
% THIS contains the mode shapes (as stated in paper)
C_tilde = C * T;  % Mode shapes in physical output coordinates

% Time step
dt = tSol(2) - tSol(1); % tSol may be oclArray, gather it first
dt = safe_gather(dt);

% Continuous-time poles (Paper, Eq. 15)
% eta_i = sigma_i ± i*omega_d_i = (1/dt) * ln(s_bar_i)
eta = log(s_bar) / dt;

% --- PROJECTION USING MODE SHAPES IN PHYSICAL COORDINATES ---

% Interpolate original measured signal and ensure it is a matrix of size [n x N_measurements]
% Use CPU data for interpolation to avoid moving huge array to GPU if not already there
tSol_cpu = safe_gather(tSol);
tspan_cpu = safe_gather(tspan);
y_original_interp = interp1(tSol_cpu, Y_measured_cpu, tspan_cpu, 'linear', 'extrap');
y_original = y_original_interp'; % Transpose to get [n x N_measurements]

% Move y_original to GPU for projection? 
% No, keep on CPU since Phi (C_tilde) is on CPU and pinv works better on CPU for complex
% [y_original, ~] = gpu_cast(y_original); 

% Output mode shapes in physical coordinates: Tilde(C) = C*T
Phi = C_tilde; % (n x N_r) - mode shapes after transformation

% Project onto first 2 modes
num_modes = 2;
Phi_subset = Phi(:, 1:num_modes);
q_subset = pinv(Phi_subset) * y_original;

% Reconstruct signal
y_reconstructed = Phi_subset * q_subset;

% Return reconstructed pitch signal on CPU
pitch_data = real(y_reconstructed(1, :))';
reconstructed_pitch = safe_gather(pitch_data);

end

function out = safe_gather(data)
    % SAFE_GATHER Safely gather GPU data to CPU, handling already-CPU data
    if isa(data, 'gpuArray') || strncmp(class(data), 'ocl', 3)
        out = gather(data);
    else
        out = data;
    end
end
