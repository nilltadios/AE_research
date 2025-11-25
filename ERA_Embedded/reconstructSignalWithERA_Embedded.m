function [reconstructed_pitch, t_recon] = reconstructSignalWithERA_Embedded(Y_embedded, t_embedded, era_params)
%reconstructSignalWithERA_Embedded Performs ERA on Time-Delay Embedded coordinates.
% 
% INPUTS:
%   Y_embedded - (dim x M) Embedded state matrix. dim=4 usually.
%   t_embedded - (M x 1) Time vector.
%   era_params - Struct with fields: r_era, s_era, n_r (system order).
% 
% OUTPUTS:
%   reconstructed_pitch - (M x 1) Reconstructed scalar pitch signal.
%   t_recon            - (M x 1) Time vector (same as t_embedded).

    % --- Parameter Extraction ---
    r_era = era_params.r_era;
    v_era = era_params.s_era;
    N_r = era_params.n_r;
    
    % Y_embedded is [n_outputs x n_samples]
    % For ERA, we need it in the format expected by build_block_hankel
    % usually build_block_hankel expects [n_samples x n_outputs] ??
    % Let's check build_block_hankel.m usage in previous file:
    % Y_segment = Y_measured(start_index:end,:); (N x n)
    % H_0 = build_block_hankel(Y_segment, ...);
    
    % So we need to Transpose Y_embedded to be (M x dim)
    Y_for_ERA = Y_embedded'; 
    
    [M, n] = size(Y_for_ERA); % M samples, n=4 dimensions
    
    % --- ERA Algorithm ---
    
    % Step 1: Build Hankel matrices
    % Note: r_era + v_era must be < M
    if (r_era + v_era + 2) > M
        warning('Data length (%d) too short for requested Hankel size (%d+%d). Truncating parameters.', M, r_era, v_era);
        r_era = floor(M/3);
        v_era = floor(M/3);
    end

    H_0 = build_block_hankel(Y_for_ERA, 1, r_era, v_era, n);
    H_1 = build_block_hankel(Y_for_ERA, 2, r_era, v_era, n);
    
    % Step 2: SVD
    [P, Z, J_T] = svd(H_0, 'econ');
    
    % Step 3: Truncate to system order N_r
    % Ensure N_r isn't larger than the available rank
    N_r = min(N_r, size(Z,1));
    
    P_Nr = P(:, 1:N_r);
    Z_Nr = Z(1:N_r, 1:N_r);
    J_Nr = J_T(:, 1:N_r);
    
    % Step 4: SVD-based formulation
    Z_sqrt = sqrt(Z_Nr);
    Z_inv_sqrt = inv(Z_sqrt);
    
    S = Z_inv_sqrt * P_Nr' * H_1 * J_Nr * Z_inv_sqrt;
    
    E_n_T = eye(n, size(P_Nr, 1));
    C = E_n_T * P_Nr * Z_sqrt;
    
    % --- Modal Transformation ---
    [T, Lambda_S] = eig(S);
    s_bar = diag(Lambda_S);
    
    % Mode shapes in embedded coordinate space (4 x N_r)
    C_tilde = C * T;
    
    % --- Modal Analysis Info ---
    dt = t_embedded(2) - t_embedded(1);
    eta = log(s_bar) / dt;
    
    fprintf('\n--- ERA Embedded Analysis (Dim=%d) ---\n', n);
    for i = 1:length(eta)
        sigma = real(eta(i));
        omega = abs(imag(eta(i)));
        freq_hz = omega / (2*pi);
        fprintf('Mode %d: %.2f Hz, Growth=%.4f\n', i, freq_hz, sigma);
    end
    
    % --- Projection and Reconstruction ---
    
    % The measured data Y_embedded (dim x M) corresponds to:
    % Y(t) = C_tilde * q(t)
    % We want to find q(t).
    
    Phi = C_tilde; % (4 x N_r) 
    
    % Project data onto modes: q(t) = pinv(Phi) * Y(t)
    % Y_embedded is (4 x M)
    q_all = pinv(Phi) * Y_embedded; % (N_r x M)
    
    % Isolate dominant modes (e.g., Mode 1 which is a pair usually)
    % We'll just take the specified N_r modes (which is q_all)
    % If we want to filter, we zero out rows in q_all.
    
    % For now, let's assume N_r=4 (2 modes) and we keep them all, 
    % or we can select the "most energetic" or simply the first mode pair.
    
    % Let's reconstruct using the *First Mode Pair* (assumed dominant)
    % Usually sorted by magnitude? eigenvalues of S are not necessarily sorted.
    % We should sort by energy or amplitude.
    
    % Simple strategy: Reconstruct using ALL N_r identified modes first.
    % If the user specifically wants the "first two modes", we should sort.
    
    % Sorting by modal amplitude contribution
    modal_energy = sum(abs(q_all).^2, 2);
    [~, sort_idx] = sort(modal_energy, 'descend');
    
    % Keep top 2 modes (usually a conjugate pair count as 2 dimensions)
    keep_indices = sort_idx(1:min(2, N_r)); 
    % Actually, if it's a conjugate pair, we need both. 
    % Let's just keep all N_r modes for the "Reconstructed" signal 
    % or strictly follow the "first two modes" instruction.
    
    % Instruction: "project the dynamics down onto the first two modes"
    % This implies N_r might be > 2, but we only use 2 for reconstruction.
    
    % Let's pick the 2 most energetic modes.
    q_subset = zeros(size(q_all));
    q_subset(keep_indices, :) = q_all(keep_indices, :);
    
    % Reconstruct
    Y_recon_embedded = Phi * q_subset; % (4 x M)
    
    % The first row corresponds to y(t) (no delay)
    reconstructed_pitch = real(Y_recon_embedded(1, :))';
    t_recon = t_embedded;
    
end
