function [U, S, V] = rsvd_gpu(A, k, p)
    % RSVD_GPU Randomized SVD optimized for MATLAB
    %
    % Uses GPU if available (NVIDIA/AMD), otherwise uses optimized CPU.
    % MATLAB's built-in BLAS/LAPACK are highly optimized for multi-core CPUs.
    %
    % Usage: [U, S, V] = rsvd_gpu(A, k, p)
    %
    % Inputs:
    %   A - Matrix (gpuArray or CPU matrix)
    %   k - Target rank (number of singular values to keep)
    %   p - Oversampling parameter (default 10)
    %
    % Outputs:
    %   U, S, V - Singular value decomposition such that A approx U*S*V'

    if nargin < 3
        p = 10;
    end

    [m, n] = size(A);
    target_rank = min(k + p, min(m, n));

    % Check if GPU is available
    [A_compute, type] = gpu_cast(A);
    using_gpu = strcmp(type, 'matlab_gpu');

    % 1. Generate Random Test Matrix Omega
    if using_gpu
        Omega = gpuArray.randn(n, target_rank);
    else
        Omega = randn(n, target_rank);
    end

    % 2. Compute Sketch Y = A * Omega
    Y = A_compute * Omega;

    % 3. QR Decomposition
    if using_gpu
        [Q, ~] = qr(Y, 0);  % MATLAB gpuArray supports QR
    else
        [Q, ~] = qr(Y, 0);
    end

    % 4. Form Reduced Matrix B = Q' * A
    B = Q' * A_compute;

    % 5. SVD of Small Matrix B
    if using_gpu
        [Uhat, S, V] = svd(gather(B), 'econ');
        Q = gather(Q);
    else
        [Uhat, S, V] = svd(B, 'econ');
    end

    % 6. Recover U = Q * Uhat
    U = Q * Uhat;

    % Truncate to k
    U = U(:, 1:k);
    S = S(1:k, 1:k);
    V = V(:, 1:k);
end
