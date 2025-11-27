function [U, S, V] = rsvd_gpu(A, k, p)
    % rsvd_gpu Randomized SVD optimized for GPU (OpenCL)
    %
    % Usage: [U, S, V] = rsvd_gpu(A, k, p)
    %
    % Inputs:
    %   A - Matrix (gpuArray/oclArray or CPU matrix)
    %   k - Target rank (number of singular values to keep)
    %   p - Oversampling parameter (default 10)
    %
    % Outputs:
    %   U, S, V - Singular value decomposition such that A approx U*S*V'
    %
    % Algorithm:
    % 1. Omega = randn(n, k+p)
    % 2. Y = A * Omega          (Computed on GPU)
    % 3. [Q, ~] = qr(Y, 0)      (Computed on CPU, Y is gathered)
    % 4. B = Q' * A             (Computed on GPU)
    % 5. [Uhat, S, V] = svd(B)  (Computed on CPU)
    % 6. U = Q * Uhat           (Computed on GPU or CPU)

    if nargin < 3
        p = 10;
    end
    
    [m, n] = size(A);
    target_rank = k + p;
    
    % Ensure A is on GPU
    [A_gpu, type] = gpu_cast(A);
    using_gpu = ~strcmp(type, 'cpu');
    
    % 1. Generate Random Test Matrix Omega
    % We generate on CPU and cast to GPU to ensure compatibility
    Omega_cpu = randn(n, target_rank);
    if using_gpu
        [Omega, ~] = gpu_cast(Omega_cpu);
    else
        Omega = Omega_cpu;
    end
    
    % 2. Compute Sketch Y = A * Omega (HEAVY LIFT - ON GPU)
    Y = A_gpu * Omega;

    % 3. QR Decomposition (on CPU - OpenCL often lacks QR)
    if using_gpu
        Y_cpu = gather(Y);
    else
        Y_cpu = Y;
    end
    [Q_cpu, ~] = qr(Y_cpu, 0);

    % 4. Form Reduced Matrix B = Q' * A (HEAVY LIFT - ON GPU)
    if using_gpu
        [Q_gpu, ~] = gpu_cast(Q_cpu);
        B = Q_gpu' * A_gpu;
    else
        B = Q_cpu' * A;
    end

    % 5. SVD of Small Matrix B (ON CPU)
    if using_gpu
        B_cpu = gather(B);
    else
        B_cpu = B;
    end
    
    [Uhat, S, V] = svd(B_cpu, 'econ');
    
    % 6. Recover U = Q * Uhat
    % We can do this on CPU or GPU. Since results U, S, V are needed for
    % subsequent steps which might be CPU based in the current pipeline,
    % let's return them as CPU matrices to be safe, or GPU if requested.
    % For this specific ERA pipeline, keeping them on CPU is safer given previous errors.
    
    U = Q_cpu * Uhat;
    
    % Truncate to k
    U = U(:, 1:k);
    S = S(1:k, 1:k);
    V = V(:, 1:k);
end
