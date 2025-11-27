function [val, type] = gpu_cast(data)
    % GPU_CAST Converts data to GPU array if available.
    % For MATLAB: Tries gpuArray (NVIDIA/AMD via Parallel Computing Toolbox)
    % Falls back to CPU with optimized threading.
    %
    % Returns the data on GPU (or CPU if not available) and the type.

    persistent device_type

    if isempty(device_type)
        device_type = 'cpu';  % Default

        % Try MATLAB's gpuArray (works with NVIDIA, and AMD in R2023a+)
        if exist('gpuDeviceCount', 'file') || exist('gpuDeviceCount', 'builtin')
            try
                count = gpuDeviceCount('available');
                if count > 0
                    d = gpuDevice(1);
                    device_type = 'matlab_gpu';
                    fprintf('  [GPU Setup] Using MATLAB gpuArray: %s\n', d.Name);
                end
            catch
                % gpuDevice not available or failed
            end
        end

        % If no GPU, use optimized CPU with threading
        if strcmp(device_type, 'cpu')
            % Ensure MATLAB uses all CPU cores
            try
                maxNumCompThreads('automatic');
            catch
            end
            fprintf('  [GPU Setup] Using optimized CPU (multi-threaded)\n');
        end
    end

    type = device_type;

    if strcmp(device_type, 'matlab_gpu')
        try
            val = gpuArray(data);
        catch ME
            warning('GPU transfer failed: %s. Using CPU.', ME.message);
            val = data;
            type = 'cpu';
        end
    else
        val = data;
    end
end
