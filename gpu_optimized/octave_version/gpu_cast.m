function [val, type] = gpu_cast(data)
    % GPU_CAST Converts data to GPU array if available (NVIDIA or OpenCL).
    % Returns the data on GPU (or CPU if not available) and the type ('nvidia', 'ocl', 'cpu').
    
    persistent device_type
    
    if isempty(device_type)
        if exist('gpuDeviceCount', 'file') || exist('gpuDeviceCount', 'builtin')
            try
                if gpuDeviceCount > 0
                    device_type = 'nvidia';
                else
                    device_type = 'cpu';
                end
            catch
                device_type = 'cpu';
            end
        elseif exist('oclArray', 'file') || exist('oclArray', 'builtin') || exist('ocl_device_count', 'file')
             try
                % Try to load the package if not already loaded (Octave specific)
                pkg load ocl; 
             catch
             end
             
             try
                 % Robust check: Try to create a small array on device
                 dummy = oclArray([1]);
                 device_type = 'ocl';
             catch
                device_type = 'cpu';
             end
        else
            device_type = 'cpu';
        end
        fprintf('  [GPU Setup] Detected compute backend: %s\n', device_type);
    end
    
    type = device_type;
    
    if strcmp(device_type, 'nvidia')
        val = gpuArray(data);
    elseif strcmp(device_type, 'ocl')
        try
            % If already OCL (class starts with 'ocl'), just ensure single precision
            if strncmp(class(data), 'ocl', 3)
                 val = data;
            else
                % Ensure data is full (not diagonal/sparse) as oclArray might not support special types
                if ismatrix(data) && ~isstruct(data)
                    data_full = full(data);
                else
                    data_full = data;
                end

                % AMD/Intel consumer GPUs often lack FP64 (double) support in Mesa/Rusticl
                % explicitly cast to single to avoid creation errors.
                if isa(data_full, 'double')
                     val = oclArray(single(data_full));
                else
                     val = oclArray(data_full);
                end
            end
        catch ME
            fprintf('\n[GPU_CAST ERROR] Failed to create oclArray.\n');
            fprintf('  Input Data Class: %s\n', class(data));
            fprintf('  Input Data Typeinfo: %s\n', typeinfo(data));
            fprintf('  Error Message: %s\n', ME.message);
            rethrow(ME);
        end
    else
        val = data;
    end
end