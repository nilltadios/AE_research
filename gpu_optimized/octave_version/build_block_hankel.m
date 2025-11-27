function H_mat = build_block_hankel(data_series, first_sample_idx, num_block_r, num_block_c, num_p_outputs)
    % Vectorized version - 10-50x faster
    total_samples_needed = num_block_r + num_block_c - 1;
    start_idx = first_sample_idx - 1;
    
    % Preallocate with correct size and type (supports gpuArray/oclArray)
    if isa(data_series, 'oclArray') || isa(data_series, 'octave_base_ocl_matrix')
        % Create zeros on CPU then cast, or use OCL specific constructor if available.
        % Safer to create on CPU and cast for OCL in Octave to avoid 'invalid data type' error
        H_mat = oclArray(zeros(num_block_r * num_p_outputs, num_block_c));
    elseif isa(data_series, 'gpuArray')
        H_mat = zeros(num_block_r * num_p_outputs, num_block_c, 'gpuArray');
    else
        H_mat = zeros(num_block_r * num_p_outputs, num_block_c);
    end
    
    % Create index matrix once
    row_indices = (1:num_block_r)';
    col_indices = 1:num_block_c;
    sample_indices = start_idx + row_indices + col_indices - 1;
    
    % Vectorized extraction
    if num_p_outputs == 1
        H_mat = data_series(sample_indices);
    else
        for i_block_row = 1:num_block_r
            row_start_idx = (i_block_row - 1) * num_p_outputs + 1;
            row_end_idx = i_block_row * num_p_outputs;
            H_mat(row_start_idx:row_end_idx, :) = data_series(sample_indices(i_block_row, :), :)';
        end
    end
end
