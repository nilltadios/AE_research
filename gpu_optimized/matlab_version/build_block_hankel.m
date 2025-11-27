function H_mat = build_block_hankel(data_series, first_sample_idx, num_block_r, num_block_c, num_p_outputs)
    % BUILD_BLOCK_HANKEL Vectorized Hankel matrix construction
    % Optimized for MATLAB - 10-50x faster than loop-based version

    start_idx = first_sample_idx - 1;

    % Preallocate with correct size and type
    if isa(data_series, 'gpuArray')
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
