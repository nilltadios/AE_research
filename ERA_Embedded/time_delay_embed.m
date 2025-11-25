function [Y_embedded, t_embedded] = time_delay_embed(y, t, dim, lag)
%TIME_DELAY_EMBED Creates a time-delay embedded matrix from a scalar signal.
%
% INPUTS:
%   y   - (N x 1) Scalar time series data (e.g., pitch).
%   t   - (N x 1) Time vector corresponding to y.
%   dim - (Integer) Embedding dimension (Target is 4).
%   lag - (Integer) Delay lag in number of samples.
%
% OUTPUTS:
%   Y_embedded - (dim x M) Embedded state matrix.
%                Column k is [y(k); y(k-tau); ...; y(k-(dim-1)tau)]
%                where tau = lag.
%   t_embedded - (M x 1) Time vector corresponding to the columns of Y_embedded.
%                (Corresponds to the time of the *most recent* sample in the vector).

    N = length(y);
    
    % Calculate the number of available embedded vectors
    % The last required index is: k - (dim-1)*lag >= 1
    % So k >= 1 + (dim-1)*lag
    
    start_idx = 1 + (dim-1)*lag;
    
    if start_idx > N
        error('Signal too short for dimension %d and lag %d.', dim, lag);
    end
    
    M = N - start_idx + 1;
    Y_embedded = zeros(dim, M);
    t_embedded = zeros(M, 1);
    
    for i = 1:M
        % The current time index (most recent value in the window)
        current_idx = start_idx + i - 1;
        
        t_embedded(i) = t(current_idx);
        
        % Build the vector: [y(t); y(t-tau); y(t-2tau); ...]
        for d = 1:dim
            delay_idx = current_idx - (d-1)*lag;
            Y_embedded(d, i) = y(delay_idx);
        end
    end
    
end
