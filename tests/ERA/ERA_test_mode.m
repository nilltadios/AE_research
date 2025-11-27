clc; clear;

% Add the path to your functions
addpath('/home/nilllinux22/Downloads/new code');

% ERA parameters
era_params.r_era = 2000;  % Number of rows in Hankel matrix
era_params.s_era = 2000;  % Number of columns in Hankel matrix
era_params.n_r = 4;       % System order (number of modes * 2)

% Define time ranges for modal analysis
time_ranges = struct(...
    'name', {'Early (0-33%)', 'Middle (33-67%)', 'Late (67-100%)', 'Full (0-100%)'}, ...
    'start_pct', {0.0, 0.33, 0.67, 0.0}, ...
    'end_pct', {0.33, 0.67, 1.0, 1.0});

% Initialize results storage
results = struct();

for dataset_idx = 1:6
    fprintf('\n========================================\n');
    fprintf('Dataset %d\n', dataset_idx);
    fprintf('========================================\n');

    % Load data
    filename = sprintf('/home/nilllinux22/Downloads/new code/ERA test/ERA_test%d.mat', dataset_idx);
    data = load(filename);

    % Handle different field names (x or y)
    if isfield(data, 'x')
        signal_matrix = data.x;
    elseif isfield(data, 'y')
        signal_matrix = data.y;
    else
        error('Dataset %d does not contain field "x" or "y"', dataset_idx);
    end

    % Select signal column
    if dataset_idx < 5
        signal = signal_matrix(:, 3);
    else
        signal = signal_matrix(:, 4);
    end

    t = data.t;
    dt = t(2) - t(1);  % Time step
    fs = 1/dt;         % Sampling frequency

    % Process each time range
    for range_idx = 1:length(time_ranges)
        fprintf('\n  Time Range: %s\n', time_ranges(range_idx).name);

        % Extract time segment
        n_total = length(signal);
        idx_start = max(1, round(time_ranges(range_idx).start_pct * n_total) + 1);
        idx_end = round(time_ranges(range_idx).end_pct * n_total);

        signal_segment = signal(idx_start:idx_end);
        t_segment = t(idx_start:idx_end);
        n_segment = length(signal_segment);

        % Build Hankel matrix using build_block_hankel
        num_outputs = 1;
        r = min(era_params.r_era, floor(n_segment/2));
        s = min(era_params.s_era, n_segment - r);

        % Use build_block_hankel if available, otherwise manual construction
        try
            H0 = build_block_hankel(signal_segment, 1, r, s, num_outputs);
            H1 = build_block_hankel(signal_segment, 2, r, s, num_outputs);
        catch
            % Manual Hankel matrix construction
            H0 = zeros(r, s);
            for row = 1:r
                for col = 1:s
                    H0(row, col) = signal_segment(row + col - 1);
                end
            end

            H1 = zeros(r, s);
            for row = 1:r
                for col = 1:s
                    idx = row + col;
                    if idx <= n_segment
                        H1(row, col) = signal_segment(idx);
                    end
                end
            end
        end

        % Singular Value Decomposition
        [U, S, V] = svd(H0, 'econ');

        % Select number of modes (system order / 2)
        n_modes = min(era_params.n_r, rank(H0));

        % Reduced matrices
        Ur = U(:, 1:n_modes);
        Sr = S(1:n_modes, 1:n_modes);
        Vr = V(:, 1:n_modes);

        % System realization matrices
        sqrt_Sr = sqrt(Sr);
        inv_sqrt_Sr = diag(1./diag(sqrt_Sr));

        % Discrete-time system matrix
        A = inv_sqrt_Sr * (Ur' * H1 * Vr) * inv_sqrt_Sr;

        % Output matrix (Observability)
        C = Ur * sqrt_Sr;

        % Extract modal parameters from eigenvalues
        [Psi_r, Lambda_r_matrix] = eig(A);
        lambda_d = diag(Lambda_r_matrix);  % Discrete eigenvalues

        % Convert to continuous time
        lambda_c = log(lambda_d) / dt;

        % Extract natural frequencies and damping
        frequencies = abs(lambda_c) / (2*pi);      % Natural frequencies (Hz)
        damping_ratios = -real(lambda_c) ./ abs(lambda_c);  % Damping ratios

        % Mode shapes (from eigenvectors projected through output matrix)
        mode_shapes = C * Psi_r;

        % Store results
        results(dataset_idx, range_idx).dataset = dataset_idx;
        results(dataset_idx, range_idx).time_range = time_ranges(range_idx).name;
        results(dataset_idx, range_idx).t_start = t_segment(1);
        results(dataset_idx, range_idx).t_end = t_segment(end);
        results(dataset_idx, range_idx).frequencies = frequencies;
        results(dataset_idx, range_idx).damping_ratios = damping_ratios;
        results(dataset_idx, range_idx).mode_shapes = mode_shapes;
        results(dataset_idx, range_idx).eigenvalues_discrete = lambda_d;
        results(dataset_idx, range_idx).eigenvalues_continuous = lambda_c;
        results(dataset_idx, range_idx).eigenvectors = Psi_r;

        % Display results
        fprintf('    Identified Modes: %d\n', length(frequencies));
        for mode = 1:length(frequencies)
            fprintf('      Mode %d: f = %.4f Hz, ζ = %.6f\n', ...
                mode, frequencies(mode), damping_ratios(mode));
        end
    end
end

% Save results
save('modal_analysis_results.mat', 'results', 'time_ranges');
fprintf('\n========================================\n');
fprintf('Results saved to: modal_analysis_results.mat\n');
fprintf('========================================\n');

% Generate summary visualization
figure('Position', [100, 100, 1400, 900]);

for dataset_idx = 1:6
    % Plot frequencies across time ranges
    subplot(2, 3, dataset_idx);

    n_ranges = length(time_ranges);
    freq_matrix = [];
    x_labels = {};

    for range_idx = 1:n_ranges
        freqs = results(dataset_idx, range_idx).frequencies;
        % Pad with NaN if different number of modes
        max_modes = max(arrayfun(@(r) length(r.frequencies), results(dataset_idx, :)));
        freq_col = [freqs; nan(max_modes - length(freqs), 1)];
        freq_matrix = [freq_matrix, freq_col];
        x_labels{range_idx} = time_ranges(range_idx).name;
    end

    bar(freq_matrix', 'grouped');
    xlabel('Time Range');
    ylabel('Natural Frequency (Hz)');
    title(sprintf('Dataset %d: Modal Frequencies', dataset_idx));
    set(gca, 'XTickLabel', x_labels, 'XTickLabelRotation', 45);
    legend(arrayfun(@(i) sprintf('Mode %d', i), 1:size(freq_matrix,1), 'UniformOutput', false), ...
        'Location', 'best', 'FontSize', 8);
    grid on;
end

sgtitle('Modal Frequency Analysis Across Time Ranges', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'modal_frequencies_summary.png');

% Create damping ratio comparison plot
figure('Position', [100, 100, 1400, 900]);

for dataset_idx = 1:6
    subplot(2, 3, dataset_idx);

    damp_matrix = [];
    x_labels = {};

    for range_idx = 1:n_ranges
        damps = results(dataset_idx, range_idx).damping_ratios;
        max_modes = max(arrayfun(@(r) length(r.damping_ratios), results(dataset_idx, :)));
        damp_col = [damps; nan(max_modes - length(damps), 1)];
        damp_matrix = [damp_matrix, damp_col];
        x_labels{range_idx} = time_ranges(range_idx).name;
    end

    bar(damp_matrix', 'grouped');
    xlabel('Time Range');
    ylabel('Damping Ratio');
    title(sprintf('Dataset %d: Modal Damping', dataset_idx));
    set(gca, 'XTickLabel', x_labels, 'XTickLabelRotation', 45);
    legend(arrayfun(@(i) sprintf('Mode %d', i), 1:size(damp_matrix,1), 'UniformOutput', false), ...
        'Location', 'best', 'FontSize', 8);
    grid on;
end

sgtitle('Modal Damping Ratio Analysis Across Time Ranges', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'modal_damping_summary.png');

fprintf('\nVisualization figures saved.\n');
