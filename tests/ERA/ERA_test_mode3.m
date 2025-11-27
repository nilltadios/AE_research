clc; clear;

% Add the path to your functions
addpath('/home/nilllinux22/Downloads/new code');
close all;

% ERA parameters
era_params.r_era = 2000;  % Number of rows in Hankel matrix
era_params.s_era = 2000;  % Number of columns in Hankel matrix
era_params.n_r = 4;       % System order (4 = 2 conjugate pairs = 2 physical modes)

% Define time ranges for modal analysis
time_ranges = struct(...
    'name', {'Early (0-33%)', 'Middle (33-67%)', 'Late (67-100%)', 'Full (0-100%)'}, ...
    'start_pct', {0.0, 0.33, 0.67, 0.0}, ...
    'end_pct', {0.33, 0.67, 1.0, 1.0});

% Initialize results storage
results = struct();

fprintf('\n========================================\n');
fprintf('ERA MODAL ANALYSIS - CONJUGATE PAIRS\n');
fprintf('========================================\n');

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

    % Create figure for amplitude vs time plots
    fig = figure('Position', [100, 100, 1600, 900], ...
        'Name', sprintf('Dataset %d: Modal Amplitude vs Time', dataset_idx));

    % Process each time range
    for range_idx = 1:length(time_ranges)
        fprintf('\n  %s\n', time_ranges(range_idx).name);
        fprintf('  %s\n', repmat('-', 1, 50));

        % Extract time segment
        n_total = length(signal);
        idx_start = max(1, round(time_ranges(range_idx).start_pct * n_total) + 1);
        idx_end = round(time_ranges(range_idx).end_pct * n_total);

        signal_segment = signal(idx_start:idx_end);
        t_segment = t(idx_start:idx_end);
        n_segment = length(signal_segment);

        % Build Hankel matrices
        r = min(era_params.r_era, floor(n_segment/2));
        s = min(era_params.s_era, n_segment - r);

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

        % Singular Value Decomposition
        [U, S, V] = svd(H0, 'econ');

        % Select number of modes
        n_modes = min(era_params.n_r, rank(H0));

        % Reduced matrices
        Ur = U(:, 1:n_modes);
        Sr = S(1:n_modes, 1:n_modes);
        Vr = V(:, 1:n_modes);

        % System realization
        sqrt_Sr = sqrt(Sr);
        inv_sqrt_Sr = diag(1./diag(sqrt_Sr));
        A = inv_sqrt_Sr * (Ur' * H1 * Vr) * inv_sqrt_Sr;
        C = Ur(1, :) * sqrt_Sr;  % First row of observability matrix

        % Extract eigenvalues and eigenvectors
        [Psi_r, Lambda_r_matrix] = eig(A);
        lambda_d = diag(Lambda_r_matrix);

        % Convert to continuous time
        lambda_c = log(lambda_d) / dt;

        % Extract modal parameters (full list)
        frequencies_all = abs(lambda_c) / (2*pi);
        damping_all = -real(lambda_c) ./ abs(lambda_c);

        % Identify conjugate pairs and extract unique physical modes
        % Complex eigenvalues come in conjugate pairs
        n_physical_modes = n_modes / 2;
        frequencies = zeros(n_physical_modes, 1);
        damping_ratios = zeros(n_physical_modes, 1);
        residues = zeros(n_physical_modes, 1);

        % Sort by frequency
        [~, sort_idx] = sort(frequencies_all);
        lambda_sorted = lambda_c(sort_idx);
        Psi_sorted = Psi_r(:, sort_idx);

        % Take first mode from each conjugate pair
        for mode = 1:n_physical_modes
            pair_idx = 2*mode - 1;  % Index of first in conjugate pair
            lambda_mode = lambda_sorted(pair_idx);
            psi_mode = Psi_sorted(:, pair_idx);

            frequencies(mode) = abs(lambda_mode) / (2*pi);
            damping_ratios(mode) = -real(lambda_mode) / abs(lambda_mode);

            % Compute residue (modal participation factor)
            residues(mode) = abs(C * psi_mode);
        end

        % Store results
        results(dataset_idx, range_idx).dataset = dataset_idx;
        results(dataset_idx, range_idx).time_range = time_ranges(range_idx).name;
        results(dataset_idx, range_idx).t_start = t_segment(1);
        results(dataset_idx, range_idx).t_end = t_segment(end);
        results(dataset_idx, range_idx).frequencies = frequencies;
        results(dataset_idx, range_idx).damping_ratios = damping_ratios;
        results(dataset_idx, range_idx).residues = residues;
        results(dataset_idx, range_idx).eigenvalues = lambda_c;

        % Print modal parameters
        fprintf('\n  Physical Modes (Conjugate pairs resolved):\n');
        for mode = 1:n_physical_modes
            fprintf('    Mode %d:\n', mode);
            fprintf('      Natural Frequency: %.6f Hz\n', frequencies(mode));
            fprintf('      Damping Ratio:     %.8f\n', damping_ratios(mode));
            fprintf('      Residue (|R|):     %.6e\n', residues(mode));
        end

        % Reconstruct modal contributions for plotting
        subplot(2, 2, range_idx);
        hold on;

        for mode = 1:n_physical_modes
            pair_idx = 2*mode - 1;
            lambda_mode = lambda_sorted(pair_idx);
            psi_mode = Psi_sorted(:, pair_idx);
            R_mode = C * psi_mode;

            % Modal response: x_mode(t) = R * exp(lambda*t)
            modal_amplitude = abs(R_mode) * exp(real(lambda_mode) * (t_segment - t_segment(1)));

            plot(t_segment, modal_amplitude, 'LineWidth', 2, ...
                'DisplayName', sprintf('Mode %d (%.3f Hz)', mode, frequencies(mode)));
        end

        % Also plot original signal for reference
        plot(t_segment, abs(signal_segment), 'k--', 'LineWidth', 1.5, ...
            'DisplayName', 'Original Signal');

        hold off;
        xlabel('Time (s)');
        ylabel('Amplitude');
        title(time_ranges(range_idx).name);
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        set(gca, 'YScale', 'log');  % Log scale to see decay
    end

    sgtitle(sprintf('Dataset %d: Modal Amplitude vs Time', dataset_idx), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig, sprintf('modal_amplitude_time_dataset_%d.png', dataset_idx));
end

% Save results
save('modal_analysis_results.mat', 'results', 'time_ranges');

fprintf('\n========================================\n');
fprintf('SUMMARY TABLE\n');
fprintf('========================================\n\n');

for dataset_idx = 1:6
    fprintf('Dataset %d:\n', dataset_idx);
    fprintf('%-20s | %-6s | %-12s | %-12s\n', 'Time Range', 'Mode', 'Freq (Hz)', 'Damping');
    fprintf('%s\n', repmat('-', 1, 65));

    for range_idx = 1:length(time_ranges)
        freqs = results(dataset_idx, range_idx).frequencies;
        damps = results(dataset_idx, range_idx).damping_ratios;

        for mode = 1:length(freqs)
            if mode == 1
                range_name = time_ranges(range_idx).name;
            else
                range_name = '';
            end
            fprintf('%-20s | %-6d | %12.6f | %12.8f\n', ...
                range_name, mode, freqs(mode), damps(mode));
        end
    end
    fprintf('\n');
end

fprintf('\nAnalysis complete. Results saved to modal_analysis_results.mat\n');
