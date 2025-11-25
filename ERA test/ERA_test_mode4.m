clc; clear;
close all
% Add the path to your functions
addpath('/home/nilllinux22/Downloads/new code');


% ERA parameters
era_params.r_era = 2000;  % Number of rows in Hankel matrix
era_params.s_era = 2000;  % Number of columns in Hankel matrix
era_params.n_r = 4;       % System order (4 = 2 conjugate pairs = 2 physical modes)

% Define time ranges for modal analysis
time_ranges = struct(...
    'name', {'Early (0-33%)', 'Middle (33-67%)', 'Late (67-100%)', 'Full (0-100%)'}, ...
    'start_pct', {0.0, 0.33, 0.67, 0.0}, ...
    'end_pct', {0.33, 0.67, 1.0, 1.0});

fprintf('\n========================================\n');
fprintf('MODE SHAPE EXTRACTION (NO NORMALIZATION)\n');
fprintf('========================================\n');

for dataset_idx = 1:6
    fprintf('\nProcessing Dataset %d...\n', dataset_idx);

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

    % Create figure for mode shapes
    fig = figure('Position', [100, 100, 1600, 900], ...
        'Name', sprintf('Dataset %d: Mode Shapes Across Time Ranges', dataset_idx));
    
    mode_shape_frommean = [];
    x_coord_mode = [];
    
    % Process each time range
    for range_idx = 1:length(time_ranges)
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

        % Observability matrix (first column gives mode shape along Hankel rows)
        C = Ur * sqrt_Sr;  % r x n_modes matrix

        % Extract eigenvalues and eigenvectors
        [Psi_r, Lambda_r_matrix] = eig(A);
        lambda_d = diag(Lambda_r_matrix);

        % Convert to continuous time
        lambda_c = log(lambda_d) / dt;

        % Extract modal parameters
        frequencies_all = abs(lambda_c) / (2*pi);
        damping_all = -real(lambda_c) ./ abs(lambda_c);

        % Sort by frequency
        [~, sort_idx] = sort(frequencies_all);
        lambda_sorted = lambda_c(sort_idx);
        Psi_sorted = Psi_r(:, sort_idx);
        freq_sorted = frequencies_all(sort_idx);

        % Mode shapes: project eigenvectors through observability matrix
        % This gives the spatial distribution along the Hankel matrix rows
        mode_shapes = C * Psi_sorted;  % r x n_modes

        % Create subplot
        subplot(2, 2, range_idx);
        hold on;

        % Normalized coordinate for x-axis (Hankel row index)
        x_coord = linspace(0, 1, r);

        % Plot all 4 modes (including conjugate pairs) WITHOUT normalization
        colors = {'b', 'r', 'y', 'm'};
        markers = {'o', 'o', 'o', 'o'};

        for mode = 1:n_modes
            % Take absolute value of complex mode shape
            if mode == 1
                mode_shape_abs = abs(mode_shapes(:, mode));
                average_mode_shape_abs = mean(log(mode_shape_abs));
                mode_shape_frommean(range_idx,:) = log(mode_shape_abs) - average_mode_shape_abs; 
                x_coord_mode(range_idx,:) = x_coord;
            end

            plot(x_coord, log(mode_shape_abs), ...
                'Color', colors{mode}, ...
                'Marker', markers{mode}, ...
                'MarkerSize', 4, ...
                'MarkerIndices', 1:50:r, ...
                'LineWidth', 2, ...
                'DisplayName', sprintf('Mode %d (%.3f Hz)', mode, freq_sorted(mode)));
        end

        hold off;
        xlabel('Normalized Spatial Coordinate');
        ylabel('|Mode Shape Amplitude|');
        title(time_ranges(range_idx).name);
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        xlim([0 1]);
   
    end

    figure()
    hold on;
    plot(x_coord_mode(1,:),mode_shape_frommean(1,:),'r-');
    plot(x_coord_mode(2,:),mode_shape_frommean(2,:),'b-');
    plot(x_coord_mode(3,:),mode_shape_frommean(3,:),'k-');
    hold off

    sgtitle(sprintf('Dataset %d: Mode Shapes Across Time Ranges', dataset_idx), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % Save figure
    % saveas(fig, sprintf('mode_shapes_dataset_%d.png', dataset_idx));
end

fprintf('\nMode shape plots saved for all datasets.\n');
