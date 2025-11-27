clc; clear;

% Load the modal analysis results
load('modal_analysis_results.mat');

fprintf('\n========================================\n');
fprintf('MODE SHAPE ANALYSIS\n');
fprintf('========================================\n');

% Create detailed mode shape visualizations
for dataset_idx = 1:6
    fprintf('\nDataset %d Mode Shapes:\n', dataset_idx);
    fprintf('----------------------------\n');

    figure('Position', [100, 100, 1400, 900], ...
        'Name', sprintf('Mode Shapes - Dataset %d', dataset_idx));

    n_ranges = length(time_ranges);

    for range_idx = 1:n_ranges
        subplot(2, 2, range_idx);

        % Extract mode shapes
        mode_shapes = results(dataset_idx, range_idx).mode_shapes;
        frequencies = results(dataset_idx, range_idx).frequencies;
        damping = results(dataset_idx, range_idx).damping_ratios;

        % Number of spatial points (rows in mode shape matrix)
        n_spatial = size(mode_shapes, 1);
        n_modes = size(mode_shapes, 2);

        % Spatial coordinate (normalized)
        spatial_coord = linspace(0, 1, n_spatial);

        % Plot each mode shape
        hold on;
        for mode = 1:n_modes
            % Normalize mode shape
            mode_shape_normalized = abs(mode_shapes(:, mode)) / max(abs(mode_shapes(:, mode)));

            plot(spatial_coord, mode_shape_normalized, '-o', ...
                'LineWidth', 2, 'DisplayName', ...
                sprintf('Mode %d (%.3f Hz)', mode, frequencies(mode)));
        end
        hold off;

        xlabel('Normalized Spatial Coordinate');
        ylabel('Normalized Mode Shape Amplitude');
        title(sprintf('%s', time_ranges(range_idx).name));
        legend('Location', 'best', 'FontSize', 7);
        grid on;

        % Print mode shape details
        fprintf('  %s:\n', time_ranges(range_idx).name);
        for mode = 1:n_modes
            fprintf('    Mode %d: f=%.4f Hz, ζ=%.6f, |φ|_max=%.4f\n', ...
                mode, frequencies(mode), damping(mode), ...
                max(abs(mode_shapes(:, mode))));
        end
    end

    sgtitle(sprintf('Dataset %d: Mode Shapes Across Time Ranges', dataset_idx), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, sprintf('mode_shapes_dataset_%d.png', dataset_idx));
end

% Create comparison of mode shape evolution across time ranges
for dataset_idx = 1:6
    figure('Position', [100, 100, 1600, 400], ...
        'Name', sprintf('Mode Shape Evolution - Dataset %d', dataset_idx));

    % Focus on first mode only for clarity
    mode_num = 1;

    for range_idx = 1:length(time_ranges)
        subplot(1, 4, range_idx);

        mode_shapes = results(dataset_idx, range_idx).mode_shapes;
        if size(mode_shapes, 2) >= mode_num
            n_spatial = size(mode_shapes, 1);
            spatial_coord = linspace(0, 1, n_spatial);

            % Plot real and imaginary parts
            plot(spatial_coord, real(mode_shapes(:, mode_num)), 'b-', 'LineWidth', 2);
            hold on;
            plot(spatial_coord, imag(mode_shapes(:, mode_num)), 'r--', 'LineWidth', 2);
            hold off;

            xlabel('Normalized Position');
            ylabel('Mode Shape Amplitude');
            title(time_ranges(range_idx).name);
            legend('Real', 'Imaginary', 'Location', 'best', 'FontSize', 8);
            grid on;
        end
    end

    sgtitle(sprintf('Dataset %d: Mode 1 Evolution Across Time Ranges', dataset_idx), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, sprintf('mode1_evolution_dataset_%d.png', dataset_idx));
end

% Create summary table
fprintf('\n========================================\n');
fprintf('MODAL PARAMETER SUMMARY TABLE\n');
fprintf('========================================\n\n');

for dataset_idx = 1:6
    fprintf('Dataset %d:\n', dataset_idx);
    fprintf('%-20s | %-10s | %-10s | %-15s\n', 'Time Range', 'Mode', 'Freq (Hz)', 'Damping Ratio');
    fprintf('%s\n', repmat('-', 1, 70));

    for range_idx = 1:length(time_ranges)
        freqs = results(dataset_idx, range_idx).frequencies;
        damps = results(dataset_idx, range_idx).damping_ratios;

        for mode = 1:length(freqs)
            if mode == 1
                range_name = time_ranges(range_idx).name;
            else
                range_name = '';
            end
            fprintf('%-20s | Mode %-4d | %10.4f | %15.6f\n', ...
                range_name, mode, freqs(mode), damps(mode));
        end
        fprintf('\n');
    end
    fprintf('\n');
end

fprintf('All mode shape visualizations saved.\n');
