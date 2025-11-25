% =========================================================================
% MAIN BODY - OPTIMIZED FOR PERFORMANCE (EMBEDDED ERA)
% =========================================================================
% Uses Time-Delay Embedding + ERA to reconstruct dynamics from scalar pitch.

clear all;
clearvars -global tspan Y;
clc;

% =========================================================================
% USER CONFIGURATION
% =========================================================================
analysis_option = 2;  % ALWAYS 2 for this script (Reconstructed via Embedded ERA)
PERFORMANCE_MODE = true; 
dt_reduction_factor = 1;

% --- Embedding Parameters ---
% "The main thing that you will need to play with is what the delay time should be."
delay_time_steps = 10; % Play with this! (Number of steps, not seconds)
embedding_dim = 4;     % Fixed at 4 as per instructions

% =========================================================================
% AEROELASTIC PARAMETERS
% =========================================================================
Sup = [1 1.5 0 2 0 0]; 
gammaA = Sup(1,1); 
gammaAAA = Sup(1,2); 
gammaAAAAA = Sup(1,3);
gammaX = Sup(1,4); 
gammaXXX = Sup(1,5); 
gammaXXXXX = Sup(1,6);

% =========================================================================
% INITIAL CONDITION GENERATION
% =========================================================================
Us_initial = 7.45;
ICM = 2;

fprintf('Constructing initial conditions (ICM=%d)...
', ICM);

ode_options_fine = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
ode_options_coarse = odeset('RelTol', 1e-3, 'AbsTol', 1e-5); 
ode_opts = ode_options_fine;
if PERFORMANCE_MODE
    ode_opts = ode_options_coarse;
end

if ICM == 1
    x0_IC = [60*pi/180; 0; 0; 0; 0; 0; 0; 0];
elseif ICM == 2
    gust_type = 2; W0 = 0.1; T = 10; dtg = 0.05; tmaxg = 2*T;
    tspang = 0:dtg:tmaxg; x0g = zeros(8,1);
    [~, xic] = ode45(@(t,x) aeroelastic_Gust(t,x,W0,T,gust_type,...
        gammaA,gammaAAA,gammaAAAAA,gammaX,gammaXXX,gammaXXXXX,Us_initial),...
        tspang, x0g, ode_opts);
    x0_IC = xic(end,:)';
elseif ICM == 3
    t_IM = 0:0.1:3000; UsIM = 8; x0IM = [5*pi/180, 0, 0, 0, 0, 0, 0, 0];
    [~, xIM] = ode45(@(t,x) aeroelastic(t,x,...
        gammaA,gammaAAA,gammaAAAAA,gammaX,gammaXXX,gammaXXXXX,UsIM),...
        t_IM, x0IM, ode_opts);
    x0_IC = xIM(end,:)';
end

fprintf('Initial conditions ready.
');

% =========================================================================
% GROWTH RATE (mu) ANALYSIS
% =========================================================================
usValues = linspace(7.5, 7.6, 3);
dt = 0.1 * dt_reduction_factor;
tmax = 10000;
tspan = 0:dt:tmax;

allGrowthData = cell(1, numel(usValues));

fprintf('\nAnalyzing growth rate (mu) using EMBEDDED ERA...
');
fprintf('  Embedding Dimension: %d
', embedding_dim);
fprintf('  Delay Steps: %d (dt=%.3f s)
', delay_time_steps, dt);

iteration_logs = cell(numel(usValues), 1);

for i_speed = 1:numel(usValues)
    thisUs = usValues(i_speed);
    log_entry = sprintf('Us = %.3f m/s: ', thisUs); 
    
    % --- Solve the ODE ---
    [tSol, xSol] = ode45(@(t,x) aeroelastic(t, x,...
        gammaA,gammaAAA,gammaAAAAA,gammaX,gammaXXX,gammaXXXXX,thisUs),...
        tspan, x0_IC, ode_opts);
    
    pitch_raw = xSol(:, 1);
    
    % --- EMBEDDING & RECONSTRUCTION ---
    log_entry = [log_entry, 'Embedding & ERA... '];
    
    % 1. Create Embedded Coordinates
    [Y_embedded, t_embedded] = time_delay_embed(pitch_raw, tSol, embedding_dim, delay_time_steps);
    
    % 2. Run ERA on Embedded System
    era_params.start_index = 50000; % Ensure this isn't out of bounds for t_embedded
    
    % Adjust start index for potentially shorter embedded signal
    if era_params.start_index > length(t_embedded)
        % If too short, use half the signal
        era_params.start_index = floor(length(t_embedded)/2);
        warning('Adjusted ERA start index to %d due to signal length.', era_params.start_index);
    end
    
    era_params.r_era = 5000;
    era_params.s_era = 5000;
    era_params.n_r = 4; % We want to project onto 2 modes, but system order is likely 4
    
    [reconstructed_pitch, t_recon] = reconstructSignalWithERA_Embedded(Y_embedded, t_embedded, era_params);
    
    % --- Envelope Analysis ---
    % Use reconstructed signal for envelope
    [tData, rData] = plotResponseAndEnvelope(t_recon, [], thisUs, 2, pitch_raw(1:length(t_recon)), true);
    
    % Hack: plotResponseAndEnvelope expects xSol matrix. 
    % Let's create a dummy matrix for it.
    xSol_recon = zeros(length(reconstructed_pitch), 8);
    xSol_recon(:, 1) = reconstructed_pitch;
    
    % Call again with correct format
    [tData, rData] = plotResponseAndEnvelope(t_recon, xSol_recon, thisUs, 2, pitch_raw(1:length(t_recon)), true);

    if isempty(tData) || isempty(rData)
        log_entry = [log_entry, 'FAILED (no envelope data)']
        allGrowthData{i_speed} = struct('r', [], 'mu', [], 'polyfit', []);
        iteration_logs{i_speed} = log_entry;
        continue;
    end
    
    rdot = diff(rData) ./ diff(tData);
    r_for_analysis = (rData(1:end-1) + rData(2:end)) / 2;
    mu_for_analysis = rdot ./ r_for_analysis;
    
    polyCoeffs_mu_r = plotMuVsR(tData, rData, thisUs, true);
    
    allGrowthData{i_speed} = struct('r', r_for_analysis,...
        'mu', mu_for_analysis, 'polyfit', polyCoeffs_mu_r);
    
    log_entry = [log_entry, 'OK'];
    iteration_logs{i_speed} = log_entry;
end

for i = 1:numel(iteration_logs)
    fprintf('%s\n', iteration_logs{i});
end

% ... (Rest of Critical Speed and Bifurcation Logic - Copied from optimized) ...

rList = 0.00:0.005:0.07;
UcList = zeros(size(rList));

fprintf('\nCalculating critical speeds (vectorized)...
');
mu_matrix = zeros(numel(usValues), numel(rList));

for i_speed = 1:numel(usValues)
    if ~isempty(allGrowthData{i_speed}.polyfit) &&... 
            ~any(isnan(allGrowthData{i_speed}.polyfit))
        mu_matrix(i_speed, :) = polyval(allGrowthData{i_speed}.polyfit, rList);
    else
        mu_matrix(i_speed, :) = NaN;
    end
end

for ir = 1:numel(rList)
    mu_vs_Us = mu_matrix(:, ir);
    valid_idx = ~isnan(mu_vs_Us);
    if sum(valid_idx) < 3, UcList(ir) = NaN; continue; end
    
    polyCoeffs = polyfit(usValues(valid_idx), mu_vs_Us(valid_idx), 2);
    a2 = polyCoeffs(1); a1 = polyCoeffs(2); a0 = polyCoeffs(3);
    discriminant = a1^2 - 4*a2*a0;
    
    if discriminant < 0
        UcList(ir) = NaN;
    else
        possibleRoots = [(-a1 + sqrt(discriminant))/(2*a2), (-a1 - sqrt(discriminant))/(2*a2)];
        validMask = (possibleRoots > 0) & (abs(imag(possibleRoots)) < 1e-10);
        positiveRoots = possibleRoots(validMask);
        if ~isempty(positiveRoots)
            [~, minIndex] = min(abs(positiveRoots - 7.6));
            UcList(ir) = positiveRoots(minIndex);
        else
            UcList(ir) = NaN;
        end
    end
end

fprintf('Critical speeds calculated.
');

% Steady State Cache (Using original logic, no ERA needed for this part strictly speaking)
cache_file = 'steady_state_amplitudes.mat';
if isfile(cache_file)
    fprintf('\nLoading cached amplitude data...
');
    load(cache_file);
else
    fprintf('\nCalculating steady-state amplitudes...
');
    usVec = linspace(7.5, 7.7, 100);
    amps_pitch = zeros(size(usVec));
    amps_plunge = zeros(size(usVec));
    for i_amp = 1:length(usVec)
        Us_current = usVec(i_amp);
        [amps_pitch(i_amp), amps_plunge(i_amp)] = amplitudeconver(Us_current, x0_IC,...
            gammaA, gammaAAA, gammaAAAAA, gammaX, gammaXXX, gammaXXXXX);
    end
    save(cache_file, 'usVec', 'amps_pitch', 'amps_plunge');
end

% Plotting
figure('Name', 'Combined Growth Rate vs. Amplitude (Embedded)', 'Color', 'w');
hold on; grid on;
colors = lines(numel(usValues));
for i = 1:numel(usValues)
    if ~isempty(allGrowthData{i}.r) && ~any(isnan(allGrowthData{i}.polyfit))
        plot(allGrowthData{i}.r, allGrowthData{i}.mu, 'o',...
            'MarkerSize', 6, 'Color', colors(i,:), 'DisplayName', sprintf('Us=%.3f', usValues(i)));
        r_fit = linspace(0, 0.1, 200);
        mu_fit = polyval(allGrowthData{i}.polyfit, r_fit);
        plot(r_fit, mu_fit, '-', 'LineWidth', 2, 'Color', colors(i,:), 'HandleVisibility', 'off');
    end
end
title('Growth Rate vs. Amplitude (Embedded ERA)');
xlabel('Amplitude r [rad]'); ylabel('Growth Rate \mu [1/s]');
xlim([0 0.1]); legend('show'); hold off;

figure('Color', 'w', 'Name', 'Bifurcation Diagram: Pitch');
hold on; grid on;
plot(UcList, rList, 'bo-', 'LineWidth', 2, 'DisplayName', 'Unstable Limit Cycle (\mu=0)');
plot(usVec, amps_pitch, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Stable Limit Cycle');
ylim([0 0.1]);
xlabel('Flow Speed Us [m/s]'); ylabel('Pitch Amplitude [rad]');
legend('show'); hold off;

fprintf('\nAnalysis complete!\n');
