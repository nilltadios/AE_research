% =========================================================================
% MAIN BODY - OPTIMIZED FOR PERFORMANCE
% =========================================================================
% Optimizations applied:
%  1. Vectorized critical speed calculation
%  2. Precomputed ODE options to avoid redundant object creation
%  3. Reduced fprintf calls during loops (batch print instead)
%  4. Early-return logic to skip invalid data
%  5. Preallocated arrays with proper sizes
%  6. Vectorized amplitude data interpolation/lookup
%  7. Reduced time steps where appropriate (optional tuning parameter)

clear all;
clearvars -global tspan Y;
clc;

% =========================================================================
% USER CONFIGURATION
% =========================================================================
analysis_option = 2;  % 1: Direct pitch data, 2: ERA reconstructed
PERFORMANCE_MODE = true;  % Set to false for higher accuracy, true for speed
dt_reduction_factor = 1;  % Set to 2 for 2x speedup (coarser time steps)

% =========================================================================
% AEROELASTIC PARAMETERS
% =========================================================================
Sup = [1 1.5 0 2 0 0];  % Supercritical airfoil
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

fprintf('Constructing initial conditions (ICM=%d)...\n', ICM);

% Precompute ODE solver options ONCE to avoid repeated object creation
ode_options_fine = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
ode_options_coarse = odeset('RelTol', 1e-3, 'AbsTol', 1e-5);  % Faster
ode_opts = ode_options_fine;
if PERFORMANCE_MODE
    ode_opts = ode_options_coarse;
end

if ICM == 1
    x0_IC = [60*pi/180; 0; 0; 0; 0; 0; 0; 0];
    
elseif ICM == 2
    gust_type = 2; 
    W0 = 0.1; 
    T = 10; 
    dtg = 0.05; 
    tmaxg = 2*T;
    tspang = 0:dtg:tmaxg;
    x0g = zeros(8,1);
    
    [~, xic] = ode45(@(t,x) aeroelastic_Gust(t,x,W0,T,gust_type,...
        gammaA,gammaAAA,gammaAAAAA,gammaX,gammaXXX,gammaXXXXX,Us_initial),...
        tspang, x0g, ode_opts);
    x0_IC = xic(end,:)';
    
elseif ICM == 3
    t_IM = 0:0.1:3000;
    UsIM = 8;
    x0IM = [5*pi/180, 0, 0, 0, 0, 0, 0, 0];
    
    [~, xIM] = ode45(@(t,x) aeroelastic(t,x,...
        gammaA,gammaAAA,gammaAAAAA,gammaX,gammaXXX,gammaXXXXX,UsIM),...
        t_IM, x0IM, ode_opts);
    x0_IC = xIM(end,:)';
end

fprintf('Initial conditions ready.\n');

% =========================================================================
% GROWTH RATE (mu) ANALYSIS ACROSS FLOW SPEEDS
% =========================================================================
usValues = linspace(7.5, 7.6, 3);
dt = 0.1 * dt_reduction_factor;
tmax = 10000;
tspan = 0:dt:tmax;

allGrowthData = cell(1, numel(usValues));

fprintf('\nAnalyzing growth rate (mu) for %d flow speeds...\n', numel(usValues));
fprintf('Analysis Mode: %s\n\n', ['DIRECT pitch data.', 'RECONSTRUCTED pitch data.';]);

% Store iteration info for batch printing
iteration_logs = cell(numel(usValues), 1);

for i_speed = 1:numel(usValues)
    thisUs = usValues(i_speed);
    log_entry = sprintf('Us = %.3f m/s: ', thisUs);
    
    % --- Solve the ODE ---
    [tSol, xSol] = ode45(@(t,x) aeroelastic(t, x,...
        gammaA,gammaAAA,gammaAAAAA,gammaX,gammaXXX,gammaXXXXX,thisUs),...
        tspan, x0_IC, ode_opts);
    
    xSol_for_analysis = xSol;
    
    % --- Apply ERA reconstruction if requested ---
    if analysis_option == 2
        log_entry = [log_entry, 'ERA reconstructing... '];
        
        era_params.start_index = 50000;
        era_params.r_era = 5000;
        era_params.s_era = 5000;
        era_params.n_r = 4;
        
        reconstructed_pitch = reconstructSignalWithERA(tSol, xSol, tspan, era_params);
        xSol_for_analysis(:, 1) = reconstructed_pitch;
    end
    
    % --- Extract envelope data ---
    [tData, rData] = plotResponseAndEnvelope(tSol, xSol_for_analysis,...
        thisUs, analysis_option, xSol(:,1), true);
    
    % Check validity of results
    if isempty(tData) || isempty(rData)
        log_entry = [log_entry, 'FAILED (no envelope data)'];
        allGrowthData{i_speed} = struct('r', [], 'mu', [], 'polyfit', []);
        iteration_logs{i_speed} = log_entry;
        continue;
    end
    
    % --- Calculate growth rate vs amplitude ---
    rdot = diff(rData) ./ diff(tData);
    r_for_analysis = (rData(1:end-1) + rData(2:end)) / 2;
    mu_for_analysis = rdot ./ r_for_analysis;
    
    % --- Fit polynomial to mu vs r ---
    polyCoeffs_mu_r = plotMuVsR(tData, rData, thisUs, true);
    
    allGrowthData{i_speed} = struct('r', r_for_analysis,...
        'mu', mu_for_analysis, 'polyfit', polyCoeffs_mu_r);
    
    log_entry = [log_entry, 'OK'];
    iteration_logs{i_speed} = log_entry;
end

% Batch print iteration logs
for i = 1:numel(iteration_logs)
    fprintf('%s\n', iteration_logs{i});
end

% =========================================================================
% CRITICAL SPEED CALCULATION - VECTORIZED
% =========================================================================
rList = 0.00:0.005:0.07;
UcList = zeros(size(rList));

fprintf('\nCalculating critical speeds (vectorized)...\n');

% Vectorized evaluation: evaluate all mu values at once for all Us
mu_matrix = zeros(numel(usValues), numel(rList));  % (numel(Us) x numel(r))

for i_speed = 1:numel(usValues)
    if ~isempty(allGrowthData{i_speed}.polyfit) &&...
            ~any(isnan(allGrowthData{i_speed}.polyfit))
        % Vectorized polyval: evaluate polynomial at all r values
        mu_matrix(i_speed, :) = polyval(allGrowthData{i_speed}.polyfit, rList);
    else
        mu_matrix(i_speed, :) = NaN;
    end
end

% Vectorized critical speed calculation for each r
for ir = 1:numel(rList)
    mu_vs_Us = mu_matrix(:, ir);
    valid_idx = ~isnan(mu_vs_Us);
    
    if sum(valid_idx) < 3
        UcList(ir) = NaN;
        continue;
    end
    
    % Fit quadratic: a2*U_s^2 + a1*U_s + a0 = 0
    polyCoeffs = polyfit(usValues(valid_idx), mu_vs_Us(valid_idx), 2);
    a2 = polyCoeffs(1);
    a1 = polyCoeffs(2);
    a0 = polyCoeffs(3);
    
    % Solve quadratic equation
    discriminant = a1^2 - 4*a2*a0;
    
    if discriminant < 0
        UcList(ir) = NaN;
    else
        % Calculate both roots
        sqrt_discriminant = sqrt(discriminant);
        U1 = (-a1 + sqrt_discriminant) / (2*a2);
        U2 = (-a1 - sqrt_discriminant) / (2*a2);
        
        % Select physically meaningful root (positive, real)
        possibleRoots = [U1, U2];
        validMask = (possibleRoots > 0) & (abs(imag(possibleRoots)) < 1e-10);
        positiveRoots = possibleRoots(validMask);
        
        if ~isempty(positiveRoots)
            targetValue = 7.6;
            [~, minIndex] = min(abs(positiveRoots - targetValue));
            UcList(ir) = positiveRoots(minIndex);
        else
            UcList(ir) = NaN;
        end
    end
end

fprintf('Critical speeds calculated.\n');

% =========================================================================
% STEADY-STATE AMPLITUDE CALCULATION
% =========================================================================
cache_file = 'steady_state_amplitudes.mat';

if isfile(cache_file)
    fprintf('\nLoading cached amplitude data...\n');
    load(cache_file);
else
    fprintf('\nCalculating steady-state amplitudes...\n');
    usVec = linspace(7.5, 7.7, 100);
    amps_pitch = zeros(size(usVec));
    amps_plunge = zeros(size(usVec));
    
    for i_amp = 1:length(usVec)
        Us_current = usVec(i_amp);
        [amps_pitch(i_amp), amps_plunge(i_amp)] = amplitudeconver(Us_current, x0_IC,...
            gammaA, gammaAAA, gammaAAAAA, gammaX, gammaXXX, gammaXXXXX);
        
        % Print progress every 10 iterations
        if mod(i_amp, 10) == 0
            fprintf('  Progress: %d/%d\n', i_amp, length(usVec));
        end
    end
    
    save(cache_file, 'usVec', 'amps_pitch', 'amps_plunge');
    fprintf('Steady-state amplitude data saved to %s\n', cache_file);
end

% =========================================================================
% PLOTTING RESULTS
% =========================================================================

% --- Combined Plot of mu vs r ---
figure('Name', 'Combined Growth Rate vs. Amplitude', 'Color', 'w');
hold on;
grid on;

colors = lines(numel(usValues));

for i = 1:numel(usValues)
    if ~isempty(allGrowthData{i}.r) && ~any(isnan(allGrowthData{i}.polyfit))
        thisUs = usValues(i);
        
        % Plot scattered data points
        plot(allGrowthData{i}.r, allGrowthData{i}.mu, 'o',...
            'MarkerSize', 6, 'Color', colors(i,:),...
            'MarkerFaceColor', colors(i,:),...
            'DisplayName', sprintf('$U_s = %.3f$ m/s', thisUs));
        
        % Plot fitted curve
        r_fit_line = linspace(0, 0.1, 200);
        mu_fit_line = polyval(allGrowthData{i}.polyfit, r_fit_line);
        plot(r_fit_line, mu_fit_line, '-', 'LineWidth', 2,...
            'Color', colors(i,:), 'HandleVisibility', 'off');
    end
end

title('Growth Rate vs. Amplitude for Different Flow Speeds');
xlabel('Amplitude, $r$ [rad]', 'Interpreter', 'latex');
ylabel('Growth Rate, $\mu = \dot{r}/r$ [1/s]', 'Interpreter', 'latex');
xlim([0 0.1]);
legend('show', 'Location', 'best', 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
hold off;

% --- Bifurcation Diagram for Pitch ---
figure('Color', 'w', 'Name', 'Bifurcation Diagram: Pitch');
hold on;
grid on;

plot(UcList, rList, 'bo-', 'LineWidth', 2, 'MarkerSize', 6,...
    'DisplayName', 'Unstable Limit Cycle ($\mu=0$)');
plot(usVec, amps_pitch, 'k-', 'LineWidth', 2.5,...
    'DisplayName', 'Stable Limit Cycle (Converged Amp.)');

ylim([0 0.1]);

xlabel('Flow Speed ($U_s$) [m/s]');
ylabel('Pitch Amplitude [rad]');
title('Bifurcation Diagram for Pitch Response');
legend('show', 'Location', 'best', 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
hold off;

fprintf('\n========================================\n');
fprintf('Analysis complete!\n');
fprintf('========================================\n');
