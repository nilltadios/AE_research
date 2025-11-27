clc; clear;
% Add the path to your functions
addpath('/home/nilllinux22/Downloads/new code');

% Clear any existing figures before starting
close all; 
% Set parameters for the ERA function
era_params.start_index = 1;
era_params.r_era = 2000;
era_params.s_era = 2000;
era_params.n_r = 4;


for i = 1:6
    % Load file
    disp(i)
    filename = "/home/nilllinux22/Downloads/new code/ERA test/ERA_test" + i + ".mat";
    [~, base_name, ~] = fileparts(filename);
    load(filename);
    
    if i< 5
        x_origin = x(:,3);
    else
        x_origin = x(:,4);
    end
    % Run the reconstruction
    x_reconstruct = reconstructSignalWithERA(t, x_origin, t, era_params);
    
    % --- First Figure: Comparison Plot with new figure name ---
    figure('Name', "Signal plot test " + i); % Set the window name here
    plot(t, x_origin, 'k-');    
    hold on;               
    plot(t, x_reconstruct, 'r--');
    hold off;                 
    legend('Original', 'Reconstructed');
    xlabel('Time (t)');
    ylabel('Signal (X)');
    title(base_name + ": Original vs. Reconstructed", 'Interpreter', 'none'); 
    grid on; 
    
    % --- Second Figure: Error Plot with new figure name ---
    figure('Name', "Error plot test " + i); % Set the window name here
    plot(t, x_reconstruct - x_origin, 'r-'); 
    xlabel('Time (t)');
    ylabel('Error (X_{reconstruct} - X_{original})');
    title(base_name + ": Reconstruction Error", 'Interpreter', 'none');
end

