function [tData, rData] = plotResponseAndEnvelope(tSol, xSol, Us, analysis_option, xSol_ori, do_plot)
%plotResponseAndEnvelope Plots the pitch response and its envelope from provided data.
%
% This function takes the results of an ODE solution, extracts the pitch
% time series, finds its upper and lower envelopes using findpeaks, and
% optionally plots the results.
%
% INPUTS:
%   tSol            - (N x 1 double) Time vector from the ODE solution.
%   xSol            - (N x M double) State solution matrix from the ODE. The
%                     first column is assumed to be the pitch data.
%   Us              - (double) The flow speed corresponding to this solution,
%                     used for the plot title.
%   analysis_option - (integer) Flag to control the plot title.
%                     1: Original Data, 2: Reconstructed Data.
%   xSol_ori        - (N x 1 double) Original pitch data (optional, used if analysis_option=2).
%   do_plot         - (boolean) Optional. If true, generates a plot. Default is true.
%
% OUTPUTS:
%   tData           - (P x 1 double) Time points of the envelope.
%   rData           - (P x 1 double) Amplitude values of the envelope.

    if nargin < 6
        do_plot = true;
    end

    % Extract the pitch data from the first column of the solution matrix
    pitchData = xSol(:,1);
    
    % Find the peaks (upper envelope) and troughs (lower envelope)
    [peakVals, peakLocs] = findpeaks(pitchData, tSol);
    [minVals, minLocs] = findpeaks(-pitchData, tSol);
    
    % Return empty arrays if not enough peaks are found for a reliable envelope
    if isempty(peakVals) || isempty(minVals) || length(peakVals) < 7 || length(minVals) < 7
        tData = []; 
        rData = [];
        warning('Not enough peaks found for Us=%.2f. Cannot create envelope.', Us);
        return;
    end
    
    % Combine peaks and troughs, then sort by time to form the envelope data
    allVals = [peakVals; -minVals];
    allLocs = [peakLocs; minLocs];
    [allLocs, sortIdx] = sort(allLocs);
    allVals = allVals(sortIdx);
    
    % Discard initial transient points to focus on the growing/stable part
    tData = allLocs(7:end);
    rData = abs(allVals(7:end));
    
    % --- Plotting ---
    if do_plot
        figure('Name', sprintf('Response and Envelope for Us=%.2f', Us));
        hold on;
        
        % Set the title based on which data source was used
        if analysis_option == 1
            plotTitle = sprintf('Original Pitch Response and Envelope ($U_s$ = %.2f m/s)', Us);
            plot(tSol, pitchData, 'r-', 'DisplayName', 'Pitch Time Series');
        elseif analysis_option == 2
            plotTitle = sprintf('Reconstructed Pitch Response and Envelope ($U_s$ = %.2f m/s)', Us);
            plot(tSol, xSol_ori, 'r-', 'DisplayName', 'Original Pitch Time Series');
            plot(tSol, pitchData, 'k--', 'DisplayName', 'Reconstructed Pitch Time Series');
        else
            plotTitle = sprintf('Pitch Response and Envelope ($U_s$ = %.2f m/s)', Us);
        end
        
        plot(tData, rData, 'bo', 'MarkerFaceColor','b', 'DisplayName', 'Detected Envelope');
        hold off;
        
        title(plotTitle, 'Interpreter', 'latex');
        xlabel('Time [s]'); 
        ylabel('Pitch Angle [rad]');
        ylim([-0.1 0.1]);
        legend('show', 'Location', 'best'); 
        grid on;
        ax = gca; 
        ax.FontSize = 12;
    end
end