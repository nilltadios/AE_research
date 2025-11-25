function [amp_pitch, amp_plunge] = amplitudeconver(Us, x0_IC, ...
                               gammaA, gammaAAA, gammaAAAAA, ...
                               gammaX, gammaXXX, gammaXXXXX)
    % Solves the ODE in chunks until the response amplitude converges.
    % Assumes the 'aeroelastic' ODE function is on the MATLAB path.
    
    odeOptions = odeset('RelTol',1e-8, 'AbsTol',1e-10);
    dtChunk = 200; tEndMax = 1000000; tCurrentEnd = 0;
    ampThreshold = 0.001; prevAmp_pitch = NaN;
    xCurrentIC = x0_IC;

    while tCurrentEnd < tEndMax
        tSpanChunk = [tCurrentEnd, tCurrentEnd + dtChunk];
        [tSol, xSol] = ode45(@(t,x) aeroelastic(t, x, ...
            gammaA, gammaAAA, gammaAAAAA, gammaX, gammaXXX, gammaXXXXX, Us), ...
            tSpanChunk, xCurrentIC, odeOptions);
        
        xCurrentIC = xSol(end,:);
        
        % Analyze the latter half of the chunk to find steady-state amplitude
        skipTime = tSpanChunk(1) + 0.5 * (tSpanChunk(2) - tSpanChunk(1));
        idxSteady = find(tSol >= skipTime);
        
        if numel(idxSteady) < 10
            tCurrentEnd = tSpanChunk(2); continue;
        end
        
        pitchSteady = xSol(idxSteady, 1);
        plungeSteady = xSol(idxSteady, 3);
        
        % Calculate pitch amplitude
        [peakVals, ~] = findpeaks(pitchSteady);
        [minVals, ~] = findpeaks(-pitchSteady);
        if isempty(peakVals) || isempty(minVals)
            tCurrentEnd = tSpanChunk(2); continue;
        end
        newAmp_pitch = 0.5 * (mean(peakVals) - mean(-minVals));
        
        % Calculate plunge amplitude
        [peakVals_p, ~] = findpeaks(plungeSteady);
        [minVals_p, ~] = findpeaks(-plungeSteady);
        if isempty(peakVals_p) || isempty(minVals_p)
            newAmp_plunge = NaN;
        else
            newAmp_plunge = 0.5 * (mean(peakVals_p) - mean(-minVals_p));
        end
        
        % Check for convergence
        if ~isnan(prevAmp_pitch)
            if abs(newAmp_pitch - prevAmp_pitch) / abs(prevAmp_pitch) < ampThreshold
                amp_pitch = newAmp_pitch;
                amp_plunge = newAmp_plunge;
                return; % Converged
            end
        end
        
        prevAmp_pitch = newAmp_pitch;
        tCurrentEnd = tSpanChunk(2);
    end
    
    % Return the last calculated amplitude if max time is reached
    amp_pitch = newAmp_pitch;
    amp_plunge = newAmp_plunge;
end