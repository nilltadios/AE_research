function polyCoeffs = plotMuVsR(tData, rData, Us, do_plot)
    % Calculates mu (rdot/r) vs. r, fits a curve to a specified domain,
    % plots it, and returns the fit coefficients.
    
    if nargin < 4
        do_plot = true;
    end

    % 1. Calculate rdot = dr/dt using finite differences
    rdot = diff(rData) ./ diff(tData);
    
    % 2. Create corresponding r values at the midpoints for plotting
    r_for_plot = (rData(1:end-1) + rData(2:end)) / 2;
    
    % 3. Calculate the growth rate mu = rdot / r
    mu = rdot ./ r_for_plot;
    
    % 4. Filter data to the specified domain [0, 0.07] for fitting
    fit_domain = [0, 0.07];
    valid_indices = r_for_plot >= fit_domain(1) & r_for_plot <= fit_domain(2);
    r_for_fit = r_for_plot(valid_indices);
    mu_for_fit = mu(valid_indices);
    
    if do_plot
        % Create plot figure
        figure('Name', sprintf('mu vs. r for Us=%.2f', Us));
        hold on;
        
        % Plot ALL calculated data points
        plot(r_for_plot, mu, 'o', 'MarkerFaceColor', 'b', 'DisplayName', 'All Calculated Data');
    end
    
    % 5. Check if there are enough points in the domain to create a reliable fit
    if numel(r_for_fit) < 3
        warning('Fewer than 3 data points in r=[%.2f, %.2f] for Us=%.2f. Cannot create fit.', fit_domain(1), fit_domain(2), Us);
        polyCoeffs = [NaN, NaN, NaN]; % Return NaN coefficients
    else
        % 6. Fit a 2nd-degree polynomial to the FILTERED data
        polyCoeffs = polyfit(r_for_fit, mu_for_fit, 2);
        
        if do_plot
            % 7. Generate points for the smooth fitted curve over the specified domain
            r_fit_line = linspace(fit_domain(1), fit_domain(2), 200);
            mu_fit_line = polyval(polyCoeffs, r_fit_line);
            
            % 8. Plot the fitted curve
            plot(r_fit_line, mu_fit_line, 'r-', 'LineWidth', 2, 'DisplayName', '2nd-Order Fit (r=[0, 0.07])');
        end
    end
    
    if do_plot
        % Add labels, title, and legend
        title(sprintf('Growth Rate vs. Amplitude ($U_s$ = %.2f m/s)', Us), 'Interpreter', 'latex');
        xlabel('Amplitude, $r$ [rad]', 'Interpreter', 'latex');
        ylabel('Growth Rate, $\mu = \dot{r}/r$ [1/s]', 'Interpreter', 'latex');
        xlim([0 0.1]);
        grid on;
        legend('show');
        set(gca, 'FontSize', 12);
        hold off;
    end
end