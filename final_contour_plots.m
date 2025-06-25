close all
% Generate contour plots for each variable in the list
psi = linspace(0, 2*pi, 360); % Azimuth angle (Psi)
r = x_blade; % Radial coordinate (r)

variables_to_plot = ["dT_BEMT", "dQ_BEMT", "Vel_BEMT",...
    "AoA_BEMT", "phi_BEMT", "In_axial_BEMT", "In_swirl_BEMT", "dT_PITT", ... 
    "dQ_PITT", "Vel_PITT", "phi_PITT", "In_axial_PITT"];

for var_idx = 1:length(variables_to_plot)
    var_name = variables_to_plot{var_idx};
    prop_plot = eval(var_name);  % Get the variable data dynamically

    % Adjust prop_plot data
    for i = 1:34
        prop_plot{i} = prop_plot{i}(1:end-1);
    end
%
    for i = 1:17
        if size(prop_plot{i}, 1) == 1  % Row array
            prop_plot_single(i,:) = [prop_plot{i} prop_plot{i+17}];
        else  % Column array
            prop_plot_single(i,:) = [prop_plot{i}' prop_plot{i+17}'];
        end
    end

    % Create meshgrid for the plot
    [Psi, R] = meshgrid(psi, r);

    % Convert polar coordinates to Cartesian
    [X, Y] = pol2cart(Psi, R);

% Rotate the data 90 degrees clockwise to match the text labels
X_rotated = Y;  % Swap X and Y to rotate
Y_rotated = -X; % Negate X to rotate 90 degrees clockwise

    % Plot the contour
    figure(var_idx);
    contourf(X_rotated, Y_rotated, prop_plot_single, 17, 'LineColor',  [0.5 0.5 0.5],'LineWidth', 0.5);  % 300 contour levels
    colorbar;
    title(var_name);
    axis equal;
    axis off
    % Save the figure
%     saveas(gcf, [var_name, '_contour_plot.png']);
end

