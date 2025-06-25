function plot_BEMT_results(J, alpha)
% plot_BEMT_results - Generates well-spaced contour plots with a dynamic title and save functionality.

disp('Generating plots from completed BEMT analysis...');

try
    % Get all necessary variables from the base workspace
    dT_BEMT = evalin('base', 'dT_BEMT');
    dQ_BEMT = evalin('base', 'dQ_BEMT');
    dM_BEMT = evalin('base', 'dM_BEMT');
    dn_BEMT = evalin('base', 'dn_BEMT');
    dT_PITT = evalin('base', 'dT_PITT');
    dQ_PITT = evalin('base', 'dQ_PITT');
    dM_PITT = evalin('base', 'dM_PITT');
    dn_PITT = evalin('base', 'dn_PITT');
    psi_BEMT = evalin('base', 'psi_BEMT');
    x_blade = evalin('base', 'x_blade');
    B = evalin('base', 'B');
catch ME
    error('Could not find required BEMT data in the workspace. Ensure BEMT analysis ran successfully. Error: %s', ME.message);
end

% --- Data Reconstruction ---
plot_titles = {'dT_BEMT (N)', 'dQ_BEMT (N*m)', 'dM_BEMT(N*m)', 'dn_BEMT(N*m)', ...
               'dT_PITT (N)', 'dQ_PITT (N*m)', 'dM_PITT(N*m)', 'dn_PITT(N*m)'};
psi_resolution = 360;
psi_full_range = linspace(0, 2*pi, psi_resolution);
[Psi_mesh, R_mesh] = meshgrid(psi_full_range, x_blade);
all_data_vars = {dT_BEMT, dQ_BEMT, dM_BEMT, dn_BEMT, dT_PITT, dQ_PITT, dM_PITT, dn_PITT};
plot_data = cell(1, length(all_data_vars));

for k_plot = 1:length(all_data_vars)
    reconstructed_matrix = nan(length(x_blade), psi_resolution);
    current_var_cells = all_data_vars{k_plot};
    if isempty(current_var_cells), continue; end
    for i = 1:length(x_blade)
        for bn = 1:B
            if i <= size(psi_BEMT,1) && bn <= size(psi_BEMT,2) && ...
               i <= size(current_var_cells,1) && bn <= size(current_var_cells,2) && ...
               ~isempty(psi_BEMT{i,bn}) && ~isempty(current_var_cells{i,bn})
                data_segment = current_var_cells{i, bn};
                psi_segment = psi_BEMT{i, bn};
                for j_seg = 1:length(psi_segment)
                    [~, nearest_idx] = min(abs(psi_full_range - psi_segment(j_seg)));
                    reconstructed_matrix(i, nearest_idx) = data_segment(j_seg);
                end
            end
        end
    end
    plot_data{k_plot} = reconstructed_matrix;
end

% --- Create a NEW figure with a Tiled Layout for better spacing ---
dynamic_title = sprintf('Sectional Properties | J = %.2f, AOI = %.1f deg', J, alpha);
fig_handle = figure('Name', dynamic_title, 'NumberTitle', 'off');

% Add a "File" menu
file_menu = uimenu(fig_handle, 'Text', '&File');
uimenu(file_menu, 'Text', '&Save Figure', 'Callback', @save_figure_callback, 'Accelerator', 'S');

% Create a 2x4 tiled layout with compact spacing
t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, dynamic_title, 'FontSize', 14, 'FontWeight', 'bold'); % Add a main title

% --- Plotting Logic ---
[X, Y] = pol2cart(Psi_mesh, R_mesh);
X_rotated = Y;
Y_rotated = -X;

for k = 1:length(plot_data)
    nexttile; % Move to the next tile in the layout
    if ~all(isnan(plot_data{k}(:)))
        contourf(X_rotated, Y_rotated, plot_data{k}, 17, 'LineColor', 'none');
    end
    colorbar;
    title(strrep(plot_titles{k}, '_', '\_'));
    axis equal;
    axis off;
end

disp('Plotting complete. Go to File -> Save Figure to save the image.');

% --- Nested function for the Save Callback ---
    function save_figure_callback(~, ~)
        % Create a folder named 'Saved_Plots' if it doesn't exist
        folder_name = 'Saved_Plots';
        if ~exist(folder_name, 'dir')
           mkdir(folder_name);
           disp(['Created folder: ./' folder_name]);
        end

        % Create a descriptive filename
        file_name = sprintf('BEMT_Plots_J_%.2f_AoA_%.1f.png', J, alpha);
        full_path = fullfile(folder_name, file_name);
        
        % Save the figure
        exportgraphics(fig_handle, full_path, 'Resolution', 300);
        disp(['Figure saved successfully to: ' full_path]);
    end

end