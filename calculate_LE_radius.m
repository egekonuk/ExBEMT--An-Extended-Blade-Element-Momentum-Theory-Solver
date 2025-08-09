function le_radius_ratio = calculate_LE_radius(dat_filepath)
% CALCULATE_LE_RADIUS - Approximates the leading edge radius from an airfoil .dat file.
%   It reads the coordinate file, finds the points closest to the leading
%   edge (x=0), and fits a circle to these points to find the radius.
%
%   Inputs:
%   dat_filepath - Full path to the airfoil coordinate .dat file.
%
%   Outputs:
%   le_radius_ratio - The calculated leading edge radius non-dimensionalized by the chord.

    if ~exist(dat_filepath, 'file')
        warning('Airfoil .dat file not found: %s', dat_filepath);
        le_radius_ratio = 0.015; % Return a default value
        return;
    end

    % --- Load and Parse Airfoil Data ---
    try
        coords = readmatrix(dat_filepath, 'FileType', 'text', 'CommentStyle', '#');
        % Some .dat files have a name in the first line, handle this
        if isempty(coords)
            opts = detectImportOptions(dat_filepath);
            opts.DataLines = [2, Inf]; % Skip the first line
            coords = readmatrix(dat_filepath, opts);
        end
    catch
        warning('Could not read airfoil file: %s', dat_filepath);
        le_radius_ratio = 0.015; % Return a default value
        return;
    end
    
    % Ensure data is in two columns
    if size(coords, 2) ~= 2
        warning('Invalid format in airfoil file: %s', dat_filepath);
        le_radius_ratio = 0.015; % Return a default value
        return;
    end

    % --- Isolate Leading Edge Points ---
    % Find the point with the minimum x-coordinate (the leading edge)
    [~, min_x_idx] = min(coords(:,1));
    
    % Select a few points around the leading edge for a stable fit.
    % Typically 3-5 points on either side is sufficient.
    num_points_fit = 5; 
    start_idx = max(1, min_x_idx - num_points_fit);
    end_idx = min(size(coords, 1), min_x_idx + num_points_fit);
    
    le_points = coords(start_idx:end_idx, :);
    
    % --- Fit a Circle to the Points (Least Squares Method) ---
    x = le_points(:,1);
    y = le_points(:,2);
    
    % Set up the linear system Ax = b to solve for circle parameters
    A = [2*x, 2*y, ones(size(x))];
    b = x.^2 + y.^2;
    
    % Solve for circle parameters [xc, yc, R^2 - xc^2 - yc^2]
    circle_params = A\b;
    
    xc = circle_params(1);
    yc = circle_params(2);
    
    % Calculate the radius from the solved parameters
    le_radius_ratio = sqrt(circle_params(3) + xc^2 + yc^2);
end
