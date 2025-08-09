function installNeuralFoilDependencies()
% installNeuralFoilDependencies Automates the installation of Python
% dependencies for NeuralFoil, including AeroSandbox, CasADi, SciPy,
% and SortedContainers.

    % Define the Python packages required
    required_packages = {
        'neuralfoil',
        'aerosandbox',
        'numpy',
        'casadi',
        'scipy',
        'sortedcontainers',
        'dill',
        'matplotlib',
        'pandas',
        'requests',
    };

    fprintf('Checking Python environment and dependencies for NeuralFoil...\n');

    try
        % 1. Check if a Python environment is configured in MATLAB
        py_env = pyenv;
        if isempty(py_env.Executable)
            warning('MATLAB''s Python environment is not configured. Attempting to auto-detect or configure default Python.');
            % Attempt to set default Python (MATLAB will try to find one)
            pyenv; % Call pyenv without arguments to auto-detect
            py_env = pyenv; % Get the updated environment
            if isempty(py_env.Executable)
                 error('Could not auto-detect a Python environment. Please configure one manually using pyenv(''Version'', ''path/to/python.exe'').');
            end
        end
        fprintf('MATLAB is using Python executable: %s\n', py_env.Executable);

        % Get the path to the pip executable for the active Python environment
        [python_dir, ~, ~] = fileparts(py_env.Executable);
        if ispc % Windows
            pip_executable = fullfile(python_dir, 'Scripts', 'pip.exe');
        else % macOS/Linux
            % For macOS/Linux, first check for 'pip' or 'pip3' in the same 
            % directory as the python executable. This handles system Python.
            pip_path_1 = fullfile(python_dir, 'pip');
            pip_path_2 = fullfile(python_dir, 'pip3');
            
            if isfile(pip_path_1)
                pip_executable = pip_path_1;
            elseif isfile(pip_path_2)
                pip_executable = pip_path_2;
            else
                % Fallback for virtual environments where pip is in a parallel bin folder
                pip_executable = fullfile(python_dir, 'bin', 'pip');
            end
        end

        % 2. Check and install each required package
        for i = 1:length(required_packages)
            pkg = required_packages{i};
            fprintf('Checking for Python package: %s...\n', pkg);

            % Use pip show to check if the package is installed
            [status, cmdout] = system(sprintf('"%s" show %s', pip_executable, pkg));

            if status == 0 % pip show command succeeded, package might be installed
                if contains(lower(cmdout), lower(sprintf('Name: %s', pkg)))
                    fprintf('%s is already installed.\n', pkg);
                else
                    % This case should be rare if status == 0, but good to be robust
                    fprintf('%s not found or installed incorrectly. Attempting installation...\n', pkg);
                    install_package(pip_executable, pkg);
                end
            else % pip show command failed, likely package not found
                fprintf('%s not found. Attempting installation...\n', pkg);
                install_package(pip_executable, pkg);
            end
        end

        fprintf('\nAll NeuralFoil Python dependencies are installed and checked.\n');
        fprintf('Please restart MATLAB if a new module is installed to ensure all new modules are loaded correctly.\n');

    catch ME
        fprintf(2, 'An error occurred during Python dependency installation: %s\n', ME.message);
        fprintf(2, 'Please ensure Python is installed and configured correctly with MATLAB, and you have internet access.\n');
        rethrow(ME); % Re-throw the error to stop execution if critical
    end

end

function install_package(pip_executable, pkg_name)
    % Helper function to execute pip install
    cmd = sprintf('"%s" install %s', pip_executable, pkg_name);
    fprintf('Executing: %s\n', cmd);
    [install_status, install_cmdout] = system(cmd);
    if install_status == 0
        fprintf('%s installed successfully.\\n', pkg_name);
    else
        fprintf(2, 'Failed to install %s:\\n%s\\n', pkg_name, install_cmdout);
        error('Failed to install %s. Please check your internet connection and permissions.', pkg_name);
    end
end
