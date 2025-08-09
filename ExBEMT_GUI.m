classdef ExBEMT_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        %% --- UI FIGURE AND PANEL COMPONENTS ---
        UIFigure          matlab.ui.Figure
        GridLayout        matlab.ui.container.GridLayout
        LeftPanel         matlab.ui.container.Panel
        RightPanel        matlab.ui.container.Panel

        %% --- RIGHT PANEL: PLOTTING AND LOGGING ---
        PlotTabGroup      matlab.ui.container.TabGroup
        LogTextAreaLabel  matlab.ui.control.Label
        LogTextArea       matlab.ui.control.HTML

        %% --- LEFT PANEL: CONTROL BUTTONS ---
        RunBEMTAnalysisButton  matlab.ui.control.Button
        StopAnalysisButton matlab.ui.control.Button
        ClearPlotsButton  matlab.ui.control.Button
        CloseAllPlotsButton matlab.ui.control.Button

        %% --- LEFT PANEL: SETTINGS PANELS ---
        AnalysisSettingsPanel matlab.ui.container.Panel
        GeneralSettingsPanel matlab.ui.container.Panel
        XFOILSettingsPanel matlab.ui.container.Panel
        NeuralFoilSettingsPanel matlab.ui.container.Panel
        AnalysisRangePanel matlab.ui.container.Panel
        CorrectionsPanel        matlab.ui.container.Panel
        RunOptionsPanel         matlab.ui.container.Panel

        %% --- SETTINGS: ANALYSIS & PROPELLER ---
        PropellerDropDown matlab.ui.control.DropDown
        AnalysisTypeDropDown matlab.ui.control.DropDown
        
        %% --- SETTINGS: GENERAL ---
        HeightmEditField  matlab.ui.control.NumericEditField
        NumberofBladesBEditField matlab.ui.control.NumericEditField
        RPMEditField      matlab.ui.control.NumericEditField
        DpsiEditField     matlab.ui.control.NumericEditField
        
        %% --- SETTINGS: XFOIL ---
        NCritEditField    matlab.ui.control.NumericEditField
        FilterPassesEditField matlab.ui.control.NumericEditField
        UpperAoAEditField matlab.ui.control.NumericEditField
        LowerAoAEditField matlab.ui.control.NumericEditField
        MaxIterationsEditField matlab.ui.control.NumericEditField
        PlotXFOILPolarsCheckBox matlab.ui.control.CheckBox

        %% --- SETTINGS: NEURALFOIL ---
        MinReEditField matlab.ui.control.NumericEditField
        MaxReEditField matlab.ui.control.NumericEditField
        MinAoAEditField matlab.ui.control.NumericEditField
        MaxAoAEditField matlab.ui.control.NumericEditField
        ModelSizeDropDown       matlab.ui.control.DropDown

        %% --- SETTINGS: DYNAMIC ANALYSIS RANGE PANELS ---
        AlphaRangePanel    matlab.ui.container.Panel
        JPrimeRangePanel   matlab.ui.container.Panel
        VelocityRangePanel matlab.ui.container.Panel
        MinAlphaEditField    matlab.ui.control.NumericEditField
        AlphaStepEditField   matlab.ui.control.NumericEditField
        MaxAlphaEditField    matlab.ui.control.NumericEditField
        MinJEditField      matlab.ui.control.NumericEditField
        JStepEditField     matlab.ui.control.NumericEditField
        MaxJEditField      matlab.ui.control.NumericEditField
        MinVEditField      matlab.ui.control.NumericEditField
        VStepEditField     matlab.ui.control.NumericEditField
        MaxVEditField      matlab.ui.control.NumericEditField

        %% --- SETTINGS: CORRECTIONS & FIDELITY ---
        UseOldPolarDataCheckBox matlab.ui.control.CheckBox
        RadialFlowCheckBox      matlab.ui.control.CheckBox
        Use3DPolarCheckBox      matlab.ui.control.CheckBox
        FidelityModelDropDown   matlab.ui.control.DropDown
        PlotAzimuthalContoursCheckBox matlab.ui.control.CheckBox

        %% --- SETTINGS: RUN OPTIONS ---
        AnalysisMethodDropDown  matlab.ui.control.DropDown
        PolarDataSourceDropDown matlab.ui.control.DropDown
        
        %% --- INTERNAL STATE PROPERTIES ---
        StopRequested     matlab.lang.OnOffSwitchState % Flag to signal the analysis loop to stop.
        RunCounter        double = 0                   % Counter for naming plot tabs.
    end

    %% ===================================================================
    %% APP CONSTRUCTION AND DESTRUCTION
    %% ===================================================================
    methods (Access = public)
        % --- CONSTRUCTOR ---
        % Initializes the app, creates components, and runs startup function.
        function app = ExBEMT_GUI, clc, createComponents(app), startupFcn(app), if nargout == 0, clear app, end, end
        
        % --- DESTRUCTOR ---
        % Called when the app is closed to delete the UI figure.
        function delete(app), delete(app.UIFigure), end
    end

    %% ===================================================================
    %% PRIVATE METHODS (CALLBACKS AND HELPERS)
    %% ===================================================================
    methods (Access = private)

        % --- STARTUP FUNCTION ---
        % Executes once after the components are created.
        function startupFcn(app)
            % --- Initial Log Messages ---
            loghtml(app, 'App starting up...', 'info');
            
            % --- Populate UI Elements and Set Initial State ---
            populatePropellers(app);        % Find available propellers.
            analysisTypeChanged(app);       % Set visibility of range panels.
            analysisMethodChanged(app);     % Set visibility of settings panels.
            oldDataCheckChanged(app);       % Set visibility based on using old data.
            loghtml(app, 'App ready. Please select settings and run analysis!', 'success');
            
            % --- Display Startup Icon ---
            try
                fid = fopen('ExBEMT_Start.png', 'rb'); bytes = fread(fid, Inf, '*uint8'); fclose(fid);
                base64_string = matlab.net.base64encode(bytes);
                img_html = sprintf('<div style="text-align:center;"><img src="data:image/png;base64,%s" width="250"></div>', base64_string);
                loghtml(app, img_html, 'info');
            catch
                loghtml(app, 'Could not load and display the icon in the log.', 'warning');
            end
            
            % --- Initial Warning ---
            loghtml(app, 'Q-Unsteady Inflow Solver can be Unstable in Low-Speed and High Incidences. Use standard Dynamic Inflow for stability!!', 'warning');
        end

        %% ====================================================================
        %% UI CALLBACK FUNCTIONS (To control visibility and state)
        %% ====================================================================

        % --- Called when RPM value is changed ---
        function rpmChanged(app, ~)
            loghtml(app, 'RPM changed. Re-run XFOIL polar analysis for accuracy by unchecking "Use Existing Polars".', 'warning');
            updateDataSourceDropDown(app); % Update for RPM-dependent XFOIL files
        end
        
        % --- Called when a new propeller is selected ---
        function propellerChanged(app, ~)
            updateDataSourceDropDown(app); % Check for data for the newly selected propeller.
        end

        % --- Called when the analysis type (J', Velocity, etc.) is changed ---
        function analysisTypeChanged(app, ~)
            % --- Get selected analysis type ---
            type = app.AnalysisTypeDropDown.Value;
            isPrompt = strcmp(type, 'Select Analysis Type...');
            
            % --- Show/Hide the entire range panel ---
            app.AnalysisRangePanel.Visible = ~isPrompt;
            if isPrompt, return; end
            
            % --- Enable/Disable specific range sub-panels ---
            app.AlphaRangePanel.Enable = 'on';
            app.JPrimeRangePanel.Enable = 'on';
            app.VelocityRangePanel.Enable = 'on';
            if strcmp(type, 'Advance Ratio')
                app.VelocityRangePanel.Enable = 'off';
            elseif strcmp(type, 'Velocity')
                app.JPrimeRangePanel.Enable = 'off';
            elseif strcmp(type, 'Angle of Incidence')
                app.VelocityRangePanel.Enable = 'off';
            end
        end

        % --- Called when the "Use Existing Polars" checkbox is changed ---
        function oldDataCheckChanged(app, ~)
            useOld = app.UseOldPolarDataCheckBox.Value;
            
            % --- Toggle visibility and state of dependent controls ---
            app.PolarDataSourceDropDown.Visible = useOld;
            
            if useOld
                % --- Using old data: Disable polar generation settings ---
                app.AnalysisMethodDropDown.Enable = 'off';
                app.ModelSizeDropDown.Enable = 'off';
                app.XFOILSettingsPanel.Enable = 'off';
                app.NeuralFoilSettingsPanel.Enable = 'off';
                updateDataSourceDropDown(app); % Find available data files.
            else
                % --- Generating new data: Enable polar generation settings ---
                app.AnalysisMethodDropDown.Enable = 'on';
                analysisMethodChanged(app); % Update visibility based on method.
            end
        end

        % --- Called when the polar analysis method (XFOIL/NeuralFoil) is changed ---
        function analysisMethodChanged(app, ~)
            is_neuralfoil = strcmp(app.AnalysisMethodDropDown.Value, 'NeuralFoil');
            
            % --- Enable the correct settings panel if generating new data ---
            if ~app.UseOldPolarDataCheckBox.Value
                app.XFOILSettingsPanel.Enable = ~is_neuralfoil;
                app.NeuralFoilSettingsPanel.Enable = is_neuralfoil;
                app.ModelSizeDropDown.Enable = is_neuralfoil;
            end
        end

        % --- Populates the "Data Source" dropdown with available .mat or .xlsx files ---
        function updateDataSourceDropDown(app)
            % Checks for existing data files and populates the dropdown.
            if ~app.UseOldPolarDataCheckBox.Value || strcmp(app.PropellerDropDown.Value, 'Select Propeller...')
                app.PolarDataSourceDropDown.Items = {'(N/A)'};
                app.PolarDataSourceDropDown.Value = '(N/A)';
                return;
            end

            % --- Get propeller directory and RPM ---
            propName = app.PropellerDropDown.Value;
            rpm = app.RPMEditField.Value;
            [app_path, ~, ~] = fileparts(mfilename('fullpath'));
            dir_prop = fullfile(app_path, 'Propeller Packages', propName);
            
            % --- Find matching data files ---
            items = {};
            % --- Check for XFOIL data (RPM dependent) ---
            xfoil_cl_filename = sprintf('ClvsAlpha_RPM%d.xlsx', rpm);
            if exist(fullfile(dir_prop, xfoil_cl_filename), 'file')
                items{end+1} = sprintf('XFOIL Polars (RPM %d)', rpm);
            end

            % --- Check for NeuralFoil data (RPM independent) by pattern ---
            nf_files = dir(fullfile(dir_prop, 'NeuralFoil_Grid_*.mat'));
            for k = 1:length(nf_files)
                items{end+1} = nf_files(k).name; % Add the actual filename to the list
            end
            
            % --- Update dropdown items ---
            if isempty(items)
                app.PolarDataSourceDropDown.Items = {'No Data Found'};
                app.PolarDataSourceDropDown.Value = 'No Data Found';
            else
                app.PolarDataSourceDropDown.Items = items;
                app.PolarDataSourceDropDown.Value = items{1};
            end
        end

        %% ====================================================================
        %% EXECUTION AND PLOTTING CONTROL
        %% ====================================================================

        % --- Callback for the "Stop" button ---
        function stopAnalysisButtonPushed(app, ~)
            app.StopRequested = 'on';
            loghtml(app, 'STOP REQUESTED. Analysis will terminate at the next check point.', 'warning');
            app.StopAnalysisButton.Enable = 'off';
            % Create a flag file to signal parallel workers to stop
            fclose(fopen('cancel_analysis.flag', 'w'));
        end

        % --- Callback for the "Clear Plots" button ---
        function clearPlotsButtonPushed(app, ~)
            selectedTab = app.PlotTabGroup.SelectedTab;
            if ~isempty(selectedTab)
                tabTitle = selectedTab.Title;
                delete(selectedTab);
                loghtml(app, ['Plot tab "' tabTitle '" cleared.'], 'info');
            else
                loghtml(app, 'No plot tab to clear.', 'info');
            end
        end

        % --- Callback for the "Close All Plots" button ---
        function closeAllPlotsButtonPushed(app, ~)
            close all; 
            loghtml(app, 'All plot windows closed.', 'info');
        end

       % --- MAIN ANALYSIS FUNCTION: Callback for the "Run Analysis" button ---
       function runBEMTAnalysisButtonPushed(app, event)
            %% 0. PRE-ANALYSIS VALIDATION AND SETUP
            % ==========================================================
            
            % --- Input Validation ---
            if strcmp(app.PropellerDropDown.Value, 'Select Propeller...'), msg = 'Please select a propeller before running.'; uialert(app.UIFigure, msg, 'Invalid Input'); loghtml(app, msg, 'error'); return; end
            if strcmp(app.AnalysisTypeDropDown.Value, 'Select Analysis Type...'), msg = 'Please select an analysis type before running.'; uialert(app.UIFigure, msg, 'Invalid Input'); loghtml(app, msg, 'error'); return; end

            % --- Initialize Plotting Tab and UI State ---
            app.RunCounter = app.RunCounter + 1;
            analysisType = app.AnalysisTypeDropDown.Value;
            fidelityModel = app.FidelityModelDropDown.Value;
            runTitle = sprintf('Run %d: %s (%s)', app.RunCounter, app.PropellerDropDown.Value, analysisType);
            newTab = uitab(app.PlotTabGroup, 'Title', runTitle);
            app.PlotTabGroup.SelectedTab = newTab;

            % --- Create Axes within the new tab ---
            plotGridLayout = uigridlayout(newTab);
            plotGridLayout.ColumnWidth = {'1x', '1x'}; plotGridLayout.RowHeight = {'1x', '1x'};
            newAxes.CT = uiaxes(plotGridLayout); title(newAxes.CT, 'Thrust Coefficient');
            newAxes.CN = uiaxes(plotGridLayout); title(newAxes.CN, 'Normal Force Coefficient');
            newAxes.CQ = uiaxes(plotGridLayout); title(newAxes.CQ, 'Torque Coefficient');
            newAxes.Efficiency = uiaxes(plotGridLayout); title(newAxes.Efficiency, 'Propulsive Efficiency');
            hold(newAxes.CT, 'on'); hold(newAxes.CN, 'on'); hold(newAxes.CQ, 'on'); hold(newAxes.Efficiency, 'on');

            % --- Set UI Controls to "Running" State ---
            app.StopRequested = 'off';
            app.RunBEMTAnalysisButton.Enable = 'off';
            app.StopAnalysisButton.Enable = 'on';
            app.ClearPlotsButton.Enable = 'off';
            loghtml(app, ['Starting ExBEMT analysis with model: <b>' fidelityModel '</b>'], 'header');
            total_tic = tic;
            drawnow;

            % --- Cleanup Stop Flag from Previous Runs ---
            if exist('cancel_analysis.flag', 'file'), delete('cancel_analysis.flag'); end

            try
                %% 1. PARALLEL POOL MANAGEMENT
                % ==========================================================
                
                % --- Determine if NeuralFoil (process-based) is needed ---
                isNeuralFoil = strcmp(app.AnalysisMethodDropDown.Value, 'NeuralFoil') || ...
                               (app.UseOldPolarDataCheckBox.Value && contains(app.PolarDataSourceDropDown.Value, 'NeuralFoil'));
                current_pool = gcp('nocreate'); % Get current pool without creating one.

                % --- Manage Pool Type (Threads vs Processes) ---
                if isNeuralFoil
                    % NeuralFoil needs a process-based pool. Restart if current pool is thread-based.
                    if ~isempty(current_pool) && isa(current_pool, 'parallel.ThreadPool')
                        loghtml(app, 'NeuralFoil requires a process-based pool. Restarting existing thread-pool...', 'info');
                        delete(current_pool);
                        current_pool = []; 
                    end
                    if isempty(current_pool), parpool('Processes'); loghtml(app, 'Initializing process-based pool...', 'info'); end
                else
                    % XFOIL runs fine on a lighter thread-based pool.
                    if isempty(current_pool), parpool('Threads'); loghtml(app, 'Initializing thread-based pool...', 'info'); end
                end

                %% 2. GATHER USER INPUTS & CONSTRUCT SETTINGS STRUCT
                % ==========================================================
                loghtml(app, 'Gathering settings from GUI...', 'info');
                settings = struct();
                
                % --- General & Analysis Settings ---
                settings.Height = app.HeightmEditField.Value; 
                settings.Np = 1;
                settings.B = app.NumberofBladesBEditField.Value;
                settings.RPM = app.RPMEditField.Value; 
                settings.dpsi_deg = app.DpsiEditField.Value;
                settings.Sel_Prop = app.PropellerDropDown.Value;
                settings.analysisType = analysisType; 
                settings.fidelityModel = fidelityModel;
                
                % --- Run Options & Corrections ---
                settings.Use_old_data_val = app.UseOldPolarDataCheckBox.Value;
                settings.radial_correct = app.RadialFlowCheckBox.Value; 
                settings.Use_3D_polar = app.Use3DPolarCheckBox.Value;
                settings.use_neuralfoil = isNeuralFoil; 
                settings.neuralfoil_model_size = app.ModelSizeDropDown.Value;
                plotAzimuthal = app.PlotAzimuthalContoursCheckBox.Value;

                % --- XFOIL-Specific Settings ---
                NCrit_val = app.NCritEditField.Value; 
                filt = app.FilterPassesEditField.Value;
                UpperA_val = app.UpperAoAEditField.Value; 
                LowerA_val = app.LowerAoAEditField.Value;
                max_iter = app.MaxIterationsEditField.Value; 
                plot_polars = app.PlotXFOILPolarsCheckBox.Value;

                % --- NeuralFoil-Specific Settings ---
                settings.nf_re_min = app.MinReEditField.Value; 
                settings.nf_re_max = app.MaxReEditField.Value;
                settings.nf_aoa_min = app.MinAoAEditField.Value; 
                settings.nf_aoa_max = app.MaxAoAEditField.Value;

                % --- Analysis Range Vectors ---
                alpha_vector = app.MinAlphaEditField.Value:app.AlphaStepEditField.Value:app.MaxAlphaEditField.Value;
                j_vector = app.MinJEditField.Value:app.JStepEditField.Value:app.MaxJEditField.Value;
                v_vector = app.MinVEditField.Value:app.VStepEditField.Value:app.MaxVEditField.Value;

                %% 3. ATMOSPHERE & PROPELLER DATA LOADING
                % ==========================================================
                
                % --- Calculate Atmospheric Properties ---
                [~, ~, ~, rho, nu, a] = atmosisa(settings.Height);
                settings.rho = rho; 
                settings.nu = nu; 
                settings.omega = convangvel(settings.RPM, 'RPM', 'rad/s');
                loghtml(app, 'Atmosphere and propeller data loaded.', 'info');

                % --- Load Propeller Geometry File ---
                [app_path, ~, ~] = fileparts(mfilename('fullpath'));
                dir_prop = fullfile(app_path, 'Propeller Packages', settings.Sel_Prop);
                settings.dir_prop = dir_prop;
                prop_filepath = fullfile(dir_prop, 'Prop_sections.xlsx');
                if ~exist(prop_filepath, 'file'), error('Prop_sections.xlsx not found for propeller: %s', settings.Sel_Prop); end

                % --- Read and Parse Geometry Data from Excel File ---
                prop_data_cell = readcell(prop_filepath);
                c_row_data = prop_data_cell(3, 2:end); 
                is_valid_section = cellfun(@(x) isnumeric(x), c_row_data);
                num_sections = find(is_valid_section, 1, 'last'); 
                last_col_idx = num_sections + 1;

                % --- Store Geometry in Settings Struct ---
                settings.x_blade = cell2mat(prop_data_cell(2, 2:last_col_idx));
                settings.c = cell2mat(prop_data_cell(3, 2:last_col_idx));
                settings.beta = cell2mat(prop_data_cell(4, 2:last_col_idx));
                settings.delta_x = cell2mat(prop_data_cell(6, 2:last_col_idx));
                settings.R = prop_data_cell{8, 1};
                settings.foils_list = string(prop_data_cell(5, 2:last_col_idx));
                settings.foils_list(ismissing(settings.foils_list)) = [];
                
                % --- Import Python Module if using NeuralFoil ---
                if settings.use_neuralfoil
                    try
                        settings.py_nf_module = py.importlib.import_module('run_neuralfoil');
                    catch ME
                        loghtml(app, 'Error: Could not import Python module ''run_neuralfoil''.', 'error'); 
                        error(ME.message); 
                    end
                end

                %% 4. POLAR GENERATION OR DATA LOADING
                % ==========================================================
                CL_P = {}; CD_P = {}; AoA_data = {}; Cla = []; alpha_0 = [];

                if ~settings.Use_old_data_val
                   % --- Case: Generate New Polar Data ---
                   loghtml(app,'New analysis requested. Generating polars...', 'info');
                   
                   if settings.use_neuralfoil
                        % --- Generate Polars using NeuralFoil ---
                        loghtml(app, 'NeuralFoil selected. Pre-computing aerodynamic grids...', 'header');
                        installNeuralFoilDependencies(); % Helper function assumed to exist
                        
                        % --- Define AoA and Reynolds Number Grids ---
                        aoa_grid = settings.nf_aoa_min:0.5:settings.nf_aoa_max;
                        re_grid_log = logspace(log10(settings.nf_re_min), log10(settings.nf_re_max), 50);

                        % --- Initialize Interpolant Cells ---
                        CL_NF_Interp = cell(1, length(settings.foils_list)); 
                        CD_NF_Interp = cell(1, length(settings.foils_list));
                        unique_foils = unique(settings.foils_list);

                        % --- Loop through unique airfoils to generate data ---
                        for i = 1:length(unique_foils)
                            foil_name = unique_foils{i}; 
                            loghtml(app, ['Generating grid for airfoil: ', foil_name], 'info');
                            cl_data = zeros(length(aoa_grid), length(re_grid_log)); 
                            cd_data = zeros(length(aoa_grid), length(re_grid_log));
                            foil_path = fullfile(settings.dir_prop, 'foils', [foil_name '.dat']);
                            
                            % --- Call Python script for each Reynolds number ---
                            for re_idx = 1:length(re_grid_log)
                                re = re_grid_log(re_idx);
                                json_string = settings.py_nf_module.get_neuralfoil_aero(char(foil_path), py.list(aoa_grid), double(re), char(settings.neuralfoil_model_size));
                                aero_results_array = jsondecode(char(json_string));
                                if iscell(aero_results_array), aero_results_array = aero_results_array{1}; end
                                cl_data(:, re_idx) = [aero_results_array.CL]'; 
                                cd_data(:, re_idx) = [aero_results_array.CD]';
                            end
                            
                            % --- Create gridded interpolants for this airfoil ---
                            [AOA_GRID, RE_GRID] = ndgrid(aoa_grid, re_grid_log);
                            CL_NF_Interp{i} = griddedInterpolant(AOA_GRID, RE_GRID, cl_data, 'linear', 'nearest');
                            CD_NF_Interp{i} = griddedInterpolant(AOA_GRID, RE_GRID, cd_data, 'linear', 'nearest');
                        end
                        
                        % --- Map unique foil interpolants to the full list of blade sections ---
                        [~, foil_map_indices] = ismember(settings.foils_list, unique_foils);
                        settings.CL_NF_Interp = CL_NF_Interp(foil_map_indices);
                        settings.CD_NF_Interp = CD_NF_Interp(foil_map_indices);
                        loghtml(app, 'NeuralFoil pre-computation complete.', 'success');
                        
                        % --- Save the newly generated data to a .mat file ---
                        nf_mat_filename = sprintf('NeuralFoil_Grid_Re[%dk-%dk]_AoA[%d-%d].mat', ...
                                                  round(settings.nf_re_min/1000), round(settings.nf_re_max/1000), ...
                                                  settings.nf_aoa_min, settings.nf_aoa_max);
                        nf_mat_filepath = fullfile(dir_prop, nf_mat_filename);
                        save(nf_mat_filepath, 'CL_NF_Interp', 'CD_NF_Interp', 'aoa_grid', 're_grid_log');
                        loghtml(app, ['NeuralFoil grid data saved to: <b>' nf_mat_filename '</b>'], 'success');

                   else 
                        % --- Generate Polars using XFOIL ---
                        loghtml(app, 'Running XFOIL with extended polar model...', 'header');
                        loghtml(app, 'Calculating CD_max for each airfoil section...', 'info');
                        cd_max_array = zeros(1, length(settings.x_blade));
                        
                        % --- Calculate CDmax for Viterna model based on LE radius ---
                        for i = 1:length(settings.x_blade)
                           try
                               le_radius_ratio = calculate_LE_radius(fullfile(dir_prop, 'foils', [settings.foils_list{i} '.dat']));
                               cd_max_array(i) = 2.0772 - 3.978 * le_radius_ratio;
                           catch
                               cd_max_array(i) = 2.0; % Default value if LE parsing fails
                               loghtml(app, ['Could not process LE for ' settings.foils_list{i} '. Using default CD_max.'], 'warning');
                           end
                        end

                       % --- Calculate Re and Mach for each blade section ---
                       Re = (settings.rho * settings.omega .* settings.x_blade * settings.R .* settings.c) / settings.nu;
                       Mach = (settings.omega .* settings.x_blade * settings.R) / a;
                       xf = cell(1, length(settings.x_blade));
                       
                       % --- Run XFOIL for each section ---
                       for iy = 1:length(settings.x_blade)
                           if app.StopRequested, error('Analysis stopped by user.'); end
                           loghtml(app, ['Running XFOIL for section ', num2str(iy), ' (', settings.foils_list{iy}, ')...'], 'info');
                           [xf{iy}] = xfoil_call(settings.foils_list{iy}, max(Re(iy),2000), min(Mach(iy),0.295), NCrit_val, max_iter, LowerA_val, UpperA_val, filt, settings.Sel_Prop, dir_prop, plot_polars);
                       end
                       
                        % --- Blend XFOIL data with Viterna post-stall model ---
                        loghtml(app, 'Blending polar data using Viterna post-stall model...', 'info');
                        for i = 1:length(settings.x_blade)
                            cl_section_data = [xf{i}.Polars{1}.Alpha, xf{i}.Polars{1}.CL];
                            cd_section_data = [xf{i}.Polars{1}.Alpha, xf{i}.Polars{1}.CD];
                            [alpha_final, cl_final, cd_final] = app.generateBlendedPolarForSection(i, cl_section_data, cd_section_data, cd_max_array, settings.foils_list);
                            if isempty(alpha_final), error(['Could not generate blended polar for section ' num2str(i)]); end
                            
                            % --- Create interpolants and calculate Cla, alpha_0 ---
                            AoA_data{i} = alpha_final;
                            CL_P{i} = fit(alpha_final, cl_final, 'linearinterp');
                            CD_P{i} = fit(alpha_final, cd_final, 'linearinterp');
                            alpha_0(i) = fzero(CL_P{i}, 0);
                            slope = polyfit(deg2rad(linspace(alpha_0(i), alpha_0(i) + 2, 100)), CL_P{i}(linspace(alpha_0(i), alpha_0(i) + 2, 100)), 1);
                            Cla(i) = slope(1);
                        end
                        settings.Cla = Cla;
                        settings.alpha_0 = alpha_0;
                        
                        % --- Save newly generated XFOIL data to Excel files ---
                        if ~app.StopRequested
                            oldfolder = cd(dir_prop);
                            cl_filename = sprintf('ClvsAlpha_RPM%d.xlsx', settings.RPM);
                            cd_filename = sprintf('CdvsAlpha_RPM%d.xlsx', settings.RPM);
                            if exist(cl_filename, 'file'), delete(cl_filename); end
                            if exist(cd_filename, 'file'), delete(cd_filename); end
                            
                            % --- Prepare data matrix for writing ---
                            max_rows = 0;
                            for i = 1:length(settings.x_blade), max_rows = max(max_rows, length(AoA_data{i})); end
                    
                            num_cols = length(settings.x_blade) * 2;
                            cl_data_matrix = NaN(max_rows, num_cols); cd_data_matrix = NaN(max_rows, num_cols);
                            cl_var_names = cell(1, num_cols); cd_var_names = cell(1, num_cols);
                    
                            % --- Populate matrix and variable names ---
                            for i = 1:length(settings.x_blade)
                                aoa_vec = AoA_data{i}; cl_vec = CL_P{i}(aoa_vec); cd_vec = CD_P{i}(aoa_vec);
                                current_rows = length(aoa_vec);
                                col_idx1 = 2*i - 1; col_idx2 = 2*i;
                                cl_data_matrix(1:current_rows, col_idx1) = aoa_vec; cl_data_matrix(1:current_rows, col_idx2) = cl_vec;
                                cd_data_matrix(1:current_rows, col_idx1) = aoa_vec; cd_data_matrix(1:current_rows, col_idx2) = cd_vec;
                                foil_name = settings.foils_list{i};
                                cl_var_names{col_idx1} = ['Alpha_Sec' num2str(i) '_' foil_name]; cl_var_names{col_idx2} = ['CL_Sec' num2str(i) '_' foil_name];
                                cd_var_names{col_idx1} = ['Alpha_Sec' num2str(i) '_' foil_name]; cd_var_names{col_idx2} = ['CD_Sec' num2str(i) '_' foil_name];
                            end
                    
                            % --- Write tables to .xlsx files ---
                            final_cl_table = array2table(cl_data_matrix, 'VariableNames', cl_var_names);
                            final_cd_table = array2table(cd_data_matrix, 'VariableNames', cd_var_names);
                            writetable(final_cl_table, cl_filename);
                            writetable(final_cd_table, cd_filename);
                            cd(oldfolder);
                            loghtml(app, 'Polar analysis complete and data saved.', 'success');
                        end
                   end
                   
                else 
                   % --- Case: Use Existing Polar Data ---
                   dataSource = app.PolarDataSourceDropDown.Value;
                   loghtml(app,['Loading existing polar data from: <b>' dataSource '</b>'], 'info');

                   if contains(dataSource, 'NeuralFoil')
                       % --- Load from NeuralFoil .mat file ---
                       nf_mat_filepath = fullfile(dir_prop, dataSource);
                       if ~exist(nf_mat_filepath, 'file'), error(['NeuralFoil .mat file not found: ' dataSource]); end
                       
                       % --- Load data and update settings to match file's ranges ---
                       loaded_data = load(nf_mat_filepath, 'CL_NF_Interp', 'CD_NF_Interp', 'aoa_grid', 're_grid_log');
                       
                       CL_NF_Interp = loaded_data.CL_NF_Interp; 
                       CD_NF_Interp = loaded_data.CD_NF_Interp;
                       
                       settings.nf_aoa_min = min(loaded_data.aoa_grid);
                       settings.nf_aoa_max = max(loaded_data.aoa_grid);
                       settings.nf_re_min = min(loaded_data.re_grid_log);
                       settings.nf_re_max = max(loaded_data.re_grid_log);
                       
                       % --- Map unique foil interpolants to full blade section list ---
                       unique_foils = unique(settings.foils_list);
                       [~, foil_map_indices] = ismember(settings.foils_list, unique_foils);
                       settings.CL_NF_Interp = CL_NF_Interp(foil_map_indices);
                       settings.CD_NF_Interp = CD_NF_Interp(foil_map_indices);

                       settings.use_neuralfoil = true;
                       loghtml(app, 'Successfully loaded NeuralFoil grid and applied its settings.', 'success');

                   elseif contains(dataSource, 'XFOIL')
                       % --- Load from XFOIL .xlsx files ---
                       cl_filename = sprintf('ClvsAlpha_RPM%d.xlsx', settings.RPM); 
                       cd_filename = sprintf('CdvsAlpha_RPM%d.xlsx', settings.RPM);
                       cl_filepath = fullfile(dir_prop, cl_filename); 
                       cd_filepath = fullfile(dir_prop, cd_filename);
                       if ~exist(cl_filepath, 'file') || ~exist(cd_filepath, 'file'), error('XFOIL .xlsx files not found!'); end

                       % --- Read tables and create interpolants ---
                       cl_table = readtable(cl_filepath); 
                       cd_table = readtable(cd_filepath);
                       for i = 1:length(settings.x_blade)
                           CL_data_file = rmmissing([cl_table.(2*i-1), cl_table.(2*i)]);
                           CD_data_file = rmmissing([cd_table.(2*i-1), cd_table.(2*i)]);
                           AoA_data{i} = CL_data_file(:,1);
                           CL_P{i} = fit(CL_data_file(:,1), CL_data_file(:,2), 'linearinterp');
                           CD_P{i} = fit(CD_data_file(:,1), CD_data_file(:,2), 'linearinterp');
                           
                           % --- Calculate Cla and alpha_0 from loaded data ---
                           alpha_0(i) = fzero(CL_P{i},0);
                           slope = polyfit(deg2rad(linspace(alpha_0(i),alpha_0(i)+2,100)), CL_P{i}(linspace(alpha_0(i),alpha_0(i)+2,100)),1);
                           Cla(i) = slope(1);
                       end
                       settings.use_neuralfoil = false;
                       loghtml(app, 'Successfully loaded XFOIL data from .xlsx files.', 'success');
                   else
                       error('No valid data source selected or found.');
                   end
                end
                
                % --- Final assignments to settings struct ---
                settings.Cla = Cla; 
                settings.alpha_0 = alpha_0;
                
                % --- Check for stop request before starting main loop ---
                if app.StopRequested, error('Analysis stopped by user.'); end

                %% 5. ExBEMT MAIN LOOP
                % ==========================================================
                loghtml(app, 'ExBEMT Subroutine is Initiated...', 'header');

                % --- Determine Outer and Inner Loop Variables based on Analysis Type ---
                if strcmp(analysisType, 'Advance Ratio'), outer_loop_vector = alpha_vector; inner_loop_vector = j_vector;
                elseif strcmp(analysisType, 'Velocity'), outer_loop_vector = alpha_vector; inner_loop_vector = v_vector;
                else, outer_loop_vector = j_vector; inner_loop_vector = alpha_vector;
                end
                colors = lines(length(outer_loop_vector));

                % --- Pre-allocate Results Matrices ---
                C_T_all = nan(length(outer_loop_vector), length(inner_loop_vector));
                C_N_all = nan(length(outer_loop_vector), length(inner_loop_vector));
                C_Q_all = nan(length(outer_loop_vector), length(inner_loop_vector));
                J_eff_all = nan(length(outer_loop_vector), length(inner_loop_vector));
                
                % --- Setup for Parallel Progress Logging and Animated Plots ---
                n = settings.RPM/60;
                progressQueue = parallel.pool.DataQueue;
                afterEach(progressQueue, @(msg) loghtml(app, msg, 'info'));
                animated_lines = struct();
                
                % --- OUTER LOOP (Iterates through series, e.g., different AoA) ---
                for i_outer = 1:length(outer_loop_vector)
                    if exist('cancel_analysis.flag', 'file'), continue; end
                    
                    % --- Create Animated Line for this Series ---
                    if strcmp(analysisType,'Velocity'), legend_prefix='\alpha = '; legend_suffix=' deg';
                    elseif strcmp(analysisType,'Advance Ratio'), legend_prefix='\alpha = '; legend_suffix=' deg';
                    else, legend_prefix='J'' = '; legend_suffix='';
                    end
                    lgd_text = sprintf('%s%.1f%s', legend_prefix, outer_loop_vector(i_outer), legend_suffix);

                    animated_lines(i_outer).CT = animatedline(newAxes.CT, 'LineStyle','--','Marker','o','Color',colors(i_outer,:),'LineWidth',2,'DisplayName',lgd_text);
                    animated_lines(i_outer).CN = animatedline(newAxes.CN, 'LineStyle','--','Marker','o','Color',colors(i_outer,:),'LineWidth',2);
                    animated_lines(i_outer).CQ = animatedline(newAxes.CQ, 'LineStyle','--','Marker','o','Color',colors(i_outer,:),'LineWidth',2);
                    animated_lines(i_outer).Efficiency = animatedline(newAxes.Efficiency, 'LineStyle','--','Marker','o','Color',colors(i_outer,:),'LineWidth',2);

                    % --- INNER LOOP (Iterates through points in a series, e.g., different J') ---
                    for i_inner = 1:length(inner_loop_vector)
                        % --- Calculate current AoA, J', and V_inf for this iteration ---
                        if strcmp(analysisType, 'Angle of Incidence')
                            current_alpha = inner_loop_vector(i_inner); 
                            current_j = outer_loop_vector(i_outer);
                            current_V_inf = current_j * (settings.omega * settings.R) / (pi * cosd(current_alpha));
                        elseif strcmp(analysisType, 'Advance Ratio')
                            current_alpha = outer_loop_vector(i_outer); 
                            current_j = inner_loop_vector(i_inner);
                            current_V_inf = current_j * (settings.omega * settings.R) / (pi * cosd(current_alpha));
                        else % Velocity
                            current_alpha = outer_loop_vector(i_outer); 
                            current_V_inf = inner_loop_vector(i_inner);
                            current_j = pi * current_V_inf * cosd(current_alpha) / (settings.omega * settings.R);
                        end
                        V_inf_prll = current_V_inf * sind(current_alpha); 
                        V_inf_perp = current_V_inf * cosd(current_alpha);
                        
                        % --- CORE CALCULATION: Call the BEMT solver ---
                        [T, N_force, Q_torque, contourPlotData] = BEMT_inflow_Averaged(V_inf_prll, V_inf_perp, settings, CL_P, CD_P);
                        
                        % --- Log iteration completion ---
                        loghtml(app, sprintf('BEMT inflow iteration complete for AOI = %.1f deg, J'' = %.2f.', current_alpha, current_j), 'info');

                        % --- Store Results and Convert to Coefficients ---
                        J_eff_all(i_outer, i_inner) = current_j;
                        C_T_all(i_outer, i_inner) =  T/(settings.rho*(n)^2*(2*settings.R)^4);
                        C_N_all(i_outer, i_inner) =  N_force/(settings.rho*(n)^2*(2*settings.R)^4);
                        C_Q_all(i_outer, i_inner) =  Q_torque/(settings.rho*(n)^2*(2*settings.R)^5);

                        % --- Add new point to animated plots ---
                        x_val = inner_loop_vector(i_inner);
                        C_T_point = C_T_all(i_outer, i_inner);
                        if C_T_point <= 1e-2, break; end % Stop if thrust is going below "0" aka windmilling
                        
                        % --- Plot contours if requested ---
                        if plotAzimuthal, plotContourResults(app, fidelityModel, current_alpha, current_j, contourPlotData, settings); end

                        % --- Calculate and plot efficiency ---
                        eta_point = (C_T_point * current_j) / (C_Q_all(i_outer, i_inner) * 2 * pi);
                        if eta_point < 0 || eta_point > 1, eta_point = NaN; end % Exclude invalid efficiency values
                        
                        % --- Add points to live plots ---
                        addpoints(animated_lines(i_outer).CT, x_val, C_T_point);
                        addpoints(animated_lines(i_outer).CN, x_val, C_N_all(i_outer, i_inner));
                        addpoints(animated_lines(i_outer).CQ, x_val, C_Q_all(i_outer, i_inner));
                        addpoints(animated_lines(i_outer).Efficiency, x_val, eta_point);
                        drawnow('limitrate'); % Update plots without overwhelming the system
                    if app.StopRequested,error('Analysis stopped by user.'); end
                    end % End of inner loop
                    
                    % --- Log completion of the series ---
                    if strcmp(analysisType, 'Angle of Incidence'), log_msg = sprintf('--- Series J''=<b>%.2f</b> complete. ---', outer_loop_vector(i_outer));
                    else, log_msg = sprintf('--- Series AOI=<b>%.1f deg</b> complete. ---', outer_loop_vector(i_outer));
                    end
                    send(progressQueue, log_msg); % Send to parallel queue for safe logging
                end % End of outer loop

                %% 6. POST-ANALYSIS FINALIZATION
                % ==========================================================
                loghtml(app, 'Computation finished. Finalizing plots and logging data...', 'success');
                if app.StopRequested, error('Analysis stopped by user post-computation.'); end
                
                % --- Finalize plot appearance ---
                axesHandles.CT = newAxes.CT; 
                axesHandles.CN = newAxes.CN;
                axesHandles.CQ = newAxes.CQ; 
                axesHandles.Efficiency = newAxes.Efficiency;
                finalizePlots(app, axesHandles);
                
                % --- Log numerical results to the text area ---
                logPlotData(app, 'RESULTS SUMMARY', 'header', [], [], []);
                logPlotData(app, 'Thrust Coefficient (CT)', analysisType, inner_loop_vector, outer_loop_vector, C_T_all);
                logPlotData(app, 'Normal Force Coefficient (CN)', analysisType, inner_loop_vector, outer_loop_vector, C_N_all);
                logPlotData(app, 'Torque Coefficient (CQ)', analysisType, inner_loop_vector, outer_loop_vector, C_Q_all);
                eta_all = (C_T_all .* J_eff_all) ./ (C_Q_all * 2 * pi);
                logPlotData(app, 'Propulsive Efficiency (eta)', analysisType, inner_loop_vector, outer_loop_vector, eta_all);

                % --- Log total run time and final message ---
                total_time = toc(total_tic);
                loghtml(app, sprintf('====== TOTAL RUN TIME: %.2f seconds ======', total_time), 'header');
                loghtml(app, 'Analysis finished.', 'success');

            catch E
                % --- GENERIC ERROR HANDLING ---
                loghtml(app, ['<b>SCRIPT ERROR:</b> ' E.message], 'error');
                for i_err = 1:length(E.stack)
                    err_loc = ['<span style="padding-left: 20px;">In ', E.stack(i_err).name, ' at line ', num2str(E.stack(i_err).line), '</span>']; 
                    loghtml(app, err_loc, 'error'); 
                end
            end
            
            % --- Final UI State Cleanup ---
            if exist('cancel_analysis.flag', 'file'), delete('cancel_analysis.flag'); end
            app.RunBEMTAnalysisButton.Enable = 'on'; 
            app.StopAnalysisButton.Enable = 'off'; 
            app.ClearPlotsButton.Enable = 'on';
        end
        %% ====================================================================
        %% GUI CREATION (Generated by App Designer)
        %% ====================================================================
        function createComponents(app)
            % --- Main Figure and Layouts ---
            app.UIFigure=uifigure('Visible','off','Icon','ExBEMT_Icon.png','Position',[100 100 1350 950],'Name','ExBEMT Analysis GUI V3');
            app.GridLayout=uigridlayout(app.UIFigure,'ColumnWidth',{360,'1x'},'RowHeight',{'1x'});
            app.LeftPanel=uipanel(app.GridLayout,'Title','Controls');
            leftPanelLayout=uigridlayout(app.LeftPanel,'RowHeight',{'fit','fit','fit','fit','fit','fit','fit','fit'},'ColumnWidth',{'1x'}, 'Scrollable', 'on');
            app.RightPanel=uipanel(app.GridLayout,'Title','Output');
            rightPanelLayout=uigridlayout(app.RightPanel,'RowHeight',{'1x',22,170},'ColumnWidth',{'1x'});
            
            % --- Right Panel Components (Plots and Log) ---
            app.PlotTabGroup = uitabgroup(rightPanelLayout); 
            app.PlotTabGroup.Layout.Row = 1;
            app.LogTextAreaLabel = uilabel(rightPanelLayout, 'Text', 'Logging:'); 
            app.LogTextAreaLabel.Layout.Row = 2;
            app.LogTextArea = uihtml(rightPanelLayout); 
            app.LogTextArea.Layout.Row = 3;
            html_style='<style>body{font-family:Segoe UI,Arial,sans-serif;font-size:10pt}h2{color:#0072BD;margin-bottom:2px;margin-top:10px}hr{border:0;border-top:1px solid #7E2F8E;margin-top:2px;margin-bottom:10px}table{border-collapse:collapse;width:95%%}th,td{border:1px solid #ddd;text-align:left;padding:4px}th{background-color:#f2f2f2}</style>';
            app.LogTextArea.HTMLSource=['<html><head>' html_style '</head><body></body></html>'];
            
            % --- Left Panel Components (All controls) ---
            row = 1; % Row counter for stacking panels
            
            % --- General Settings Panel ---
            app.GeneralSettingsPanel = uipanel(leftPanelLayout, 'Title', 'General Settings');
            app.GeneralSettingsPanel.Layout.Row = row; row=row+1;
            gsPanelLayout=uigridlayout(app.GeneralSettingsPanel,'RowHeight',{'fit','fit'},'ColumnWidth',{'fit','2x','fit','1.5x'});
            uilabel(gsPanelLayout,'Text','Altitude (m)');
            app.HeightmEditField=uieditfield(gsPanelLayout,'numeric','Value',0);
            uilabel(gsPanelLayout,'Text','Blade Count (B)');
            app.NumberofBladesBEditField=uieditfield(gsPanelLayout,'numeric','Value',2);
            uilabel(gsPanelLayout,'Text','RPM');
            app.RPMEditField=uieditfield(gsPanelLayout,'numeric','Value',6000,'ValueDisplayFormat', '%.0f','ValueChangedFcn',createCallbackFcn(app,@rpmChanged,true));
            uilabel(gsPanelLayout,'Text','Azimuthal Step (deg)');
            app.DpsiEditField=uieditfield(gsPanelLayout,'numeric','Value',1, 'Limits', [0.1, 20], 'ValueDisplayFormat', '%.1f');

            % --- Analysis Settings Panel ---
            app.AnalysisSettingsPanel = uipanel(leftPanelLayout, 'Title', 'Analysis Settings');
            app.AnalysisSettingsPanel.Layout.Row = row; row=row+1;
            asPanelLayout=uigridlayout(app.AnalysisSettingsPanel,'RowHeight',{'fit','fit'},'ColumnWidth',{'fit','1x'});
            uilabel(asPanelLayout,'Text','Propeller');
            app.PropellerDropDown=uidropdown(asPanelLayout, 'ValueChangedFcn',createCallbackFcn(app,@propellerChanged,true));
            uilabel(asPanelLayout,'Text','Analysis Type (x-axis)');
            app.AnalysisTypeDropDown=uidropdown(asPanelLayout,'Items',{'Select Analysis Type...','Advance Ratio','Velocity','Angle of Incidence'},'Value','Select Analysis Type...','ValueChangedFcn',createCallbackFcn(app,@analysisTypeChanged,true)); 
            
            % --- Analysis Range Panel ---
            app.AnalysisRangePanel = uipanel(leftPanelLayout, 'Title', 'Analysis Range');
            app.AnalysisRangePanel.Layout.Row = row; row=row+1; app.AnalysisRangePanel.Visible = 'off';
            arPanelLayout = uigridlayout(app.AnalysisRangePanel);
            arPanelLayout.RowHeight = {'fit', 'fit', 'fit'}; arPanelLayout.ColumnWidth = {'1x'};
            app.AlphaRangePanel = uipanel(arPanelLayout); app.JPrimeRangePanel = uipanel(arPanelLayout);
            app.VelocityRangePanel = uipanel(arPanelLayout);
            app.AlphaRangePanel.Title = 'Angle of Incidence (deg)';
            alphaLayout = uigridlayout(app.AlphaRangePanel);
            alphaLayout.RowHeight = {'fit'}; alphaLayout.ColumnWidth = {'fit', '1x', 'fit', '1x', 'fit', '1x'};
            uilabel(alphaLayout, 'Text', 'Min'); app.MinAlphaEditField = uieditfield(alphaLayout, 'numeric', 'Value', 0);
            uilabel(alphaLayout, 'Text', 'Step'); app.AlphaStepEditField = uieditfield(alphaLayout, 'numeric', 'Value', 5);
            uilabel(alphaLayout, 'Text', 'Max'); app.MaxAlphaEditField = uieditfield(alphaLayout, 'numeric', 'Value', 40);
            app.JPrimeRangePanel.Title = 'Advance Ratio (J'') Range';
            jLayout = uigridlayout(app.JPrimeRangePanel);
            jLayout.RowHeight = {'fit'}; jLayout.ColumnWidth = {'fit', '1x', 'fit', '1x', 'fit', '1x'};
            uilabel(jLayout, 'Text', 'Min'); app.MinJEditField = uieditfield(jLayout,'numeric','Value',0.1);
            uilabel(jLayout, 'Text', 'Step'); app.JStepEditField = uieditfield(jLayout,'numeric','Value',0.05);
            uilabel(jLayout, 'Text', 'Max'); app.MaxJEditField = uieditfield(jLayout,'numeric','Value',0.7);
            app.VelocityRangePanel.Title = 'Velocity (m/s) Range';
            vLayout = uigridlayout(app.VelocityRangePanel);
            vLayout.RowHeight = {'fit'}; vLayout.ColumnWidth = {'fit', '1x', 'fit', '1x', 'fit', '1x'};
            uilabel(vLayout, 'Text', 'Min'); app.MinVEditField = uieditfield(vLayout,'numeric','Value',2);
            uilabel(vLayout, 'Text', 'Step'); app.VStepEditField = uieditfield(vLayout,'numeric','Value',2);
            uilabel(vLayout, 'Text', 'Max'); app.MaxVEditField = uieditfield(vLayout,'numeric','Value',10);

            % --- Corrections & Fidelity Panel ---
            app.CorrectionsPanel = uipanel(leftPanelLayout, 'Title', 'Corrections & Fidelity');
            app.CorrectionsPanel.Layout.Row = row; row=row+1;
            corrPanelLayout=uigridlayout(app.CorrectionsPanel,'RowHeight',{'fit','fit'},'ColumnWidth',{'1x','1x'});
            app.RadialFlowCheckBox=uicheckbox(corrPanelLayout,'Text','Radial Flow','Value',true);
            app.Use3DPolarCheckBox=uicheckbox(corrPanelLayout,'Text','3D Polar','Value',true); 
            app.Use3DPolarCheckBox.Layout.Column = 2;
            app.FidelityModelDropDown=uidropdown(corrPanelLayout,'Items',{'BEMT','BEMT+ Dynamic Inflow Model','BEMT+ Q-Unsteady Dynamic Inflow Model'},'Value','BEMT');
            app.FidelityModelDropDown.Layout.Row = 2; app.FidelityModelDropDown.Layout.Column = [1 2];
            
            % --- Run Options Panel ---
            app.RunOptionsPanel = uipanel(leftPanelLayout, 'Title', 'Run Options');
            app.RunOptionsPanel.Layout.Row = row; row=row+1;
            roPanelLayout=uigridlayout(app.RunOptionsPanel,'RowHeight',{'fit','fit','fit'},'ColumnWidth',{'fit','1x'});
            app.UseOldPolarDataCheckBox=uicheckbox(roPanelLayout,'Text','Use Existing Polars','Value',true,'ValueChangedFcn',createCallbackFcn(app,@oldDataCheckChanged,true));
            app.PlotAzimuthalContoursCheckBox=uicheckbox(roPanelLayout,'Text','Azimuthal Contours','Value',false); 
            app.PlotAzimuthalContoursCheckBox.Layout.Column = 2;
            uilabel(roPanelLayout,'Text','Analysis Method');
            app.AnalysisMethodDropDown=uidropdown(roPanelLayout,'Items',{'XFOIL','NeuralFoil'},'Value','XFOIL','ValueChangedFcn',createCallbackFcn(app,@analysisMethodChanged,true));
            uilabel(roPanelLayout,'Text','Data Source');
            app.PolarDataSourceDropDown=uidropdown(roPanelLayout,'Items',{'(N/A)'},'Value','(N/A)','Visible',true);
            
            % --- XFOIL Settings Panel ---
            app.XFOILSettingsPanel = uipanel(leftPanelLayout, 'Title', 'XFOIL Settings');
            app.XFOILSettingsPanel.Layout.Row = row; row=row+1;
            xfoilLayout=uigridlayout(app.XFOILSettingsPanel,'RowHeight',{'fit','fit','fit'},'ColumnWidth',{'fit','1x','fit','1x'});
            uilabel(xfoilLayout,'Text','Upper AoA');
            app.UpperAoAEditField=uieditfield(xfoilLayout,'numeric','Value',20);
            uilabel(xfoilLayout,'Text','Lower AoA');
            app.LowerAoAEditField=uieditfield(xfoilLayout,'numeric','Value',-15);
            uilabel(xfoilLayout,'Text','NCrit');
            app.NCritEditField=uieditfield(xfoilLayout,'numeric','Value',9);
            uilabel(xfoilLayout,'Text','Filter Passes');
            app.FilterPassesEditField=uieditfield(xfoilLayout,'numeric','Value',10);
            uilabel(xfoilLayout,'Text','Max Iterations');
            app.MaxIterationsEditField=uieditfield(xfoilLayout,'numeric','Value',250);
            app.PlotXFOILPolarsCheckBox=uicheckbox(xfoilLayout,'Text','Plot XFOIL Polars'); 
            app.PlotXFOILPolarsCheckBox.Layout.Column = [3 4];

            % --- NeuralFoil Settings Panel ---
            app.NeuralFoilSettingsPanel = uipanel(leftPanelLayout, 'Title', 'NeuralFoil Settings');
            app.NeuralFoilSettingsPanel.Layout.Row = row; row=row+1;
            nfLayout=uigridlayout(app.NeuralFoilSettingsPanel,'RowHeight',{'fit','fit','fit'},'ColumnWidth',{'fit','1x','fit','1x'});
            uilabel(nfLayout,'Text','Min AoA (deg)');
            app.MinAoAEditField=uieditfield(nfLayout,'numeric','Value',-90);
            uilabel(nfLayout,'Text','Max AoA (deg)');
            app.MaxAoAEditField=uieditfield(nfLayout,'numeric','Value',90);
            uilabel(nfLayout,'Text','Min Reynolds');
            app.MinReEditField=uieditfield(nfLayout,'numeric','Value',10000);
            uilabel(nfLayout,'Text','Max Reynolds');
            app.MaxReEditField=uieditfield(nfLayout,'numeric','Value',500000);
            uilabel(nfLayout,'Text','Model Size');
            app.ModelSizeDropDown=uidropdown(nfLayout,'Items',{'xxsmall','xsmall','small','medium','large','xlarge','xxlarge','xxxlarge'},'Value','large');
            app.ModelSizeDropDown.Layout.Column = [2 3];
            
            % --- Main Control Button Panel ---
            buttonPanel = uipanel(leftPanelLayout, 'BorderType', 'none');
            buttonPanel.Layout.Row = row; row=row+1;
            buttonLayout = uigridlayout(buttonPanel);
            buttonLayout.ColumnWidth = {'1.5x', '1x', '1.2x', '1.2x'}; buttonLayout.RowHeight = {'fit'};
            app.RunBEMTAnalysisButton = uibutton(buttonLayout, 'push', 'Text', 'Run Analysis');
            app.RunBEMTAnalysisButton.BackgroundColor = [0.07, 0.62, 1.00]; app.RunBEMTAnalysisButton.FontColor = [1 1 1];
            app.RunBEMTAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @runBEMTAnalysisButtonPushed, true);
            app.StopAnalysisButton = uibutton(buttonLayout, 'push', 'Text', 'Stop', 'Enable', 'off');
            app.StopAnalysisButton.BackgroundColor = [0.85, 0.33, 0.10]; app.StopAnalysisButton.FontColor = [1 1 1];
            app.StopAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @stopAnalysisButtonPushed, true);
            app.ClearPlotsButton = uibutton(buttonLayout, 'push', 'Text', 'Clear Plots');
            app.ClearPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @clearPlotsButtonPushed, true);
            app.CloseAllPlotsButton = uibutton(buttonLayout, 'push', 'Text', 'Close All');
            app.CloseAllPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @closeAllPlotsButtonPushed, true);           
            
            % --- Make the Figure Visible After All Components Are Created ---
            app.UIFigure.Visible = 'on';
        end
        
        %% ====================================================================
        %% HELPER AND UTILITY FUNCTIONS
        %% ====================================================================
        
        % --- Blends XFOIL data with Viterna post-stall model ---
        function [alpha_final, cl_final, cd_final] = generateBlendedPolarForSection(app, section_idx, cl_section_data, cd_section_data, cd_max_array, ~)
            % --- Data Cleaning ---
            cl_section_data = rmmissing(cl_section_data);
            cd_section_data = rmmissing(cd_section_data);
        
            % --- Check for valid input ---
            if isempty(cl_section_data) || isempty(cd_section_data)
                alpha_final=[]; cl_final=[]; cd_final=[];
                loghtml(app,['No valid XFOIL data for section ' num2str(section_idx) ' to blend.'],'warning');
                return;
            end
            
            % --- Extract XFOIL data and find stall points ---
            aoa_xfoil = cl_section_data(:,1);
            cl_xfoil_data = cl_section_data(:,2);
            cd_xfoil_data = cd_section_data(:,2);
            pos_indices = aoa_xfoil >= 0;
            [~, max_idx] = max(cl_xfoil_data(pos_indices));
            aoa_pos = aoa_xfoil(pos_indices);
            stall_pos_deg = aoa_pos(max_idx);
            cl_stall_pos = interp1(aoa_xfoil, cl_xfoil_data, stall_pos_deg);
            cd_stall_pos = interp1(aoa_xfoil, cd_xfoil_data, stall_pos_deg);
            neg_indices = aoa_xfoil < 0;
            [~, min_idx] = min(cl_xfoil_data(neg_indices));
            aoa_neg = aoa_xfoil(neg_indices);
            stall_neg_deg = aoa_neg(min_idx);
            cl_stall_neg = interp1(aoa_xfoil, cl_xfoil_data, stall_neg_deg);
            cd_stall_neg = interp1(aoa_xfoil, cd_xfoil_data, stall_neg_deg);
            cd_max = cd_max_array(section_idx);
            
            % --- Define blending functions (cosine-based) ---
            blend_start_pos = stall_pos_deg - 3;
            blend_end_pos = stall_pos_deg + 10;
            blend_fn_pos = @(a) (1 + cos(pi * (a - blend_start_pos) / (blend_end_pos - blend_start_pos))) / 2;
            blend_start_neg = stall_neg_deg + 3;
            blend_end_neg = stall_neg_deg - 10;
            blend_fn_neg = @(a) (1 + cos(pi * (a - blend_start_neg) / (blend_end_neg - blend_start_neg))) / 2;
            
            % --- Initialize output and smooth input data ---
            alpha_final = (-180:0.5:180)'; %180 degrees
            cl_final = zeros(size(alpha_final));
            cd_final = zeros(size(alpha_final));
            cl_xfoil_smoothed = smoothdata(cl_xfoil_data,'sgolay',5,'Degree',3);
            cd_xfoil_smoothed = smoothdata(cd_xfoil_data,'sgolay',5,'Degree',3);
            cl_xfoil_interp = fit(aoa_xfoil, cl_xfoil_smoothed, 'linearinterp');
            cd_xfoil_interp = fit(aoa_xfoil, cd_xfoil_smoothed, 'linearinterp');
            
            % --- Loop through all angles and apply blending logic ---
            for j = 1:length(alpha_final)
                a = alpha_final(j);
                [cl_sep, cd_sep] = extend_polar(deg2rad(a), deg2rad(stall_pos_deg), cl_stall_pos, cd_stall_pos, deg2rad(stall_neg_deg), cl_stall_neg, cd_stall_neg, cd_max);
                if a >= blend_start_pos && a <= blend_end_pos
                    % Positive stall region blending
                    f = blend_fn_pos(a);
                    cl_final(j) = f*cl_xfoil_interp(a) + (1-f)*cl_sep;
                    cd_final(j) = f*cd_xfoil_interp(a) + (1-f)*cd_sep;
                elseif a <= blend_start_neg && a >= blend_end_neg
                    % Negative stall region blending
                    f = blend_fn_neg(a);
                    cl_final(j) = f*cl_xfoil_interp(a) + (1-f)*cl_sep;
                    cd_final(j) = f*cd_xfoil_interp(a) + (1-f)*cd_sep;
                elseif a > blend_end_pos || a < blend_end_neg
                    % Fully post-stall region (use Viterna model only)
                    cl_final(j) = cl_sep;
                    cd_final(j) = cd_sep;
                else
                    % Pre-stall region (use XFOIL data only)
                    cl_final(j) = cl_xfoil_interp(a);
                    cd_final(j) = cd_xfoil_interp(a);
                end
            end
            
            % --- Final cleanup of any NaNs or Infs ---
            invalid_indices = isnan(cl_final) | isinf(cl_final) | isnan(cd_final) | isinf(cd_final);
            alpha_final(invalid_indices)=[];
            cl_final(invalid_indices)=[];
            cd_final(invalid_indices)=[];
        end

        % --- Appends a styled message to the HTML log area ---
        function loghtml(app, message, style)
            % --- Select HTML style based on message type ---
            switch style
                case 'error',html_line=sprintf('<p style="color:red;">%s</p>',message);
                case 'warning',html_line=sprintf('<p style="color:#D95319;"><b>Warning:</b> %s</p>',message); 
                case 'success',html_line=sprintf('<p style="color:green;">%s</p>',message);
                case 'header',html_line=sprintf('<h2>%s</h2><hr>',message); 
                case 'table',html_line=['<pre>' message '</pre>'];
                otherwise,html_line=sprintf('<p>%s</p>',message); 
            end
            
            % --- Javascript to auto-scroll the log to the bottom ---
            js_scroll='<script>window.scrollTo(0,document.body.scrollHeight);</script>';
            
            % --- Insert the new line and scroll script into the HTML ---
            current_html=app.LogTextArea.HTMLSource;
            new_html=strrep(current_html,'</body>',[html_line '</body>']);
            new_html=strrep(new_html,'</html>',[js_scroll '</html>']);
            app.LogTextArea.HTMLSource=new_html;
            
            % --- Refresh the UI to show the new message ---
            drawnow('limitrate');
        end
        
        % --- Finds propeller folders and populates the dropdown menu ---
        function populatePropellers(app)
            try
                propPackagesDir=fullfile(pwd,'Propeller Packages');
                loghtml(app,['Searching for propellers in: ' propPackagesDir],'info');
                
                % --- Check if directory exists ---
                if ~exist(propPackagesDir,'dir')
                    loghtml(app,'"Propeller Packages" folder not found.','error');
                    app.PropellerDropDown.Items={'Folder not found'};
                    return;
                end
                
                % --- Find all subdirectories ---
                files=dir(propPackagesDir);
                propNames_temp=files([files.isdir] & ~ismember({files.name},{'.','..'}));
                
                % --- Update dropdown items ---
                if isempty(propNames_temp)
                    app.PropellerDropDown.Items={'None found'};
                else
                    propNames={propNames_temp.name};
                    app.PropellerDropDown.Items=[{'Select Propeller...'},propNames];
                end
            catch ME
                loghtml(app,['Error populating propellers: ' ME.message],'error');
                app.PropellerDropDown.Items={'Error'};
            end
        end
        
        % --- Formats and logs a table of results to the HTML log area ---
        function logPlotData(app, titleStr, analysisType, x_vector, series_vector, y_matrix)
            try
                % --- Handle header-only case ---
                if isempty(x_vector)||isempty(series_vector)||isempty(y_matrix)
                    loghtml(app,titleStr,'header');
                    return;
                end
                
                % --- Define table labels based on analysis type ---
                if strcmp(analysisType,'Angle of Incidence'),x_label='AOI(deg)';series_label_prefix='J''=';
                else,x_label='J''/V';series_label_prefix='AOI=';end
                if ~ismatrix(y_matrix),error('Input data is not a 2D matrix.');end
                y_matrix_log=y_matrix'; % Transpose for logging

                % --- Define formats for header and data rows ---
                header_format = ['%-10s' repmat('%-10s', 1, length(series_vector)) '<br>'];
                row_format = ['%-10.2f' repmat('%-10.4f', 1, length(series_vector)) '<br>'];
                
                % --- Build header string ---
                header_labels = {x_label};
                for i=1:length(series_vector)
                    header_labels{end+1} = [series_label_prefix num2str(series_vector(i))];
                end
                
                % --- Assemble full table string ---
                table_str = sprintf(header_format, header_labels{:});
                table_str = [table_str, repmat('-', 1, 10*(length(series_vector)+1)), '<br>'];
                
                for i=1:length(x_vector)
                    row_data = [x_vector(i), y_matrix_log(i,:)];
                    table_str = [table_str, sprintf(row_format, row_data)];
                end
                
                % --- Send to log ---
                loghtml(app,['<h3>' titleStr '</h3>'],'info');
                loghtml(app,table_str,'table');
            catch ME_log
                loghtml(app,['<b>Could not generate results table for "' titleStr '".</b> Error: ' ME_log.message],'error');
            end
        end
        
        % --- Finalizes plot appearance (grids, legends, etc.) after run ---
        function finalizePlots(~, axesHandles)
            % --- Get all axes handles ---
            all_axes=[axesHandles.CT,axesHandles.CN,axesHandles.CQ,axesHandles.Efficiency];
            
            % --- Apply grid and box to all axes ---
            for ax=all_axes
                hold(ax,'off');
                grid(ax,'on');
                box(ax,'on');
            end
            
            % --- Create a single, shared legend ---
            lgd=legend(axesHandles.CT,'Location','northeastoutside');
            title(lgd,'Analysis Series');
        end
        
        % --- Generates and displays azimuthal contour plots ---
        function plotContourResults(app, fidelityModel, alpha, J_prime, plotData, settings)
            loghtml(app,'Generating contour plots...','info');
            try
                % --- Get necessary settings ---
                x_blade = settings.x_blade;
                B = settings.B;
                
                % --- Determine which variables to plot based on fidelity model ---
                if strcmp(fidelityModel,'BEMT')
                    var_list={'dT_BEMT','dQ_BEMT','Vel_BEMT','AoA_BEMT'};
                    subplot_dims=[1,4];
                else
                    var_list={'dT_BEMT','dQ_BEMT','Vel_BEMT','AoA_BEMT','dT_PITT','dQ_PITT','Vel_PITT','AoA_PITT'};
                    subplot_dims=[2,4];
                end
                
                % --- Create new figure for contour plots ---
                fig=figure('Name','Contour Plots','NumberTitle','off','Position',[50,50,1600,800]);
                sgtitle(fig,sprintf('Contour Plots for AOI=%.1f deg, J''=%.2f',alpha,J_prime),'FontSize',14,'FontWeight','bold');
                
                % --- Add vertical text label for the BEMT row ---
                annotation(fig, 'textbox', [0.07, 0.7, 0, 0], ...
                    'String', 'BEMT', ...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', 14, ...
                    'FontWeight', 'bold', ...
                    'Rotation', 90);

                % --- Conditionally add vertical text label for the Dynamic Inflow row ---
                if ~strcmp(fidelityModel,'BEMT')
                    annotation(fig, 'textbox', [0.07, 0.25, 0, 0], ...
                        'String', 'Dynamic Inflow', ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 14, ...
                        'FontWeight', 'bold', ...
                        'Rotation', 90);
                end

                % --- Loop through variables and create each subplot ---
                for i=1:length(var_list)
                    var_name=var_list{i};
                    if ~isfield(plotData, var_name), continue; end
                    
                    prop_plot_cell = plotData.(var_name);
                    if isempty(prop_plot_cell), continue; end

                    % --- Assemble data from all blades into a single matrix ---
                    num_azi_points = size(prop_plot_cell{1,1},2);
                    prop_plot_matrix = zeros(length(x_blade), (num_azi_points-1)*B);
                    for sec=1:length(x_blade)
                        full_rot_data=[];
                        for blade=1:B
                            blade_data=prop_plot_cell{sec,blade};
                            full_rot_data=[full_rot_data, blade_data(1:end-1)];
                        end
                        prop_plot_matrix(sec,:)=full_rot_data;
                    end
                    
                    % --- Make data periodic for seamless plotting ---
                    prop_plot_matrix_periodic = [prop_plot_matrix, prop_plot_matrix(:,1)];
                    psi_vector = linspace(0, 2*pi, size(prop_plot_matrix_periodic, 2));
                    
                    % --- Convert from polar to cartesian coordinates for plotting ---
                    [Psi, R_grid] = meshgrid(psi_vector, x_blade);
                    [X,Y] = pol2cart(Psi, R_grid);
                    
                    % --- Rotate 90 degrees for standard propeller view ---
                    X_rot = Y; Y_rot = -X;
                    
                    % --- Create subplot and plot data ---
                    ax=subplot(subplot_dims(1),subplot_dims(2),i);
                    contourf(ax,X_rot,Y_rot,prop_plot_matrix_periodic,20,'LineStyle','none');
                    colorbar(ax);
                    
                    % --- Set appropriate title ---
                    if contains(var_name,'dT'),t_str='Thrust (dT)';
                    elseif contains(var_name,'dQ'),t_str='Torque (dQ)';
                    elseif contains(var_name,'Vel'),t_str='Relative Velocity (m/s)';
                    else,t_str='Angle of Attack (deg)';
                    end
                    
                    % --- Finalize subplot appearance ---
                    title(ax,t_str);
                    if contains(var_name,'PITT'),title(ax,t_str,'Color','r');end
                    axis(ax,'equal','off');
                    drawnow
                end
            catch ME
                loghtml(app,['<b>Contour Plotting Error:</b> ' ME.message],'error');
            end
        end
    end
end