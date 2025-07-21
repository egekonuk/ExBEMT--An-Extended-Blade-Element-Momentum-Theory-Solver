classdef ExBEMT_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure          matlab.ui.Figure
        GridLayout        matlab.ui.container.GridLayout
        LeftPanel         matlab.ui.container.Panel
        RightPanel        matlab.ui.container.Panel
        LogTextAreaLabel  matlab.ui.control.Label
        LogTextArea       matlab.ui.control.HTML
        RunBEMTAnalysisButton  matlab.ui.control.Button
        StopAnalysisButton matlab.ui.control.Button
        ClearPlotsButton  matlab.ui.control.Button
        CloseAllPlotsButton matlab.ui.control.Button % New button
        AnalysisSettingsPanel matlab.ui.container.Panel
        PropellerDropDown matlab.ui.control.DropDown
        AnalysisTypeDropDown matlab.ui.control.DropDown
        GeneralSettingsPanel matlab.ui.container.Panel
        HeightmEditField  matlab.ui.control.NumericEditField
        NumberofBladesBEditField matlab.ui.control.NumericEditField
        RPMEditField      matlab.ui.control.NumericEditField
        XFOILSettingsPanel matlab.ui.container.Panel
        NCritEditField    matlab.ui.control.NumericEditField
        FilterPassesEditField matlab.ui.control.NumericEditField
        UpperAoAEditField matlab.ui.control.NumericEditField
        LowerAoAEditField matlab.ui.control.NumericEditField
        MaxIterationsEditField matlab.ui.control.NumericEditField % ADD THIS LINE
        PlotXFOILPolarsCheckBox matlab.ui.control.CheckBox       % ADD THIS LINE
        
        % Dynamic Range Panels
        AnalysisRangePanel matlab.ui.container.Panel
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
        
        % Corrections Panel
        CorrectionsPanel        matlab.ui.container.Panel
        UseOldPolarDataCheckBox matlab.ui.control.CheckBox
        RadialFlowCheckBox      matlab.ui.control.CheckBox
        Use3DPolarCheckBox      matlab.ui.control.CheckBox
        FidelityModelDropDown   matlab.ui.control.DropDown
        PlotAzimuthalContoursCheckBox matlab.ui.control.CheckBox

        % Run Options Panel
        RunOptionsPanel         matlab.ui.container.Panel

        % Plotting Area
        PlotTabGroup      matlab.ui.container.TabGroup
        
        % Control Properties
        StopRequested     matlab.lang.OnOffSwitchState
        RunCounter        double = 0
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ExBEMT_GUI
            clc
            % Create UIFigure and components
            createComponents(app)
            
            % Run startup logic after all components are created
            startupFcn(app);
            
            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)
            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
    
    % Callbacks and Helper methods
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            loghtml(app, 'App starting up...', 'info');
            populatePropellers(app);
            analysisTypeChanged(app); % Set initial visibility of range panels
            oldDataCheckChanged(app); % Set initial state of XFOIL panel
            loghtml(app, 'App ready. Please select settings and run analysis!', 'success');
             try
                % Read the image file into a byte array
                fid = fopen('ExBEMT_Start.png', 'rb');
                bytes = fread(fid, Inf, '*uint8');
                fclose(fid);
                
                % Convert to Base64 string
                base64_string = matlab.net.base64encode(bytes);
                
                % Create the HTML for the image, centered
                img_html = sprintf('<div style="text-align:center;"><img src="data:image/png;base64,%s" width="250"></div>', base64_string);
                loghtml(app, img_html, 'info');
            catch ME
                loghtml(app, 'Could not load and display the icon in the log.', 'warning');
             end
             %
            msg = 'Q-Unsteady Inflow Solver can be Unstable in slow and High Incidences! (Check warnings for in Command Window!)';
            loghtml(app, msg, 'warning');
        end
        
        % Callback for RPM field change
        function rpmChanged(app, ~)
            msg = 'RPM changed. Re-run XFOIL analysis for accurate polars by setting "Use Existing Polar Data" unchecked!!';
            loghtml(app, msg, 'warning');
        end
        
        % Callback for the Analysis Type dropdown
        function analysisTypeChanged(app, ~)
            type = app.AnalysisTypeDropDown.Value;
            
            isPrompt = strcmp(type, 'Select Analysis Type...');
            app.AnalysisRangePanel.Visible = ~isPrompt;
            if isPrompt, return; end
            
            % Enable all, then disable the unused one.
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
        
        % Callback for the "Use Old Polar Data" checkbox
        function oldDataCheckChanged(app, ~)
            if app.UseOldPolarDataCheckBox.Value
                % If using old data, disable XFOIL settings
                app.XFOILSettingsPanel.Enable = 'off';
            else
                % If not using old data, enable XFOIL settings
                app.XFOILSettingsPanel.Enable = 'on';
            end
        end
        
        % Callback for the Stop Button
        function stopAnalysisButtonPushed(app, ~)
            app.StopRequested = 'on';
            loghtml(app, 'STOP REQUESTED. Analysis will terminate at the next check point.', 'warning');
            app.StopAnalysisButton.Enable = 'off';
        end
        
        % Callback for the Clear Plots Button
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

        % New callback for the Close All Plots button
        function closeAllPlotsButtonPushed(app, ~)
            close all;
            loghtml(app, 'All plot windows closed.', 'info');
        end

        % Button pushed function: RunBEMTAnalysisButton
% Button pushed function: RunBEMTAnalysisButton
        function runBEMTAnalysisButtonPushed(app, event)
            % --- Input Validation ---
            if strcmp(app.PropellerDropDown.Value, 'Select Propeller...')
                msg = 'Please select a propeller before running.';
                loghtml(app, msg, 'error');
                uialert(app.UIFigure, msg, 'Invalid Input');
                return;
            end
            if strcmp(app.AnalysisTypeDropDown.Value, 'Select Analysis Type...')
                msg = 'Please select an analysis type before running.';
                loghtml(app, msg, 'error');
                uialert(app.UIFigure, msg, 'Invalid Input');
                return;
            end
            
            app.RunCounter = app.RunCounter + 1;
            analysisType = app.AnalysisTypeDropDown.Value;
            fidelityModel = app.FidelityModelDropDown.Value;

            % Create a single new tab for this entire run
            runTitle = sprintf('Run %d: %s (%s)', app.RunCounter, app.PropellerDropDown.Value, analysisType);
            newTab = uitab(app.PlotTabGroup, 'Title', runTitle);
            app.PlotTabGroup.SelectedTab = newTab;
            
            plotGridLayout = uigridlayout(newTab);
            plotGridLayout.ColumnWidth = {'1x', '1x'};
            plotGridLayout.RowHeight = {'1x', '1x'};

            newAxes.CT = uiaxes(plotGridLayout); title(newAxes.CT, 'Thrust Coefficient');
            newAxes.CN = uiaxes(plotGridLayout); title(newAxes.CN, 'Normal Force Coefficient');
            newAxes.CQ = uiaxes(plotGridLayout); title(newAxes.CQ, 'Torque Coefficient');
            newAxes.Efficiency = uiaxes(plotGridLayout); title(newAxes.Efficiency, 'Propulsive Efficiency');
            
            % Set all axes to hold on for live plotting
            hold(newAxes.CT, 'on'); hold(newAxes.CN, 'on');
            hold(newAxes.CQ, 'on'); hold(newAxes.Efficiency, 'on');

            % Reset stop flag and manage buttons
            app.StopRequested = 'off';
            app.RunBEMTAnalysisButton.Enable = 'off';
            app.StopAnalysisButton.Enable = 'on';
            app.ClearPlotsButton.Enable = 'off';
            
            loghtml(app, ['Starting ExBEMT analysis with model: <b>' fidelityModel '</b>'], 'header');
            total_tic = tic; % Start total timer
            drawnow;

            % Declare global variables used by the BEMT functions
            global c B R omega CL_P CD_P k_w Np x_blade delta_x rho RPM mu iy beta UpperA LowerA AoA_data radial_correct Use_3D_polar Sel_Prop
            
            try
                %% Get User Inputs from GUI
                plotData.Height = app.HeightmEditField.Value;
                Np = 1; % Hardcode Number of Propellers to 1
                plotData.Np = Np;
                B = app.NumberofBladesBEditField.Value;
                plotData.B = B;
                RPM = app.RPMEditField.Value;
                plotData.Sel_Prop = app.PropellerDropDown.Value;
                plotData.analysisType = analysisType;
                
                Use_old_data_val = app.UseOldPolarDataCheckBox.Value;
                radial_correct = app.RadialFlowCheckBox.Value;
                Use_3D_polar = app.Use3DPolarCheckBox.Value;
                
                % MODIFIED: Read new XFOIL settings from GUI
                if ~Use_old_data_val
                    NCrit_val = app.NCritEditField.Value;
                    filt = app.FilterPassesEditField.Value;
                    UpperA_val = app.UpperAoAEditField.Value;
                    LowerA_val = app.LowerAoAEditField.Value;
                    max_iter = app.MaxIterationsEditField.Value;
                    plot_polars = app.PlotXFOILPolarsCheckBox.Value;
                else
                    % Set default values if using old data (these won't be used anyway)
                    NCrit_val = 3; filt = 10; UpperA_val = 25; LowerA_val = -20;
                    max_iter = 200; plot_polars = true;
                end

                min_alpha = app.MinAlphaEditField.Value;
                step_alpha = app.AlphaStepEditField.Value;
                max_alpha = app.MaxAlphaEditField.Value;
                alpha_vector = min_alpha:step_alpha:max_alpha;
                
                min_j = app.MinJEditField.Value;
                step_j = app.JStepEditField.Value;
                max_j = app.MaxJEditField.Value;
                j_vector = min_j:step_j:max_j;
                
                min_v = app.MinVEditField.Value;
                step_v = app.VStepEditField.Value;
                max_v = app.MaxVEditField.Value;
                v_vector = min_v:step_v:max_v;
                
                k_w = 1; 
                
                %% Atmosphere Calculation
                [~, a, ~, rho, mu, ~] = atmosisa(plotData.Height);
                loghtml(app, 'Atmosphere calculated.', 'info');
                
                %% Propeller Data Loading
                [app_path, ~, ~] = fileparts(mfilename('fullpath'));
                dir_prop = fullfile(app_path, 'Propeller Packages', plotData.Sel_Prop);
                prop_filepath = fullfile(dir_prop, 'Prop_sections.xlsx');

                if ~exist(prop_filepath, 'file')
                    error('Prop_sections.xlsx not found in the selected propeller package directory.');
                end
                
                prop_data_cell = readcell(prop_filepath);

                c_row_data = prop_data_cell(3, 2:end);
                is_valid_section = cellfun(@(x) isnumeric(x) && ~isnan(x) && isreal(x), c_row_data);
                num_sections = find(is_valid_section, 1, 'last');

                if isempty(num_sections)
                    error('Could not determine the number of propeller sections. Check for valid numeric data in Row 3 of Prop_sections.xlsx.');
                end
                last_col_idx = num_sections + 1;
                
                x_blade = cell2mat(prop_data_cell(2, 2:last_col_idx));
                c       = cell2mat(prop_data_cell(3, 2:last_col_idx));
                beta    = cell2mat(prop_data_cell(4, 2:last_col_idx));
                delta_x = cell2mat(prop_data_cell(6, 2:last_col_idx));
                R       = prop_data_cell{8, 1};

                s = dir(fullfile(dir_prop, 'ClvsAlpha.xlsx'));
                if isempty(s) || s.bytes < 7000 || ~Use_old_data_val
                   loghtml(app,'Polar file not found or new analysis requested. Running XFOIL.', 'info');
                   Use_old_data = 2;
                else
                   loghtml(app,'Old CL,CD data found and will be used.', 'info');
                   Use_old_data = 1;
                end
                
                Foil_Use = 1; 
                n = RPM/60;
                omega = convangvel(RPM,'RPM','rad/s');
                Vel = omega * (x_blade * R);
                Re = (rho .* Vel .* c) ./ mu;
                Mach = Vel./a;
                loghtml(app, 'Propeller blade data loaded and conditions calculated.', 'info');

                %% XFOIL Calculations or Data Loading from File
                if Foil_Use && Use_old_data == 2
                    Foils = string(prop_data_cell(5, 2:last_col_idx));
                    Foils(ismissing(Foils)) = []; 

                    % REMOVED: Scaling algorithm for AoA and NCrit
                    
                    for iy = 1:length(x_blade)
                        if app.StopRequested, error('Analysis stopped by user.'); end
                        
                        loghtml(app, ['Running XFOIL for section ', num2str(iy), ' (', Foils{iy}, ')...'], 'info');
                        
                        % MODIFIED: Call xfoil_call with direct user inputs and new options
                        [xf{iy}] = xfoil_call(Foils{iy}, Re(iy), Mach(iy), NCrit_val, ...
                                              max_iter, round(LowerA_val), round(UpperA_val), ...
                                              filt, plotData.Sel_Prop, plot_polars);
                    end
                    
                    if ~app.StopRequested
                        oldfolder = cd(dir_prop);
                        recycle('on');
                        if exist('ClvsAlpha.xlsx', 'file'), delete('ClvsAlpha.xlsx'); end
                        if exist('CdvsAlpha.xlsx', 'file'), delete('CdvsAlpha.xlsx'); end

                        spec = length(x_blade)*2-2;
                        for j_xfoil = 0:2:spec
                            ind = round((j_xfoil+1)/2);
                            AoA_data{ind} = xf{1, ind}.Polars{1, 1}.Alpha;
                            CL_data{ind} = xf{1, ind}.Polars{1, 1}.CL;
                            CD_data{ind} = xf{1, ind}.Polars{1, 1}.CD;

                            if j_xfoil>24, asciiChars = char('A',(j_xfoil-26)+1+'A'-1,'2')'; asciiChars2 = char('A',(j_xfoil-26)+1+'A'-1,'1')';
                            else, asciiChars = char(j_xfoil+1+'A'-1,'2')'; asciiChars2 = char(j_xfoil+1+'A'-1,'1')';
                            end

                            writematrix(["Alpha" , "CL@" + ind], 'ClvsAlpha.xlsx', 'Range', asciiChars2);
                            writematrix(["Alpha" , "CD@" + ind], 'CdvsAlpha.xlsx', 'Range', asciiChars2);
                            writematrix([AoA_data{ind} CL_data{ind}], 'ClvsAlpha.xlsx', 'Range', asciiChars);
                            writematrix([AoA_data{ind} CD_data{ind}], 'CdvsAlpha.xlsx', 'Range', asciiChars);

                            CL_P{ind} = fit(AoA_data{ind},CL_data{ind},'linearinterp');
                            CD_P{ind} = fit(AoA_data{ind},CD_data{ind},'linearinterp');
                            alpha_0(ind) = fzero(CL_P{ind},0);
                            slope = polyfit(deg2rad(linspace(alpha_0(ind),alpha_0(ind)+2,100)), CL_P{ind}(linspace(alpha_0(ind),alpha_0(ind)+2,100)),1);
                            Cla(ind) = slope(1);
                        end
                        cd(oldfolder);
                        loghtml(app, 'XFOIL analysis complete and data saved.', 'success');
                    end
                else
                    spec = length(x_blade)*2-2;
                    j_read = 1;
                    for i_read = 0:2:spec
                        if i_read>24, asciiChars = char('A',(i_read-26)+1+'A'-1,':','A',(i_read-26)+2+'A'-1)';
                        else, asciiChars = char(i_read+1+'A'-1,':',i_read+2+'A'-1)';
                        end
                        CL_data_prop{j_read} = readmatrix(fullfile(dir_prop, 'ClvsAlpha.xlsx'),'Range',asciiChars);
                        CD_data_prop{j_read} = readmatrix(fullfile(dir_prop, 'CdvsAlpha.xlsx'),'Range',asciiChars);
                        CL_data_prop{j_read} = rmmissing(CL_data_prop{j_read});
                        CD_data_prop{j_read} = rmmissing(CD_data_prop{j_read});
                        
                        CL_data_file=CL_data_prop{j_read}; CD_data_file = CD_data_prop{j_read};
                        AoA_data{j_read} = CL_data_file(:,1); AoA_data2{j_read} = CD_data_file(:,1);
                        CL_P{j_read} = fit(AoA_data{j_read},CL_data_file(:,2),'linearinterp');
                        CD_P{j_read} = fit(AoA_data2{j_read},CD_data_file(:,2),'linearinterp');
                        
                        alpha_0(j_read) = fzero(CL_P{j_read},0);
                        slope = polyfit(deg2rad(linspace(alpha_0(j_read),alpha_0(j_read)+2,100)), CL_P{j_read}(linspace(alpha_0(j_read),alpha_0(j_read)+2,100)),1);
                        Cla(j_read) = slope(1);
                        j_read = j_read+1;
                    end
                    loghtml(app, 'Loaded existing polar data.', 'info');
                end
                
                if app.StopRequested, error('Analysis stopped by user.'); end
                
                %% ExBEMT Main Loop
                loghtml(app, 'ExBEMT Subroutine is Initiated...', 'header');
                
                if strcmp(analysisType, 'Advance Ratio')
                    outer_loop_vector = alpha_vector;
                    inner_loop_vector = j_vector;
                    colors = lines(length(outer_loop_vector));
                elseif strcmp(analysisType, 'Velocity')
                    outer_loop_vector = alpha_vector;
                    inner_loop_vector = v_vector;
                    colors = lines(length(outer_loop_vector));
                else % Angle of Incidence
                    outer_loop_vector = j_vector;
                    inner_loop_vector = alpha_vector;
                    colors = lines(length(outer_loop_vector));
                end
                
                % Pre-allocate full result matrices
                num_outer = length(outer_loop_vector);
                num_inner = length(inner_loop_vector);
                C_T_all = nan(num_outer, num_inner); C_N_all = nan(num_outer, num_inner); C_Q_all = nan(num_outer, num_inner);
                J_eff_all = nan(num_outer, num_inner);
                                
                % Main calculation loops
                for i_outer = 1:length(outer_loop_vector)
                    if app.StopRequested, error('Analysis stopped by user.'); end
                    iter_tic = tic; % Start iteration timer
                    
                    keep_calculating = true;
                    
                    for i_inner = 1:length(inner_loop_vector)
                        drawnow; % Makes the Stop button responsive
                        if app.StopRequested, error('Analysis stopped by user.'); end
                        
                        if keep_calculating
                            if strcmp(analysisType, 'Angle of Incidence')
                                current_alpha_for_calc = inner_loop_vector(i_inner);
                                current_j_for_calc = outer_loop_vector(i_outer);
                                current_V_inf = current_j_for_calc * (omega*R) / (pi*cosd(current_alpha_for_calc));
                            else 
                                current_alpha_for_calc = outer_loop_vector(i_outer);
                                if strcmp(analysisType, 'Advance Ratio')
                                    current_j_for_calc = inner_loop_vector(i_inner);
                                    current_V_inf = current_j_for_calc * (omega*R) / (pi*cosd(current_alpha_for_calc));
                                else % Velocity
                                    current_V_inf = inner_loop_vector(i_inner);
                                    current_j_for_calc = pi * current_V_inf * cosd(current_alpha_for_calc) / (omega * R);
                                end
                            end
                            
                            V_inf_prll = current_V_inf * sind(current_alpha_for_calc);
                            V_inf_perp = current_V_inf * cosd(current_alpha_for_calc);
                            
                            [T, N_force, Q_torque, ~, ~, ~] = BEMT_inflow_Averaged(V_inf_prll, V_inf_perp, current_alpha_for_calc, beta, Cla, alpha_0, fidelityModel);
                            
                            log_msg = sprintf('Calculating: AOI = <b>%.1f deg</b>, J'' = <b>%.2f</b>', current_alpha_for_calc, current_j_for_calc);
                            loghtml(app, log_msg, 'info');
                            
                            if app.PlotAzimuthalContoursCheckBox.Value
                                plotContourResults(app, fidelityModel, current_alpha_for_calc, current_j_for_calc);
                            end
                                                        
                            J_eff_all(i_outer, i_inner) = current_j_for_calc;
                            C_T_all(i_outer,i_inner) =  T/(plotData.Np*rho*(n)^2*(2*R)^4);
                            C_N_all(i_outer,i_inner) =  N_force/(plotData.Np*rho*(n)^2*(2*R)^4);
                            C_Q_all(i_outer,i_inner) =  Q_torque/(plotData.Np*rho*(n)^2*(2*R)^5);
                            
                            if C_T_all(i_outer,i_inner) <= 1e-2
                                keep_calculating = false;
                            end
                        end
                    end
                    
                    % After completing an inner loop series, plot it
                    plotDataForLine.analysisType = analysisType;
                    plotDataForLine.C_T_line = C_T_all(i_outer,:); 
                    plotDataForLine.C_N_line = C_N_all(i_outer,:);
                    plotDataForLine.C_Q_line = C_Q_all(i_outer,:);
                    plotDataForLine.J_eff_line = J_eff_all(i_outer, :);
                    plotDataForLine.x_vector = inner_loop_vector;
                    
                    plotSingleCurve(app, plotDataForLine, newAxes, outer_loop_vector(i_outer), colors(i_outer,:));
                    iter_time = toc(iter_tic);
                    
                    if strcmp(analysisType, 'Angle of Incidence')
                        log_msg = sprintf('--- Iteration for Series J'' = <b>%.2f</b> complete. Time: <b>%.2f seconds</b> ---', outer_loop_vector(i_outer), iter_time);
                        loghtml(app, log_msg, 'info');
                    else
                        log_msg = sprintf('--- Iteration for Series AOI = <b>%.1f deg</b> complete. Time: <b>%.2f seconds</b> ---', outer_loop_vector(i_outer), iter_time);
                        loghtml(app, log_msg, 'info');
                    end
                    drawnow;
                end
                
                if app.StopRequested, error('Analysis stopped by user.'); end
                                
                %% Finalize plots and log data
                finalizePlots(app, newAxes);
                
                logPlotData(app, 'RESULTS SUMMARY', 'header', [], [], []);
                logPlotData(app, 'Thrust Coefficient (CT)', analysisType, inner_loop_vector, outer_loop_vector, C_T_all);
                logPlotData(app, 'Normal Force Coefficient (CN)', analysisType, inner_loop_vector, outer_loop_vector, C_N_all);
                logPlotData(app, 'Torque Coefficient (CQ)', analysisType, inner_loop_vector, outer_loop_vector, C_Q_all);
                eta_all = (C_T_all .* J_eff_all) ./ (C_Q_all * 2 * pi);
                logPlotData(app, 'Propulsive Efficiency (eta)', analysisType, inner_loop_vector, outer_loop_vector, eta_all);
                
                total_time = toc(total_tic);
                loghtml(app, sprintf('====== TOTAL RUN TIME: %.2f seconds ======', total_time), 'header');
                loghtml(app, 'Analysis finished.', 'success');

            catch E
                loghtml(app, ['<b>SCRIPT ERROR:</b> ' E.message], 'error');
                for i_err = 1:length(E.stack)
                    err_loc = ['<span style="padding-left: 20px;">In ', E.stack(i_err).name, ' at line ', num2str(E.stack(i_err).line), '</span>'];
                    loghtml(app, err_loc, 'error');
                end
            end
            
            % Re-enable run button and disable stop button
            app.RunBEMTAnalysisButton.Enable = 'on';
            app.StopAnalysisButton.Enable = 'off';
            app.ClearPlotsButton.Enable = 'on';
        end

        % Create UIFigure and components
% Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide it until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Icon = 'ExBEMT_Icon.png';
            app.UIFigure.Position = [100 100 1350 950]; % Increased width
            app.UIFigure.Name = 'ExBEMT Analysis GUI V2 (072125)';
            
            % Main App Grid
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {360, '1x'};
            app.GridLayout.RowHeight = {'1x'};

            % Left Panel
            app.LeftPanel = uipanel(app.GridLayout, 'Title', 'Controls');
            leftPanelLayout = uigridlayout(app.LeftPanel);
            leftPanelLayout.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
            leftPanelLayout.ColumnWidth = {'1x'};

            % Right Panel
            app.RightPanel = uipanel(app.GridLayout, 'Title', 'Output');
            rightPanelLayout = uigridlayout(app.RightPanel);
            rightPanelLayout.RowHeight = {'1x', 22, 120};
            rightPanelLayout.ColumnWidth = {'1x'};
            
            % --- Right Panel Components ---
            app.PlotTabGroup = uitabgroup(rightPanelLayout);
            app.PlotTabGroup.Layout.Row = 1;
            app.PlotTabGroup.Layout.Column = 1;

            app.LogTextAreaLabel = uilabel(rightPanelLayout); app.LogTextAreaLabel.Layout.Row = 2;
            app.LogTextAreaLabel.Text = 'Log';
            
            app.LogTextArea = uihtml(rightPanelLayout); 
            app.LogTextArea.Layout.Row = 3;
            html_style = '<style>body { font-family: Segoe UI, Arial, sans-serif; font-size: 10pt; } h2 { color: #0072BD; margin-bottom: 2px; margin-top: 10px;} hr { border: 0; border-top: 1px solid #7E2F8E; margin-top: 2px; margin-bottom: 10px;} table { border-collapse: collapse; width: 95%%; } th, td { border: 1px solid #dddddd; text-align: left; padding: 4px; } th { background-color: #f2f2f2; } </style>';
            app.LogTextArea.HTMLSource = ['<html><head>' html_style '</head><body></body></html>'];
            
            % --- Left Panel Components ---
            
            % General Settings Panel
            app.GeneralSettingsPanel = uipanel(leftPanelLayout, 'Title', 'General Settings');
            app.GeneralSettingsPanel.Layout.Row = 1;
            gsPanelLayout = uigridlayout(app.GeneralSettingsPanel);
            gsPanelLayout.RowHeight = {'fit', 'fit'}; gsPanelLayout.ColumnWidth = {'fit', '1x', 'fit', '1x'};
            uilabel(gsPanelLayout, 'Text', 'Height (m)');
            app.HeightmEditField = uieditfield(gsPanelLayout, 'numeric', 'Value', 0);
            uilabel(gsPanelLayout, 'Text', 'Num. Blades (B)');
            app.NumberofBladesBEditField = uieditfield(gsPanelLayout, 'numeric', 'Value', 2);
            uilabel(gsPanelLayout, 'Text', 'RPM');
            app.RPMEditField = uieditfield(gsPanelLayout, 'numeric', 'Value', 6000);
            app.RPMEditField.ValueChangedFcn = createCallbackFcn(app, @rpmChanged, true);
            
            % Analysis Settings Panel
            app.AnalysisSettingsPanel = uipanel(leftPanelLayout, 'Title', 'Analysis Settings');
            app.AnalysisSettingsPanel.Layout.Row = 2;
            asPanelLayout = uigridlayout(app.AnalysisSettingsPanel);
            asPanelLayout.RowHeight = {'fit', 'fit'}; asPanelLayout.ColumnWidth = {'fit', '1x'};
            uilabel(asPanelLayout, 'Text', 'Propeller');
            app.PropellerDropDown = uidropdown(asPanelLayout); app.PropellerDropDown.Layout.Column = 2;
            uilabel(asPanelLayout, 'Text', 'Analysis Type');
            app.AnalysisTypeDropDown = uidropdown(asPanelLayout, 'Items', {'Select Analysis Type...', 'Advance Ratio', 'Velocity', 'Angle of Incidence'}, 'Value', 'Select Analysis Type...'); 
            app.AnalysisTypeDropDown.Layout.Column = 2;
            app.AnalysisTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @analysisTypeChanged, true);
            
            % Corrections Panel
            app.CorrectionsPanel = uipanel(leftPanelLayout, 'Title', 'Corrections & Models');
            app.CorrectionsPanel.Layout.Row = 3;
            corrPanelLayout = uigridlayout(app.CorrectionsPanel);
            corrPanelLayout.RowHeight = {'fit', 'fit'};
            corrPanelLayout.ColumnWidth = {'1x', '1x'};

            app.RadialFlowCheckBox = uicheckbox(corrPanelLayout, 'Text', 'Radial Flow', 'Value', true);
            app.RadialFlowCheckBox.Layout.Row = 1;
            app.RadialFlowCheckBox.Layout.Column = 1;

            app.Use3DPolarCheckBox = uicheckbox(corrPanelLayout, 'Text', '3D Polar', 'Value', true);
            app.Use3DPolarCheckBox.Layout.Row = 1;
            app.Use3DPolarCheckBox.Layout.Column = 2;
            
            app.FidelityModelDropDown = uidropdown(corrPanelLayout, 'Items', {'BEMT', 'BEMT+ Dynamic Inflow Model', 'BEMT+ Q-Unsteady Dynamic Inflow Model'}, 'Value', 'BEMT');
            app.FidelityModelDropDown.Layout.Row = 2; 
            app.FidelityModelDropDown.Layout.Column = [1 2];
            
            % Run Options Panel
            app.RunOptionsPanel = uipanel(leftPanelLayout, 'Title', 'Run Options');
            app.RunOptionsPanel.Layout.Row = 4;
            roPanelLayout = uigridlayout(app.RunOptionsPanel);
            roPanelLayout.RowHeight = {'fit'};
            roPanelLayout.ColumnWidth = {'1x', '1x'};
            app.UseOldPolarDataCheckBox = uicheckbox(roPanelLayout, 'Text', 'Use Existing Polars', 'Value', true);
            app.UseOldPolarDataCheckBox.Layout.Column = 1;
            app.UseOldPolarDataCheckBox.ValueChangedFcn = createCallbackFcn(app, @oldDataCheckChanged, true);
            app.PlotAzimuthalContoursCheckBox = uicheckbox(roPanelLayout, 'Text', 'Plot Contour Plots', 'Value', false);
            app.PlotAzimuthalContoursCheckBox.Layout.Column = 2;
            
            % Analysis Range Panel
            app.AnalysisRangePanel = uipanel(leftPanelLayout, 'Title', 'Analysis Range');
            app.AnalysisRangePanel.Layout.Row = 5;
            app.AnalysisRangePanel.Visible = 'off';
            arPanelLayout = uigridlayout(app.AnalysisRangePanel);
            arPanelLayout.RowHeight = {'fit', 'fit', 'fit'}; arPanelLayout.ColumnWidth = {'1x'};
            app.AlphaRangePanel = uipanel(arPanelLayout);
            app.JPrimeRangePanel = uipanel(arPanelLayout);
            app.VelocityRangePanel = uipanel(arPanelLayout);
            
            app.AlphaRangePanel.Title = 'Angle of Incidence (deg)';
            alphaLayout = uigridlayout(app.AlphaRangePanel);
            alphaLayout.RowHeight = {'fit'}; alphaLayout.ColumnWidth = {'fit', '1x', 'fit', '1x', 'fit', '1x'};
            uilabel(alphaLayout, 'Text', 'Min');
            app.MinAlphaEditField = uieditfield(alphaLayout, 'numeric', 'Value', 0);
            uilabel(alphaLayout, 'Text', 'Step');
            app.AlphaStepEditField = uieditfield(alphaLayout, 'numeric', 'Value', 5);
            uilabel(alphaLayout, 'Text', 'Max');
            app.MaxAlphaEditField = uieditfield(alphaLayout, 'numeric', 'Value', 40);
            
            app.JPrimeRangePanel.Title = 'Advance Ratio (J'') Range';
            jLayout = uigridlayout(app.JPrimeRangePanel);
            jLayout.RowHeight = {'fit'}; jLayout.ColumnWidth = {'fit', '1x', 'fit', '1x', 'fit', '1x'};
            uilabel(jLayout, 'Text', 'Min');
            app.MinJEditField = uieditfield(jLayout,'numeric','Value',0.1);
            uilabel(jLayout, 'Text', 'Step');
            app.JStepEditField = uieditfield(jLayout,'numeric','Value',0.05);
            uilabel(jLayout, 'Text', 'Max');
            app.MaxJEditField = uieditfield(jLayout,'numeric','Value',0.7);

            app.VelocityRangePanel.Title = 'Velocity (m/s) Range';
            vLayout = uigridlayout(app.VelocityRangePanel);
            vLayout.RowHeight = {'fit'}; vLayout.ColumnWidth = {'fit', '1x', 'fit', '1x', 'fit', '1x'};
            uilabel(vLayout, 'Text', 'Min');
            app.MinVEditField = uieditfield(vLayout,'numeric','Value',2);
            uilabel(vLayout, 'Text', 'Step');
            app.VStepEditField = uieditfield(vLayout,'numeric','Value',2);
            uilabel(vLayout, 'Text', 'Max');
            app.MaxVEditField = uieditfield(vLayout,'numeric','Value',10);
            
            % XFOIL Settings Panel (MODIFIED)
            app.XFOILSettingsPanel = uipanel(leftPanelLayout, 'Title', 'XFOIL Settings');
            xfoilLayout = uigridlayout(app.XFOILSettingsPanel);
            xfoilLayout.RowHeight = {'fit', 'fit', 'fit'}; % 3 rows
            xfoilLayout.ColumnWidth = {'fit', '1x', 'fit', '1x'};
            
            uilabel(xfoilLayout, 'Text', 'Upper AoA');
            app.UpperAoAEditField = uieditfield(xfoilLayout, 'numeric', 'Value', 20);
            uilabel(xfoilLayout, 'Text', 'Lower AoA');
            app.LowerAoAEditField = uieditfield(xfoilLayout, 'numeric', 'Value', -15);
            
            uilabel(xfoilLayout, 'Text', 'NCrit');
            app.NCritEditField = uieditfield(xfoilLayout, 'numeric', 'Value', 3);
            uilabel(xfoilLayout, 'Text', 'Filter Passes');
            app.FilterPassesEditField = uieditfield(xfoilLayout, 'numeric', 'Value', 10);

            uilabel(xfoilLayout, 'Text', 'Max Iterations');
            app.MaxIterationsEditField = uieditfield(xfoilLayout, 'numeric', 'Value', 200);
            app.PlotXFOILPolarsCheckBox = uicheckbox(xfoilLayout, 'Text', 'Plot XFOIL Polars');
            app.PlotXFOILPolarsCheckBox.Layout.Column = [3 4]; % Span last two columns

            % Button Panel
            buttonPanel = uipanel(leftPanelLayout, 'BorderType', 'none');
            buttonPanel.Layout.Row = 7;
            buttonLayout = uigridlayout(buttonPanel);
            buttonLayout.ColumnWidth = {'1.5x', '1x', '1.2x', '1.2x'}; 
            buttonLayout.RowHeight = {'fit'};
            
            app.RunBEMTAnalysisButton = uibutton(buttonLayout, 'push', 'Text', 'Run Analysis');
            app.RunBEMTAnalysisButton.BackgroundColor = [0.07, 0.62, 1.00];
            app.RunBEMTAnalysisButton.FontColor = [1 1 1];
            app.RunBEMTAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @runBEMTAnalysisButtonPushed, true);

            app.StopAnalysisButton = uibutton(buttonLayout, 'push', 'Text', 'Stop', 'Enable', 'off');
            app.StopAnalysisButton.BackgroundColor = [0.85, 0.33, 0.10];
            app.StopAnalysisButton.FontColor = [1 1 1];
            app.StopAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @stopAnalysisButtonPushed, true);

            app.ClearPlotsButton = uibutton(buttonLayout, 'push', 'Text', 'Clear Plots');
            app.ClearPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @clearPlotsButtonPushed, true);
            
            app.CloseAllPlotsButton = uibutton(buttonLayout, 'push', 'Text', 'Close All');
            app.CloseAllPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @closeAllPlotsButtonPushed, true);
            
            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
        
        % UPDATED HTML Logging function with JavaScript scrolling
        function loghtml(app, message, style)
            
            % Define styles
            switch style
                case 'error'
                    html_line = sprintf('<p style="color:red;">%s</p>', message);
                case 'warning'
                    html_line = sprintf('<p style="color:#D95319;"><b>Warning:</b> %s</p>', message);
                case 'success'
                    html_line = sprintf('<p style="color:green;">%s</p>', message);
                case 'header'
                    html_line = sprintf('<h2>%s</h2><hr>', message);
                case 'table'
                    html_line = ['<pre>' message '</pre>']; % Use preformatted text for tables
                otherwise % 'info'
                    html_line = sprintf('<p>%s</p>', message);
            end

            % JavaScript to scroll to the bottom
            js_scroll = '<script>window.scrollTo(0, document.body.scrollHeight);</script>';
            
            % Append new line and scroll script to the existing HTML body
            current_html = app.LogTextArea.HTMLSource;
            new_html = strrep(current_html, '</body>', [html_line '</body>']);
            new_html = strrep(new_html, '</html>', [js_scroll '</html>']);
            
            app.LogTextArea.HTMLSource = new_html;
            drawnow('limitrate');
        end
        
        % Populate the propeller dropdown
        function populatePropellers(app)
            try
                current_dir = pwd;
                propPackagesDir = fullfile(current_dir, 'Propeller Packages');
                
                loghtml(app, ['Searching for propellers in: ', propPackagesDir], 'info');

                if ~exist(propPackagesDir, 'dir')
                    loghtml(app, '"Propeller Packages" folder not found in current directory.', 'error');
                    loghtml(app, 'Please navigate to your project folder in MATLAB before running the app.', 'warning');
                    app.PropellerDropDown.Items = {'Folder not found'};
                    return;
                end
                
                files = dir(propPackagesDir);
                propNames_temp = files([files.isdir] & ~ismember({files.name}, {'.','..'}));
                
                if isempty(propNames_temp)
                    loghtml(app, 'No sub-folders found in "Propeller Packages".', 'warning');
                    app.PropellerDropDown.Items = {'None found'};
                else
                    propNames = cell(1, numel(propNames_temp));
                    for k = 1:numel(propNames_temp)
                        propNames{k} = propNames_temp(k).name;
                    end
                    app.PropellerDropDown.Items = [{'Select Propeller...'}, propNames];
                    app.PropellerDropDown.Value = 'Select Propeller...';
                end
            catch ME
                loghtml(app, ['Error populating propellers: ', ME.message], 'error');
                app.PropellerDropDown.Items = {'Error'};
            end
        end
        
        % SIMPLIFIED Plotting helper function
        function plotSingleCurve(app, plotData, axesHandles, series_val, color)
            
            if strcmp(plotData.analysisType, 'Velocity')
                x_vector = plotData.x_vector;
                x_label_text = 'Velocity (m/s)';
                legend_prefix = '\alpha = ';
                legend_suffix = ' deg';
            elseif strcmp(plotData.analysisType, 'Advance Ratio')
                x_vector = plotData.x_vector;
                x_label_text = 'Effective Advance Ratio (J'')';
                legend_prefix = '\alpha = ';
                legend_suffix = ' deg';
            else % Angle of Incidence
                x_vector = plotData.x_vector;
                x_label_text = 'Angle of Incidence (deg)';
                legend_prefix = 'J'' = ';
                legend_suffix = '';
            end
            
            linewidth = 1.5;
            markerstyle = 'o';
            markersize = 4;
            linestyle = '-';

            lgd_text = sprintf('%s%.1f%s', legend_prefix, series_val, legend_suffix);

            % Plot CT
            plot(axesHandles.CT, x_vector, plotData.C_T_line, 'LineStyle', linestyle, 'Marker', markerstyle, 'MarkerSize', markersize, 'Color', color, 'LineWidth', linewidth, 'DisplayName', lgd_text);
            xlabel(axesHandles.CT, x_label_text); ylabel(axesHandles.CT, 'C_T');

            % Plot CN
            plot(axesHandles.CN, x_vector, plotData.C_N_line, 'LineStyle', linestyle, 'Marker', markerstyle, 'MarkerSize', markersize, 'Color', color, 'LineWidth', linewidth, 'HandleVisibility', 'off');
            xlabel(axesHandles.CN, x_label_text); ylabel(axesHandles.CN, 'C_N');
            
            % Plot CQ
            plot(axesHandles.CQ, x_vector, plotData.C_Q_line, 'LineStyle', linestyle, 'Marker', markerstyle, 'MarkerSize', markersize, 'Color', color, 'LineWidth', linewidth, 'HandleVisibility', 'off');
            xlabel(axesHandles.CQ, x_label_text); ylabel(axesHandles.CQ, 'C_Q');
            
            % Plot Efficiency
            eta = (plotData.C_T_line .* plotData.J_eff_line) ./ (plotData.C_Q_line * 2 * pi);
            eta(eta<0 | eta > 1.2) = NaN; % Clean up efficiency values
            plot(axesHandles.Efficiency, x_vector, eta, 'LineStyle', linestyle, 'Marker', markerstyle, 'MarkerSize', markersize, 'Color', color, 'LineWidth', linewidth, 'HandleVisibility', 'off');
            xlabel(axesHandles.Efficiency, x_label_text); ylabel(axesHandles.Efficiency, '\eta');
        end

        % CORRECTED Data logging helper function
        function logPlotData(app, titleStr, analysisType, x_vector, series_vector, y_matrix)
            try
                if isempty(x_vector) || isempty(series_vector) || isempty(y_matrix)
                    loghtml(app, titleStr, 'header');
                    return;
                end

                if strcmp(analysisType, 'Angle of Incidence')
                    x_label = 'AOI(deg)';
                    series_label_prefix = 'J''=';
                else
                    x_label = 'J''/V';
                    series_label_prefix = 'AOI=';
                end
                
                if ~ismatrix(y_matrix)
                    error('Input data is not a 2D matrix.');
                end
                y_matrix_log = y_matrix';

                header_format = ['%-10s' repmat('%-10s', 1, length(series_vector)) '\n'];
                % header_format = strjoin(header_format, '');
                row_format = ['%-10.2f' repmat('%-10.4f', 1, length(series_vector)) '\n'];
                
                header_labels = {x_label};
                for i = 1:length(series_vector)
                    header_labels{end+1} = [series_label_prefix num2str(series_vector(i))];
                end
                
                table_str = sprintf(header_format, header_labels{:});
                table_str = [table_str, repmat('-', 1, 10*(length(series_vector)+1)), '\n'];
                
                for i = 1:length(x_vector)
                    row_data = [x_vector(i), y_matrix_log(i,:)];
                    table_str = [table_str, sprintf(row_format, row_data)];
                end
                
                loghtml(app, ['<h3>' titleStr '</h3>'], 'info');
                loghtml(app, table_str, 'table');

            catch ME_log
                loghtml(app, ['<b>Could not generate results table for "' titleStr '".</b> Error: ' ME_log.message], 'error');
            end
        end
        
        % SIMPLIFIED Finalize Plots function
        function finalizePlots(app, axesHandles)
            all_axes = [axesHandles.CT, axesHandles.CN, axesHandles.CQ, axesHandles.Efficiency];
            for ax = all_axes
                hold(ax, 'off');
                grid(ax, 'on');
                box(ax, 'on');
            end
            
            lgd = legend(axesHandles.CT, 'Location', 'northeastoutside');
            title(lgd, 'Analysis Series')
        end
        
        % UPDATED: Contour plotting helper function
        function plotContourResults(app, fidelityModel, alpha, J_prime)
            loghtml(app, 'Generating contour plots...', 'info');
            
            try
                % Fetch data from base workspace
                x_blade = evalin('base', 'x_blade');
                B = evalin('base', 'B');

                % Define plot layout and variables based on fidelity model
                if strcmp(fidelityModel, 'BEMT')
                    figure_title = 'BEMT Contour Plots';
                    var_list = {'dT_BEMT', 'dQ_BEMT', 'Vel_BEMT', 'AoA_BEMT'};
                    subplot_dims = [1, 4];
                else
                    figure_title = 'BEMT vs. Dynamic Inflow Contour Plots';
                    var_list = {'dT_BEMT', 'dQ_BEMT', 'Vel_BEMT', 'AoA_BEMT', 'dT_PITT', 'dQ_PITT', 'Vel_PITT', 'AoA_PITT'};
                    subplot_dims = [2, 4];
                end

                fig = figure('Name', figure_title, 'NumberTitle', 'off', 'Position', [50, 50, 1600, 800]);
                sgtitle(fig, sprintf('Contour Plots for AOI = %.1f deg, J'' = %.2f', alpha, J_prime), 'FontSize', 14, 'FontWeight', 'bold');
                
                num_sections = length(x_blade);

                for i = 1:length(var_list)
                    var_name = var_list{i};
                    
                    if ~evalin('base', ['exist(''' var_name ''', ''var'')'])
                        loghtml(app, ['Contour plot skipped: Variable ' var_name ' not found.'], 'warning');
                        ax = subplot(subplot_dims(1), subplot_dims(2), i);
                        text(ax, 0.5, 0.5, [strrep(var_name, '_', '\_'), ' not found'], 'HorizontalAlignment', 'center');
                        axis(ax, 'off');
                        continue;
                    end
                    
                    prop_plot_cell = evalin('base', var_name);
                    
                    num_azi_points_per_blade = size(prop_plot_cell{1,1}, 2);
                    
                    prop_plot_matrix = zeros(num_sections, (num_azi_points_per_blade-1) * B);
                    
                    for sec_idx = 1:num_sections
                        full_rotation_data = [];
                        for blade_idx = 1:B
                            if size(prop_plot_cell, 1) >= sec_idx && size(prop_plot_cell, 2) >= blade_idx && ~isempty(prop_plot_cell{sec_idx, blade_idx})
                                blade_data = prop_plot_cell{sec_idx, blade_idx};
                                blade_data = blade_data(1:end-1);
                                if size(blade_data,1) > 1, blade_data = blade_data'; end
                                full_rotation_data = [full_rotation_data, blade_data];
                            end
                        end
                        prop_plot_matrix(sec_idx, 1:length(full_rotation_data)) = full_rotation_data;
                    end
                    
                    [Psi, R_grid] = meshgrid(linspace(0, 2*pi, size(prop_plot_matrix, 2)), x_blade);

                    [X, Y] = pol2cart(Psi, R_grid);
                    X_rotated = Y;  
                    Y_rotated = -X;

                    ax = subplot(subplot_dims(1), subplot_dims(2), i);
                    contourf(ax, X_rotated, Y_rotated, prop_plot_matrix, 20, 'LineStyle', 'none');
                    colorbar(ax);
                    
                    plot_title = '';
                    if contains(var_name, 'dT'), plot_title = 'Thrust (dT)';
                    elseif contains(var_name, 'dQ'), plot_title = 'Torque (dQ)';
                    elseif contains(var_name, 'Vel'), plot_title = 'Relative Velocity (m/s)';
                    elseif contains(var_name, 'AoA'), plot_title = 'Angle of Attack (deg)';
                    else, plot_title = strrep(var_name, '_', '\_');
                    end

                    title(ax, plot_title);
                    if contains(var_name, 'PITT'), title(ax, plot_title, 'Color', 'r'); end
                    
                    axis(ax, 'equal');
                    axis(ax, 'off');
                end
                
                % Add row labels for the 2x4 plot using annotation
                if subplot_dims(1) == 2 && ishandle(fig)
                    try
                        ax1 = subplot(2,4,1);
                        pos1 = get(ax1, 'Position');
                        annotation(fig, 'textbox', [pos1(1)-0.07, pos1(2), 0.05, pos1(4)], ...
                            'String', 'BEMT', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');

                        ax5 = subplot(2,4,5);
                        pos5 = get(ax5, 'Position');
                        annotation(fig, 'textbox', [pos5(1)-0.07, pos5(2), 0.05, pos5(4)], ...
                            'String', {'Dynamic', 'Inflow'}, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
                    catch ME_annot
                        loghtml(app, 'Could not create row labels for contour plot.', 'warning');
                    end
                end
                
                drawnow; % Force the graphics queue to be processed
            catch ME
                loghtml(app, ['<b>Contour Plotting Error:</b> ' ME.message], 'error');
                for i_err = 1:length(ME.stack)
                    err_loc = ['<span style="padding-left: 20px;">In ', ME.stack(i_err).name, ' at line ', num2str(E.stack(i_err).line), '</span>'];
                    loghtml(app, err_loc, 'error');
                end
            end
        end
    end
end
