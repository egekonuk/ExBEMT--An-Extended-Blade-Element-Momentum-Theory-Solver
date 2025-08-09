function [TT_final, TN_final, TQ_final, plotData] = BEMT_inflow_Averaged(V_inf_prll, V_inf_perp, settings, CL_P, CD_P)
%BEMT_inflow_Averaged Performs BEMT analysis with optional dynamic inflow models.
% This function calculates propeller performance over a full rotation for
% each blade section, averages the results, and can apply a second pass
% using a Pitt-Peters dynamic inflow model (steady or unsteady) to correct
% the induced velocity.
%
% VERSION 3 

%% 1. INITIALIZATION AND SETUP
% =========================================================================

% --- Extract Key Parameters from Settings Struct ---
% This improves readability by using local variables within the function.
B = settings.B; R = settings.R; omega = settings.omega; rho = settings.rho;
Np = settings.Np; delta_x = settings.delta_x; c = settings.c;
x_blade = settings.x_blade; radial_correct = settings.radial_correct;
Use_3D_polar = settings.Use_3D_polar; fidelityModel = settings.fidelityModel;
beta = settings.beta;

% --- General Initializations ---
sep = 2 * pi / B;           % Azimuthal separation between blades in radians.
Ua = V_inf_perp;            % Axial free-stream velocity (perpendicular to rotor disk).
Uy = V_inf_prll;            % In-plane free-stream velocity (parallel to rotor disk).
dpsi_deg = settings.dpsi_deg; % Azimuthal step size for the integration, in degrees.
dpsi_rad = deg2rad(dpsi_deg); % Convert step size to radians.
delta_time = dpsi_rad / omega;  % Time step corresponding to the azimuthal step.

% --- Determine which Calculation Path to Use ---
run_di_calcs = ~strcmp(fidelityModel, 'BEMT'); % Flag to run dynamic inflow (DI) calculations.
is_unsteady = strcmp(fidelityModel, 'BEMT+ Q-Unsteady Dynamic Inflow Model'); % Flag for unsteady DI model.
is_steady = strcmp(fidelityModel, 'BEMT+ Dynamic Inflow Model'); % Flag for steady DI model.

% --- Pre-allocate Cell Arrays for Plotting Data ---
% Pre-allocation is crucial for performance, especially inside loops.
% Cells are used because the number of azimuthal steps can vary.
num_sections = length(x_blade);
% --- BEMT Pass Results ---
dT_BEMT_plot = cell(num_sections, B); dQ_BEMT_plot = cell(num_sections, B);
Vel_BEMT_plot = cell(num_sections, B); AoA_BEMT_plot = cell(num_sections, B);
phi_BEMT_plot = cell(num_sections, B); In_axial_BEMT_plot = cell(num_sections, B);
In_swirl_BEMT_plot = cell(num_sections, B); psi_BEMT_plot = cell(num_sections, B);

% --- Dynamic Inflow (Pitt-Peters) Pass Results ---
dT_PITT_plot = cell(num_sections, B); dQ_PITT_plot = cell(num_sections, B);
Vel_PITT_plot = cell(num_sections, B); AoA_PITT_plot = cell(num_sections, B);
phi_PITT_plot = cell(num_sections, B); In_axial_PITT_plot = cell(num_sections, B);

% --- Pre-allocate Blade-Level Totals ---
TT_blade = zeros(1, B); TN_blade = zeros(1, B); TQ_blade = zeros(1, B);

%% 2. MAIN BEMT LOOP (PER BLADE)
% =========================================================================
% This loop iterates through each blade of the propeller.
for bn = 1:B
    % --- Pre-allocate section-level averages for the current blade ---
    T_section_avg = zeros(num_sections, 1); % Average Thrust
    N_section_avg = zeros(num_sections, 1); % Average Normal Force
    Q_section_avg = zeros(num_sections, 1); % Average Torque

    %% 2A. INNER LOOP (PER BLADE SECTION)
    % This loop iterates through each radial station along the current blade.
    for i = 1:num_sections
        % --- Calculate local solidity for the current section ---
        sigmal = B*c(i)/(2*pi*(x_blade(i)*R));
        num_sector_steps = round((2*pi/B/omega) / delta_time) + 1;
        
        % --- Pre-allocate arrays for the azimuthal sweep ---
        psi = zeros(1, num_sector_steps);   % Azimuth angle
        dT = zeros(1, num_sector_steps);    % Thrust differential
        dQ = zeros(1, num_sector_steps);    % Torque differential
        dN = zeros(1, num_sector_steps);    % Normal force differential
        dM = zeros(1, num_sector_steps);    % Pitching moment differential
        dn = zeros(1, num_sector_steps);    % Yawing moment differential
        d2M = zeros(1, num_sector_steps);   % Second harmonic pitching moment
        d2n = zeros(1, num_sector_steps);   % Second harmonic yawing moment
        w_a = zeros(1, num_sector_steps);   % Induced axial velocity
        
        vel_rev = zeros(1, num_sector_steps); % Relative velocity
        aoa_rev = zeros(1, num_sector_steps); % Angle of attack
        phi_rev = zeros(1, num_sector_steps); % Inflow angle
        in_axial_rev = zeros(1, num_sector_steps); % Axial inflow component
        in_swirl_rev = zeros(1, num_sector_steps); % Swirl inflow component
        
        j = 0; % Azimuthal step counter

        %% 2B. AZIMUTHAL SWEEP (BEMT PASS)
        % This loop integrates forces over one loading sector of the rotor
        % disk. aka 180 deg for 2 blades 120 deg for 3 blades etc.
        for t = 0:delta_time:(2*pi/B/omega)
            j = j + 1;
            psi(j) = sep*(bn-1) + t*omega; % Current azimuth angle
            
            % --- Calculate local velocity components at the blade section ---
            UT_j = Uy*sin(psi(j)); % Tangential component from free-stream
            if radial_correct, UR_j = Uy*cos(psi(j)); else, UR_j = 0; end % Radial component
            Ut_pr_j = omega*x_blade(i)*R + UT_j; % Total tangential velocity
            
            % --- CORE ALGORITHM: Solve for inflow angle and induced velocity using Momentum Theory ---
            [H, VA_W, VT_W, CL_final, CD_final, phi_j, b_alpha_j, ~] = phi_algorithm_ISAE(Ua, UT_j, UR_j, Ut_pr_j, beta(i), i, sigmal, settings, CL_P, CD_P);
            
            % --- Calculate resultant velocities and forces ---
            V_rel_j = Ut_pr_j/H;
            V_rel_y_j = sqrt(V_rel_j^2 + UR_j^2);
            lambda_y_j = atand(UR_j/(V_rel_j*cosd(phi_j))); % Correction for radial flow
            
            % --- Store BEMT results for this azimuthal step ---
            w_a(j) = V_rel_j * VA_W; % Axial induced velocity from Momentum Conservation
            vel_rev(j) = V_rel_y_j;
            aoa_rev(j) = b_alpha_j;
            phi_rev(j) = phi_j;
            in_axial_rev(j) = w_a(j);
            in_swirl_rev(j) = V_rel_j * VT_W;
            
            % --- Calculate Lift and Drag differentials ---
            dL = 0.5*rho*V_rel_j^2*c(i)*CL_final;
            dD = 0.5*rho*V_rel_y_j^2*c(i)*CD_final;
            
            % --- Decompose forces into Thrust, Torque, etc. ---
            dT(j) = (dL*cosd(phi_j) - dD*sind(phi_j))*delta_x(i);
            dM(j) = -dT(j)*x_blade(i)*R*cos(psi(j));
            dn(j) = -dT(j)*x_blade(i)*R*sin(psi(j));
            d2M(j) = -dT(j)*(x_blade(i)*R)^2*cos(2*psi(j));
            d2n(j) = -dT(j)*(x_blade(i)*R)^2*sin(2*psi(j));
            dQ_r_j = ((dL*sind(phi_j) + dD*cosd(phi_j))*cosd(lambda_y_j)*delta_x(i));
            dN(j) = dQ_r_j*sin(psi(j));
            dQ(j) = dQ_r_j*x_blade(i)*R;
        end % End of azimuthal sweep
        
        % --- Store full azimuthal data for plotting (BEMT Pass) ---
        psi_BEMT_plot{i,bn} = rad2deg(psi);
        dT_BEMT_plot{i,bn} = dT; dQ_BEMT_plot{i,bn} = dQ;
        Vel_BEMT_plot{i,bn} = vel_rev; AoA_BEMT_plot{i,bn} = aoa_rev;
        phi_BEMT_plot{i,bn} = phi_rev; In_axial_BEMT_plot{i,bn} = in_axial_rev;
        In_swirl_BEMT_plot{i,bn} = in_swirl_rev;

        % --- Calculate the average forces for this section from the BEMT pass ---
        T_section_avg(i) = mean(dT);
        N_section_avg(i) = mean(dN);
        Q_section_avg(i) = mean(dQ);
        
        %% 3. DYNAMIC INFLOW (PITT-PETERS) CORRECTION
        % =========================================================================
        % This section runs only if a dynamic inflow model was selected.
        % It uses the BEMT results to calculate a more accurate induced velocity.
        if run_di_calcs
            % --- Calculate average moments from the BEMT pass ---
            M_avg = mean(dM);
            n_avg = mean(dn);
            M2_avg = mean(d2M);
            n2_avg = mean(d2n);
            wa_avg_sec = mean(w_a); % Average induced velocity from BEMT
            
            % --- Calculate STEADY Pitt-Peters coefficients (if applicable) ---
            if is_steady
                [~, nu_coeff] = Peters_dynamic_inflow(T_section_avg(i), n_avg, M_avg, n2_avg, M2_avg, Ua, Uy, wa_avg_sec, omega, R, rho);
                nu_0 = nu_coeff(1); nu_s = nu_coeff(2); nu_c = nu_coeff(3); nu_2s = nu_coeff(4); nu_2c = nu_coeff(5);
            end
            
            % --- Pre-allocate arrays for the DI pass ---
            dT_P = zeros(1,j); dQ_P = zeros(1,j); dN_P = zeros(1,j);
            w_a_final = zeros(1,j); % Final corrected induced velocity
            
            %% 3A. AZIMUTHAL SWEEP (DYNAMIC INFLOW PASS)
            % This loop re-calculates forces using the corrected induced velocity.
            for j_di = 1:j
                % --- Calculate UNSTEADY Pitt-Peters coefficients (if applicable) ---
                if is_unsteady
                    % --- Setup finite difference stencil for time derivatives ---
                    if j_di == 1, diff_idx=1; stencil_range = j_di:j_di+1;       % Forward difference
                    elseif j_di == j, diff_idx=2; stencil_range = j_di-1:j_di;  % Backward difference
                    else, diff_idx=2; stencil_range = j_di-1:j_di+1;          % Central difference
                    end
                    dTs=dT(stencil_range); dns=dn(stencil_range); dMs=dM(stencil_range); 
                    d2ns=d2n(stencil_range); d2Ms=d2M(stencil_range); w_as=w_a(stencil_range);
                    
                    [~, nu_coeff_unsteady] = Peters_dynamic_unsteady_inflow(dTs,dns,dMs,d2ns,d2Ms,Ua,Uy,w_as,diff_idx,delta_time, omega, R, rho);
                    nu_0 = nu_coeff_unsteady(1); nu_s = nu_coeff_unsteady(2); nu_c = nu_coeff_unsteady(3); 
                    nu_2s = nu_coeff_unsteady(4); nu_2c = nu_coeff_unsteady(5);
                end
                
                % --- Calculate the final, corrected induced velocity at this azimuth ---
                psi_current = psi(j_di);
                w_a_final(j_di) = (nu_0 + nu_s*x_blade(i)*sin(psi_current) + nu_c*x_blade(i)*cos(psi_current) ...
                    + nu_2s*x_blade(i)^2*sin(2*psi_current) + nu_2c*x_blade(i)^2*cos(2*psi_current))*omega*R;
                
                % --- Re-calculate local flow conditions with the new induced velocity ---
                Ua_eff = Ua + w_a_final(j_di);
                Ut_pr_j = omega*x_blade(i)*R + Uy*sin(psi_current);
                UR_j = Uy*cos(psi_current);
                
                % --- Re-run the core algorithm with corrected inflow ---
                [H, ~, ~, CL_final, CD_final, phi_j, b_alpha_j, ~] = phi_algorithm_ISAE(Ua_eff, Ut_pr_j, UR_j, Ut_pr_j, beta(i), i, sigmal, settings, CL_P, CD_P);
                V_rel_j = Ut_pr_j/H;
                V_rel_y_j = sqrt(V_rel_j^2 + UR_j^2);
                
                % --- Re-calculate Lift and Drag with corrected flow ---
                dL_P = 0.5*rho*V_rel_j^2*c(i)*CL_final;
                dD_P = 0.5*rho*V_rel_y_j^2*c(i)*CD_final;
                
                % --- Decompose and store final forces ---
                dT_P(j_di) = (dL_P*cosd(phi_j) - dD_P*sind(phi_j))*delta_x(i);
                dQ_r_P = ((dL_P*sind(phi_j) + dD_P*cosd(phi_j))*delta_x(i))*cosd(atand(UR_j/(V_rel_j*cosd(phi_j))));
                dN_P(j_di) = dQ_r_P*sin(psi_current);
                dQ_P(j_di) = dQ_r_P*x_blade(i)*R;

                % --- Store full azimuthal data for plotting (DI Pass) ---
                Vel_PITT_plot{i,bn}(j_di) = V_rel_y_j;
                AoA_PITT_plot{i,bn}(j_di) = b_alpha_j;
                phi_PITT_plot{i,bn}(j_di) = phi_j;
                In_axial_PITT_plot{i,bn}(j_di) = w_a_final(j_di);
            end % End of DI azimuthal sweep
            
            % --- Store final corrected thrust and torque for plotting ---
            dT_PITT_plot{i,bn} = dT_P;
            dQ_PITT_plot{i,bn} = dQ_P;
            
            % --- OVERWRITE the average forces with the corrected DI results ---
            T_section_avg(i) = mean(dT_P);
            N_section_avg(i) = mean(dN_P);
            Q_section_avg(i) = mean(dQ_P);
        end % End of if(run_di_calcs)
    end % End of section loop
    
    % --- Sum the averaged section forces to get the total for this blade ---
    TT_blade(bn) = sum(T_section_avg);
    TN_blade(bn) = sum(N_section_avg);
    TQ_blade(bn) = sum(Q_section_avg);
end % End of blade loop

%% 4. FINAL AGGREGATION AND OUTPUT
% =========================================================================

% --- Sum the forces from all blades for the final propeller totals ---
TT_final = sum(TT_blade);
TN_final = sum(TN_blade);
TQ_final = sum(TQ_blade);

% --- Package all plotting variables into a struct to be returned ---
plotData = struct();
% --- BEMT Pass Data ---
plotData.dT_BEMT = dT_BEMT_plot;
plotData.dQ_BEMT = dQ_BEMT_plot;
plotData.Vel_BEMT = Vel_BEMT_plot;
plotData.AoA_BEMT = AoA_BEMT_plot;
plotData.phi_BEMT = phi_BEMT_plot;
plotData.In_axial_BEMT = In_axial_BEMT_plot;
plotData.In_swirl_BEMT = In_swirl_BEMT_plot;
plotData.psi_BEMT = psi_BEMT_plot;

% --- Dynamic Inflow Pass Data (if run) ---
if run_di_calcs
    plotData.dT_PITT = dT_PITT_plot;
    plotData.dQ_PITT = dQ_PITT_plot;
    plotData.Vel_PITT = Vel_PITT_plot;
    plotData.AoA_PITT = AoA_PITT_plot;
    plotData.phi_PITT = phi_PITT_plot;
    plotData.In_axial_PITT = In_axial_PITT_plot;
end

end % End of function