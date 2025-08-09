function [H,VA_W,VT_W,CL_final,CD_final,phi,b_alpha,b_alpha_y] = phi_algorithm_ISAE(Ua,UT,UR,Ut_pr,beta_i,i,sigmal,settings,CL_P,CD_P)
%phi_algorithm_ISAE Solves for the inflow angle and resulting aerodynamic
%coefficients for a single blade element using an iterative approach.
% VERSION 2.1 - Uses GUI settings for NF ranges and improved comments.

%% EXTRACT DATA FROM SETTINGS STRUCT
omega = settings.omega; R = settings.R; c = settings.c;
x_blade = settings.x_blade; Use_3D_polar = settings.Use_3D_polar;
use_neuralfoil = settings.use_neuralfoil;

%% DETERMINE LIFT-CURVE SLOPE (Cla) and ZERO-LIFT AOA (alpha_0)
% This is required for the 3D corrections. The method depends on the polar source.
if use_neuralfoil
    % --- On-the-fly calculation for NeuralFoil ---
    % For NeuralFoil, we can determine Cla and alpha_0 dynamically from the
    % pre-computed interpolant for the current flight condition.
    V_rel_approx = sqrt(Ua^2 + Ut_pr^2);
    reynolds_number = (settings.rho * V_rel_approx * settings.c(i)) / settings.nu;
    
    % Use the min/max Re from the GUI settings for clamping
    re_min = settings.nf_re_min; 
    re_max = settings.nf_re_max;
    reynolds_number_clamped = max(re_min, min(re_max, reynolds_number));
    
    % Use the AoA domain from the GUI settings to find alpha_0 and Cla
    aoa_domain_row = -20:0.5:20;
    aoa_domain_col = aoa_domain_row'; 
    
    % Query the interpolant to get the CL curve at this Reynolds number
    re_query_array = ones(size(aoa_domain_col)) * reynolds_number_clamped;
    cl_at_re = settings.CL_NF_Interp{i}(aoa_domain_col, re_query_array);
    interp_fit = fit(aoa_domain_col, cl_at_re, 'linearinterp');
    
    % Find the zero-lift angle of attack
    alpha_0 = fzero(interp_fit, 0);
    
    % Calculate the lift-curve slope in the linear region near alpha_0
    slope_range_deg = linspace(alpha_0, alpha_0 + 2, 20);
    slope_range_rad = deg2rad(slope_range_deg);
    cl_for_slope = interp_fit(slope_range_deg);
    p = polyfit(slope_range_rad, cl_for_slope, 1);
    Cla = p(1);
else
    % --- Use pre-calculated global values for XFOIL ---
    Cla = settings.Cla(i);
    alpha_0 = settings.alpha_0(i);
end

%% ITERATIVE SOLVER FOR INFLOW ANGLE (phi)
% This loop iteratively solves for the inflow angle 'phi' that balances
% the geometric angle of attack with the induced velocities.
phi_1 = 0; phi_2 = 150; itr = 1; max_itr = 100;
lambda_y = 0; % Initialize radial flow angle

while abs(phi_1-phi_2) > 1.0e-6
    phi = (phi_1+phi_2)/2;
    b_alpha = beta_i - phi;
    b_alpha_y = b_alpha * cosd(lambda_y); % Correct AoA for radial flow

    %% STEP 1: GET BASE 2D AERO COEFFICIENTS
    if use_neuralfoil
        % For NeuralFoil, interpolate on the pre-computed grid.
        V_rel_approx = sqrt(Ua^2 + Ut_pr^2);
        reynolds_number = (settings.rho * V_rel_approx * c(i)) / settings.nu;
        
        % Clamp Reynolds number to the pre-computed grid bounds
        re_min = settings.nf_re_min; re_max = settings.nf_re_max;
        reynolds_number_clamped = max(re_min, min(re_max, reynolds_number));

        CL_2D = settings.CL_NF_Interp{i}(b_alpha, reynolds_number_clamped);
        CD_2D = settings.CD_NF_Interp{i}(b_alpha, reynolds_number_clamped);
        CD_0 = settings.CL_NF_Interp{i}(0, reynolds_number_clamped);
    else
        % For XFOIL, use the fitted curves from the loaded .xlsx data.
        CL_2D = CL_P{i}(b_alpha);
        CD_2D = CD_P{i}(b_alpha_y);
        CD_0 = CD_P{i}(alpha_0);
    end

    %% STEP 2: APPLY 3D ROTATIONAL CORRECTIONS (Du & Selig Model)
    if Use_3D_polar == true && x_blade(i) < 0.8
        J_loc = 2*pi*Ua/(omega*2*R+2*pi*UT);
        R0 = 1/(1+J_loc^2)*(x_blade(i)*R/c(i));
        fL = tanh(R0^-2)^3; fD = fL/2;

        dCl = Cla*deg2rad((b_alpha-alpha_0)) - CL_2D;
        dCd = CD_2D - CD_0;
        
        CL_final = CL_2D + fL*dCl;
        CD_final = CD_2D + fD*dCd;
    else
        CL_final = CL_2D;
        CD_final = CD_2D;
    end
    
    %% STEP 3: CALCULATE INDUCED FLOW AND CHECK CONVERGENCE
    % Calculate axial (VA_W) and tangential (VT_W) induction factors.
    VA_W = sigmal/(4)*(CL_final*cosd(phi) - CD_final*sind(phi))*cscd(phi);
    VT_W = sigmal/(4)*(CL_final*sind(phi) + CD_final*cosd(phi))*cscd(phi);
    
    G = sind(phi) - VA_W;
    H = cosd(phi) + VT_W;
    
    % The function 'I_new' goes to zero when the correct phi is found.
    I_new = sind(phi)*(Ut_pr*G-Ua*H);   
    
    % Bisection method to converge on the correct phi
    if I_new < 0, phi_1 = phi; else, phi_2 = phi; end
    
    % Update radial flow angle for the next iteration
    V_rel = Ut_pr/H; 
    V_t = V_rel*cosd(phi);
    lambda_y = atand(UR/V_t);
    
    itr = itr+1;
    if itr > max_itr, warning('Inflow Iteration in phi_algorithm_ISAE exceeds limit itr=%d',itr); break; end
end

end
