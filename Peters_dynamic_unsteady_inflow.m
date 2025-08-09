function [w_a_new, nu_coeff] = Peters_dynamic_unsteady_inflow(dT, dn, dM, d2n, dM2, Ua, Uy, w_a, j, delta_time, omega, R, rho)
%PETERS_DYNAMIC_UNSTEADY_INFLOW Calculates induced velocity using Peters' unsteady model.
% CORRECTED AND ROBUST VERSION with debugging.
tolerance = 1.0e-7;
max_iter = 100;
relaxation_factor = 0.1; % Use a small value for stability

% --- Stagnation Detection Parameters ---
stagnation_patience = 30; % How many iterations to wait before checking for stagnation
stagnation_tolerance = 1.0e-5; % How small a change is considered "stuck"

best_norm = inf;
stagnation_counter = 0;

%% 2. Approximate Time Derivatives of Inflow Coefficients
% This logic remains the same for estimating the derivative term.
num_points = length(w_a);
nu_set_initial = zeros(5, num_points);

lambda = Ua / (omega * R);
mu = Uy / (omega * R);

for i = 1:num_points
    nu_i = w_a(i) / (omega * R);
    
    load_coeffs_hist = [
        dT(i) / (rho * (omega/(2*pi))^2 * (2*R)^4);
        dn(i) / (rho * (omega/(2*pi))^2 * (2*R)^5);
        dM(i) / (rho * (omega/(2*pi))^2 * (2*R)^5);
        d2n(i) / (rho * (omega/(2*pi))^2 * (2*R)^6);
        dM2(i) / (rho * (omega/(2*pi))^2 * (2*R)^4);
    ];

    V_T_hist = sqrt((lambda + nu_i)^2 + mu^2);
    V_massflow_hist = (mu^2 + (2*nu_i + lambda) * (nu_i + lambda)) / V_T_hist;
    V_mat_initial = diag([V_T_hist, V_massflow_hist, V_massflow_hist, V_massflow_hist, V_massflow_hist]);

    if mu < 1e-6, chi_hist = pi/2; else, chi_hist = atan(abs(lambda + nu_i) / mu); end
    [L_initial, ~] = inflowgains(chi_hist, 2, 'false');
    nu_set_initial(:, i) = L_initial * (V_mat_initial \ load_coeffs_hist);
end

if num_points == 2
    nu_derivatives = (nu_set_initial(:, 2) - nu_set_initial(:, 1)) / delta_time;
elseif num_points == 3
    nu_derivatives = (nu_set_initial(:, 3) - nu_set_initial(:, 1)) / (2 * delta_time);
else
    error('Input history for unsteady inflow must have 2 or 3 points.');
end

%% 3. Iterative Solution for Current Time Step
load_coeffs_current = [
    dT(j) / (rho * (omega/(2*pi))^2 * (2*R)^4);
    dn(j) / (rho * (omega/(2*pi))^2 * (2*R)^5);
    dM(j) / (rho * (omega/(2*pi))^2 * (2*R)^5);
    d2n(j) / (rho * (omega/(2*pi))^2 * (2*R)^6);
    dM2(j) / (rho * (omega/(2*pi))^2 * (2*R)^4);
];

nu = w_a(j) / (omega * R); % Initial guess for scalar average inflow
nu_old = zeros(5, 1);     % Initial guess for previous coefficient vector

for iter = 1:max_iter
    % Calculate physics based on the current best guess for 'nu'
    V_T = sqrt((lambda + nu)^2 + mu^2);
    if V_T < 1e-8; V_T = 1e-8; end % Prevent division by zero
    
    V_massflow = (mu^2 + (2*nu + lambda) * (nu + lambda)) / V_T;
    V_mat = diag([V_T, V_massflow, V_massflow, V_massflow, V_massflow]);
    
    if mu < 1e-6, chi = pi/2; else, chi = atan(abs(lambda + nu) / mu); end
    [L, M] = inflowgains(chi, 2, 'false');
    
    % Solve for a trial solution of the new inflow coefficient vector
    nu_new_trial = L * (V_mat \ (load_coeffs_current - M * nu_derivatives));
    
    % =================== CORRECTED & ROBUST UPDATE ===================
    % Apply relaxation to the entire vector to dampen oscillations
    nu_new = (1 - relaxation_factor) * nu_old + relaxation_factor * nu_new_trial;
    
    % Update the scalar 'nu' for the next iteration using the most recent value
    nu = nu_new(1);
    % =================================================================
    
    % Check for convergence
    convergence_norm = norm(nu_new - nu_old);
    if convergence_norm < tolerance, break; end
    
    nu_old = nu_new; % Update for the next loop
    
    % =================== STAGNATION CHECK ===================
    if iter > stagnation_patience
        if abs(convergence_norm - best_norm) < stagnation_tolerance
            stagnation_counter = stagnation_counter + 1;
        else
            stagnation_counter = 0; % Reset if there's progress
        end
        best_norm = min(best_norm, convergence_norm);
        
        if stagnation_counter >= 10 % If stuck for 10 iterations
            % disp('Unsteady inflow solver stagnated. Switching to steady solution for this step.');
            [w_a_new, nu_coeff] = Peters_dynamic_inflow(dT(j), dn(j), dM(j), d2n(j), dM2(j), Ua, Uy, w_a(j), omega, R, rho);
            return;
        end
    end
    % =========================================================

    if iter == max_iter
        disp('Unsteady inflow solver did not converge. Switching to steady solution.');
        [w_a_new, nu_coeff] = Peters_dynamic_inflow(dT(j), dn(j), dM(j), d2n(j), dM2(j), Ua, Uy, w_a(j), omega, R, rho);
        return;
    end
end

%% 4. Finalize Outputs
nu_coeff = nu_new;
w_a_new = nu * omega * R;

end