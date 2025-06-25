function [w_a_new, nu_coeff] = Peters_dynamic_unsteady_inflow(dT, dn, dM, d2n, d2M, Ua, Uy, w_a, j, delta_time)
% Peters_dynamic_unsteady_inflow Calculates the new axial induced velocity
% using Peters' finite second-order state dynamic inflow model.
%
% This function iteratively solves for the components of the induced flow coefficients (nu)
% based on the aerodynamic loads on the propeller blade.

% Declare global variables from the main BEMT script
global omega R rho

%% 1. Initialization
% Initialize inflow coefficient sets and iteration counter.
nu_set = {zeros(5,1)};
iter = 1;

%% 2. Non-dimensionalization of Constant Parameters
% Convert dimensional flight conditions into non-dimensional ratios.
lambda = Ua / (omega * R); % Inflow ratio (axial)
mu = Uy / (omega * R);     % Advance ratio (tangential)
lambda_T = sqrt(Ua^2 + Uy^2) / (omega * R); % Total advance ratio

%% 3. Finite Difference Derivative Calculation
% Calculate the time derivatives of the induced flow coefficients using
% finite differences on the provided aerodynamic loads.

for i = 1:length(w_a)
    % Convert the current induced velocity to a non-dimensional coefficient
    nu_induced_initial(i) = w_a(i) / (omega * R);
    
    % Non-dimensionalize the aerodynamic loads for the current time step
    CT = dT(i) / (rho * (omega/(2*pi))^2 * (2*R)^4);
    Cn = dn(i) / (rho * (omega/(2*pi))^2 * (2*R)^5);
    CM = dM(i) / (rho * (omega/(2*pi))^2 * (2*R)^5);
    C2n = d2n(i) / (rho * (omega/(2*pi))^2 * (2*R)^6);
    C2M = d2M(i) / (rho * (omega/(2*pi))^2 * (2*R)^4);
    Coeffs_initial{i} = [CT; Cn; CM; C2n; C2M];
    
    % Calculate flow parameters based on this initial induced flow
    chi_initial = atan(abs(lambda + nu_induced_initial(i)) / mu);  %wake skew
    chi_initial_eff = atan((lambda_T*sin(chi_initial))/(lambda_T*cos(chi_initial)+nu_induced_initial(i))); % effective skew angle He[1989]
    V_T_initial = sqrt((lambda + nu_induced_initial(i))^2 + mu^2);  %inflow ratio
    V_massflow_initial = (mu^2 + (2*nu_induced_initial(i) + lambda) * (nu_induced_initial(i) + lambda)) / V_T_initial; %mass-flow parameters
%Alternatives
    %    V_T = sqrt(lambda_T^2*sin(chi_d)^2+(lambda_T*cos(chi_d)+nu_d(i))^2);
    %    V_massflow = lambda_T^2*sin(chi_d)^2+(lambda_T*cos(chi_d)+nu_d(i))*(lambda_T*cos(chi_d)+2*nu_d(i))/V_T_d;



    % Get the gain matrix L based on the initial skew angle
    [L_initial, ~] = inflowgains(chi_initial_eff, 2, 'false');
    
    % Form the V matrix
    V_mat_initial = diag([V_T_initial, V_massflow_initial, V_massflow_initial, V_massflow_initial, V_massflow_initial]);
    
    % Calculate the initial set of induced flow coefficients
    nu_set_initial{i} = L_initial * inv(V_mat_initial) * Coeffs_initial{i};
end

% Approximate time derivatives using central or forward/backward differencing
if i == 2 % Forward/backward difference for 2 points
    nu0_star = (nu_set_initial{2}(1) - nu_set_initial{1}(1)) / delta_time;
    nus_star = (nu_set_initial{2}(2) - nu_set_initial{1}(2)) / delta_time;
    nuc_star = (nu_set_initial{2}(3) - nu_set_initial{1}(3)) / delta_time;
    nu2s_star = (nu_set_initial{2}(4) - nu_set_initial{1}(4)) / delta_time;
    nu2c_star = (nu_set_initial{2}(5) - nu_set_initial{1}(5)) / delta_time;
elseif i == 3 % Central difference for 3 points
    nu0_star = (nu_set_initial{3}(1) - nu_set_initial{1}(1)) / (2 * delta_time);
    nus_star = (nu_set_initial{3}(2) - nu_set_initial{1}(2)) / (2 * delta_time);
    nuc_star = (nu_set_initial{3}(3) - nu_set_initial{1}(3)) / (2 * delta_time);
    nu2s_star = (nu_set_initial{3}(4) - nu_set_initial{1}(4)) / (2 * delta_time);
    nu2c_star = (nu_set_initial{3}(5) - nu_set_initial{1}(5)) / (2 * delta_time);
end
nu_derivatives = [nu0_star; nus_star; nuc_star; nu2s_star; nu2c_star];

%% 4. Iterative Inflow Solution
% Iteratively solve for the final induced flow coefficients at the current time step `j`.

% Initialize force/moment coefficients at the current time step
CT = dT(j) / (rho * (omega/(2*pi))^2 * (2*R)^4);
Cn = dn(j) / (rho * (omega/(2*pi))^2 * (2*R)^5);
CM = dM(j) / (rho * (omega/(2*pi))^2 * (2*R)^5);
C2n = d2n(j) / (rho * (omega/(2*pi))^2 * (2*R)^6);
C2M = d2M(j) / (rho * (omega/(2*pi))^2 * (2*R)^4);
Co_Mat = [CT; Cn; CM; C2n; C2M];

% Current calue for the non-dimensional induced velocity
nu = w_a(j) / (omega * R);

% Loop until convergence
Total_err = 100;
while Total_err > 1.0e-6
    
    % Calculate flow parameters for the current iteration
    V_T = sqrt((lambda + nu)^2 + mu^2);
    V_massflow = (mu^2 + (2*nu + lambda) * (nu + lambda)) / V_T;
    %Alternatives
    %    V_T = sqrt(lambda_T^2*sin(chi_d)^2+(lambda_T*cos(chi_d)+nu_d(i))^2);
    %    V_massflow = lambda_T^2*sin(chi_d)^2+(lambda_T*cos(chi_d)+nu_d(i))*(lambda_T*cos(chi_d)+2*nu_d(i))/V_T_d;

    chi = atan(abs(lambda + nu) / mu);
    chi_eff = atan((lambda_T*sin(chi))/(lambda_T*cos(chi)+nu)); % effective skew angle He[1989]
    % Get gain matrices L and M
    [L, M] = inflowgains(chi_eff, 2, 'false');
    
    % Form the V matrix
    V_mat = diag([V_T, V_massflow, V_massflow, V_massflow, V_massflow]);
    
    % Solve the linear system for the next set of inflow coefficients
    % Suggestion: For better performance and stability, consider using the
    % backslash operator: nu_set{iter+1} = L * (V_mat \ (Co_Mat - M * nu_derivatives));
    nu_set{iter+1} = L * inv(V_mat) * (Co_Mat - M * nu_derivatives);
    
    % Update the average axial induced velocity based on Peters 1988, Eq. 28
    % Suggestion: For better performance, consider nu = (0.5 * L(1,:) * nu_set{iter+1});
    nu = (0.5 * [1 0 0 0 0] * inv(L) * nu_set{iter+1});
    
    % Check for convergence
    err = nu_set{iter+1} - nu_set{iter};
    Total_err = sum(abs(err));
    
    % Increment iteration counter and handle non-convergence (Switch to
    % Steady Soution)
    iter = iter + 1;
    if iter > 200
        warning("\n Dyn. Inflow iteration exceeds 200. Residual = %f \n Switching to the Steady solution!! \n", Total_err)
        [w_a_new, nu_set{end}] = Peters_dynamic_inflow(dT(j), dn(j), dM(j), d2n(j), d2M(j), Ua, Uy, w_a(j));
        break
    end
end

% Set the final output values
nu_coeff = nu_set{end};
w_a_new = nu * omega * R;

end
