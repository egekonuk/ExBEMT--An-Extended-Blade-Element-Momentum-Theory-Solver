function [va_avg, nu_coeff] = Peters_dynamic_inflow(dT, dn, dM, d2n, dM2, Ua, Uy, va_avg_initial, omega, R, rho)
%PETERS_DYNAMIC_INFLOW Calculates steady-state induced inflow using Peters' finite state model.
% OPTIMIZED VERSION: omega, R, rho passed as arguments.

%% 1. Initial Setup
tolerance = 1.0e-8; max_iter = 200;

%% 2. Non-dimensionalize Aerodynamic Loads
CT = dT / (rho * (omega / (2 * pi))^2 * (2 * R)^4);
Cn = dn / (rho * (omega / (2 * pi))^2 * (2 * R)^5);
CM = dM / (rho * (omega / (2 * pi))^2 * (2 * R)^5);
C2n = d2n / (rho * (omega / (2 * pi))^2 * (2 * R)^6);
CM2 = dM2 / (rho * (omega / (2 * pi))^2 * (2 * R)^4);
load_coeffs = [CT; Cn; CM; C2n; CM2];

%% 3. Non-dimensionalize Flight Conditions
lambda = Ua / (omega * R);
mu = Uy / (omega * R);
nu = va_avg_initial / (omega * R);

%% 4. Iterative Inflow Solution
nu_old = zeros(5, 1);
for iter = 1:max_iter
    V_T = sqrt((lambda + nu)^2 + mu^2);
    V_massflow = (mu^2 + (2 * nu + lambda) * (nu + lambda)) / V_T;
    V_mat = diag([V_T, V_massflow, V_massflow, V_massflow, V_massflow]);
    
    if mu < 1e-6, chi = pi / 2; else, chi = atan(abs(lambda + nu) / mu); end
    [L, ~] = inflowgains(chi, 2, 'false');
    
    nu_new = L * (V_mat \ load_coeffs);
    nu = nu_new(1); %(0.5 * ([1 0 0 0 0] / L) * nu_new);
    
    if norm(nu_new - nu_old) < tolerance, break; end
    nu_old = nu_new;
    if iter == max_iter, warning('Peters inflow solver did not converge within %d iterations.', max_iter); end
end

%% 5. Finalize Outputs
nu_coeff = nu_old;
va_avg = nu * omega * R;

end