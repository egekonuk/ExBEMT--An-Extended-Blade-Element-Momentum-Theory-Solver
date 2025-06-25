function [TT_final,TN_final,TQ_final,CTT_final,CNT_final,CQT_final]=BEMT_inflow_Averaged(V_inf_prll,V_inf_perp,alpha,beta,Cla,alpha_0, fidelityModel)
%BEMT_inflow_Averaged Performs BEMT analysis with optional dynamic inflow models.
%
%   Inputs:
%       V_inf_prll      - Parallel freestream velocity component
%       V_inf_perp      - Perpendicular freestream velocity component
%       alpha           - Angle of incidence (degrees)
%       beta            - Blade pitch angle distribution (degrees)
%       Cla             - Lift curve slope for each blade section
%       alpha_0         - Zero-lift angle of attack for each blade section
%       fidelityModel   - String specifying model: 'BEMT',
%                         'BEMT+ Dynamic Inflow Model', or
%                         'BEMT+ Q-Unsteady Dynamic Inflow Model'
%
%   Outputs:
%       TT_final, TN_final, TQ_final - Total Thrust, Normal Force, and Torque
%       CTT_final, CNT_final, CQT_final - Corresponding coefficients (from BEMT)

%% CODE STARTS
global B R omega rho k_w Np delta_x c x_blade CL_P CD_P AoA_data radial_correct Use_3D_polar Sel_Prop
%% initialization BEMT
n = omega/(2*pi);
sep = 2*pi/B;
Ua = V_inf_perp;
Uy = V_inf_prll;
J = pi*Ua/(omega*R);
dpsi = 1; %deg
%% Preallocation
phicheck = zeros(181,1); UT = zeros(181,1); UR = zeros(181,1); Ut_pr = zeros(181,1);
phi = zeros(181,1); b_alpha = zeros(181,1); V_rel = zeros(181,1); V_p = zeros(181,1);
V_t = zeros(181,1); w_a = zeros(181,1); w_t = zeros(181,1); DL = zeros(181,1); DD = zeros(181,1);
DCL = zeros(181,1); DCD = zeros(181,1); dCT = zeros(181,1); dT = zeros(181,1); dM = zeros(181,1);
d2n = zeros(181,1); d2M = zeros(181,1); dCN = zeros(181,1); dQ_r = zeros(181,1);
dN = zeros(181,1); dQ = zeros(181,1); lambda = zeros(181,1); lambda_i = zeros(181,1);
xi = zeros(181,1); xi_i = zeros(181,1); V_inf_eff = zeros(181,1); psi = zeros(181,1);

% Determine which inflow model to run
run_di_calcs = ~strcmp(fidelityModel, 'BEMT');
is_unsteady = strcmp(fidelityModel, 'BEMT+ Q-Unsteady Dynamic Inflow Model');
is_steady = strcmp(fidelityModel, 'BEMT+ Dynamic Inflow Model');

%% Main Loop
for bn = 1:B
    for i = 1:length(x_blade)
        sigmal(i) = B*c(i)/(2*pi*(x_blade(i)*R));
        j = 0;
        delta_time = deg2rad(dpsi)/omega;
        for t = 0:delta_time:2*pi/B/omega
            j = j+1;
            psi(j) = sep*(bn-1) + t*omega;
            UT(j) = Uy*sin(psi(j));
            if radial_correct, UR(j) = Uy*cos(psi(j)); else, UR(j) = 0; end
            Ut_pr(j) = omega*x_blade(i)*R+UT(j);
            
            [H,VA_W,VT_W,CL_3D_BEMT(j),CD_3D_BEMT(j),phi(j),b_alpha(j),b_alpha_y(j)] = phi_algorithm_ISAE(Ua,UT(j),UR(j),Ut_pr(j),beta(i),i,sigmal(i),Cla);
            
            if b_alpha(j) > max(AoA_data{i}) || b_alpha(j) < min(AoA_data{i})
                fprintf('HIGH BLADE ANGLE DETECTED!! @ BEMT. B_alpha = %g @blade = #%d, Vel = %g\n',b_alpha(j),i,Ua/cosd(alpha));
            end
            
            F(j) = (2/pi)*acos(exp(B*(x_blade(i)-1)/(2*x_blade(i)*sind(phi(j)))));
            V_rel(j) = Ut_pr(j)/H;
            V_p(j) = V_rel(j)*sind(phi(j));
            V_t(j) = V_rel(j)*cosd(phi(j));
            lambda_y(j) = atand(UR(j)/V_t(j));
            V_rel_y(j) =  sqrt(V_rel(j)^2+UR(j)^2);
            w_a(j) = V_rel(j)*VA_W;
            w_t(j) = V_rel(j)*VT_W;
            
            DL(j) = 0.5*rho*V_rel(j)^2*c(i)*CL_3D_BEMT(j);
            DD(j) = 0.5*rho*V_rel_y(j)^2*c(i)*CD_3D_BEMT(j);
            DCL(j) = 0.5*sigmal(i)*V_rel(j)^2/(Ut_pr(j))^2*CL_3D_BEMT(j);
            DCD(j) = 0.5*sigmal(i)*V_rel_y(j)^2/(Ut_pr(j))^2*CD_3D_BEMT(j);
            
            dT(j) = (DL(j)*cosd(phi(j)) - DD(j)*sind(phi(j)))*delta_x(i);
            dM(j) = -dT(j)*x_blade(i)*R*cos(psi(j));
            dn(j) = -dT(j)*x_blade(i)*R*sin(psi(j));
            d2M(j) = -dT(j)*(x_blade(i)*R)^2*cos(2*psi(j));
            d2n(j) = -dT(j)*(x_blade(i)*R)^2*sin(2*psi(j));
            dQ_r(j) =  ((DL(j)*sind(phi(j)) + DD(j)*cosd(phi(j)))*cosd(lambda_y(j))*delta_x(i));
            dN(j) = dQ_r(j)*sin(psi(j));
            dQ(j) = dQ_r(j)*x_blade(i)*R;
            dCT(j) = (DCL(j)*cosd(phi(j)) - DCD(j)*sind(phi(j)))*delta_x(i)/R;
            dCQ_r(j) = (DCL(j)*sind(phi(j)) + DCD(j)*cosd(phi(j)))*delta_x(i)/R;
            dCN(j) = dCQ_r(j)*sin(psi(j));
            dCQ(j) = dCQ_r(j)*x_blade(i);
        end
        
        % Store BEMT results for contour plotting
        psi_BEMT{i,bn} = psi;
        dT_BEMT{i,bn} = dT;
        dQ_BEMT{i,bn} = dQ;
        Vel_BEMT{i,bn} = V_rel_y;
        AoA_BEMT{i,bn} = b_alpha_y;
        
        % BEMT Averages
        T_avg(i,bn) = B/(2*pi)*sum(dT)*Np*deg2rad(dpsi);
        M_avg(i,bn) = B/(2*pi)*sum(dM)*Np*deg2rad(dpsi);
        n_avg(i,bn) = B/(2*pi)*sum(dn)*Np*deg2rad(dpsi);
        N_avg(i,bn) = B/(2*pi)*sum(dN)*Np*deg2rad(dpsi);
        Q_avg(i,bn) = B/(2*pi)*sum(dQ)*Np*deg2rad(dpsi);
        M2_avg(i,bn) = B/(2*pi)*sum(d2M)*Np*deg2rad(dpsi);
        n2_avg(i,bn) = B/(2*pi)*sum(d2n)*Np*deg2rad(dpsi);
        CT_avg(i,bn) = B/(2*pi)*sum(dCT)*Np*deg2rad(dpsi);
        CN_avg(i,bn) = B/(2*pi)*sum(dCN)*Np*deg2rad(dpsi);
        CQ_avg(i,bn) = B/(2*pi)*sum(dCQ)*Np*deg2rad(dpsi);
        wa_avg = B/(2*pi)*sum(w_a)*Np*deg2rad(dpsi);
        
        if run_di_calcs
            if is_steady
                [~,nu_coeff{i}] = Peters_dynamic_inflow(T_avg(i,bn),n_avg(i,bn),M_avg(i,bn),n2_avg(i,bn),M2_avg(i,bn),Ua,Uy,wa_avg);
                nu_0 = nu_coeff{i}(1); nu_s = nu_coeff{i}(2); nu_c = nu_coeff{i}(3); nu_2s = nu_coeff{i}(4); nu_2c = nu_coeff{i}(5);
            end
            
            jj = j;
            for j_di = 1:jj
                if is_unsteady
                    if j_di == 1, diff=1; dTs=[dT(j_di) dT(j_di+1)]; dns=[dn(j_di) dn(j_di+1)]; dMs=[dM(j_di) dM(j_di+1)]; d2ns=[d2n(j_di) d2n(j_di+1)]; d2Ms=[d2M(j_di) d2M(j_di+1)]; w_as=[w_a(j_di) w_a(j_di+1)];
                    elseif j_di == jj, diff=2; dTs=[dT(j_di-1) dT(j_di)]; dns=[dn(j_di-1) dn(j_di)]; dMs=[dM(j_di-1) dM(j_di)]; d2ns=[d2n(j_di-1) d2n(j_di)]; d2Ms=[d2M(j_di-1) d2M(j_di)]; w_as=[w_a(j_di-1) w_a(j_di)];
                    else, diff=2; dTs=[dT(j_di-1) dT(j_di) dT(j_di+1)]; dns=[dn(j_di-1) dn(j_di) dn(j_di+1)]; dMs=[dM(j_di-1) dM(j_di) dM(j_di+1)]; d2ns=[d2n(j_di-1) d2n(j_di) d2n(j_di+1)]; d2Ms=[d2M(j_di-1) d2M(j_di) d2M(j_di+1)]; w_as=[w_a(j_di-1) w_a(j_di) w_a(j_di+1)];
                    end
                    [~,nu_coeff_unsteady] = Peters_dynamic_unsteady_inflow(dTs,dns,dMs,d2ns,d2Ms,Ua,Uy,w_as,diff,delta_time*omega);
                    nu_0 = nu_coeff_unsteady(1); nu_s = nu_coeff_unsteady(2); nu_c = nu_coeff_unsteady(3); nu_2s = nu_coeff_unsteady(4); nu_2c = nu_coeff_unsteady(5);
                end
                
                psi_current = psi(j_di);
                w_a_final(j_di) = (nu_0 + nu_s*x_blade(i)*sin(psi_current) + nu_c*x_blade(i)*cos(psi_current) ...
                    + nu_2s*x_blade(i)^2*sin(2*psi_current) + nu_2c*x_blade(i)^2*cos(2*psi_current))*omega*R;
                
                V_t_P(j_di) = Ut_pr(j_di) - w_t(j_di);
                V_p_P(j_di) = Ua + w_a_final(j_di);
                phi_P(j_di) = atand(V_p_P(j_di)/V_t_P(j_di));
                b_alpha_P(j_di) =  beta(i) - phi_P(j_di);
                W_P(j_di) = sqrt(V_t_P(j_di)^2+V_p_P(j_di)^2);
                lambda_y_p(j_di) = atand(UR(j_di)/(V_t_P(j_di)));
                b_alpha_P_y(j_di) = b_alpha_P(j_di) * cosd(lambda_y_p(j_di));
                W_P_y(j_di) = sqrt(V_t_P(j_di)^2 + UR(j_di)^2 + V_p_P(j_di)^2);

                if b_alpha_P(j_di) > max(AoA_data{i}) || b_alpha_P(j_di) < min(AoA_data{i})
                     fprintf('HIGH BLADE ANGLE DETECTED!! @ DI. B_alpha = %g @blade = #%d, Vel = %g\n',b_alpha_P(j_di),i,Ua/cosd(alpha));
                end
                
                if Use_3D_polar && x_blade(i) < 0.8
                    J_loc = 2*pi*Ua/(omega*2*R+2*pi*UT(j_di)); R0 = 1/(1+J_loc^2)*(x_blade(i)*R/c(i));
                    fL = tanh(R0^-2)^3; fD = fL/2;
                    CD_0 = CD_P{i}(alpha_0(i));
                    dCl = Cla(i)*deg2rad((b_alpha_P(j_di)-alpha_0(i))) - CL_P{i}(b_alpha_P(j_di));
                    dCd = CD_P{i}(b_alpha_P_y(j_di))-CD_0;
                else, fL=0;fD=0;dCl=0;dCd=0;
                end
                CL_3D_DI(j_di) = CL_P{i}(b_alpha_P(j_di)) + fL*dCl;
                CD_3D_DI(j_di) = CD_P{i}(b_alpha_P_y(j_di)) + fD*dCd;
                
                DL_P(j_di) = 0.5*rho*W_P(j_di)^2*c(i)*CL_3D_DI(j_di);
                DD_P(j_di) = 0.5*rho*W_P_y(j_di)^2*c(i)*CD_3D_DI(j_di);
                dT_P(j_di) = (DL_P(j_di)*cosd(phi_P(j_di)) - DD_P(j_di)*sind(phi_P(j_di)))*delta_x(i);
                dQ_r_P(j_di) =  ((DL_P(j_di)*sind(phi_P(j_di)) + DD_P(j_di)*cosd(phi_P(j_di)))*delta_x(i))*cosd(lambda_y_p(j_di));
                dN_P(j_di) = dQ_r_P(j_di)*sin(psi_current);
                dQ_P(j_di) = dQ_r_P(j_di)*x_blade(i)*R;
            end
            
            % Store DI results for contour plotting
            dT_PITT{i,bn} = dT_P;
            dQ_PITT{i,bn} = dQ_P;
            Vel_PITT{i,bn} = W_P_y;
            AoA_PITT{i,bn} = b_alpha_P_y;
            
            % Overwrite BEMT averages with DI averages
            T_avg(i,bn) = B/(2*pi)*sum(dT_P)*Np*deg2rad(dpsi);
            N_avg(i,bn) = B/(2*pi)*sum(dN_P)*Np*deg2rad(dpsi);
            Q_avg(i,bn) = B/(2*pi)*sum(dQ_P)*Np*deg2rad(dpsi);
        end
    end
    
    TT_blade(bn) = sum(T_avg(:,bn));
    TN_blade(bn) = sum(N_avg(:,bn));
    TQ_blade(bn) = sum(Q_avg(:,bn));
    CT_blade(bn) = sum(CT_avg(:,bn));
    CN_blade(bn) = sum(CN_avg(:,bn));
    CQ_blade(bn) = sum(CQ_avg(:,bn));
end

%% Final Total Propeller sum
TT_final = sum(TT_blade);
TN_final = sum(TN_blade);
TQ_final = sum(TQ_blade);
CTT_final = sum(CT_blade);
CNT_final = sum(CN_blade);
CQT_final = sum(CQ_blade);

%% Assign Plotting Variables to Base Workspace
disp('Assigning results to base workspace for potential plotting...');
try
    assignin('base', 'dT_BEMT', dT_BEMT);
    assignin('base', 'dQ_BEMT', dQ_BEMT);
    assignin('base', 'Vel_BEMT', Vel_BEMT);
    assignin('base', 'AoA_BEMT', AoA_BEMT);
    if run_di_calcs
        assignin('base', 'dT_PITT', dT_PITT);
        assignin('base', 'dQ_PITT', dQ_PITT);
        assignin('base', 'Vel_PITT', Vel_PITT);
        assignin('base', 'AoA_PITT', AoA_PITT);
    end
    assignin('base', 'psi_BEMT', psi_BEMT);
    assignin('base', 'x_blade', x_blade);
    assignin('base', 'B', B);
    assignin('base', 'J_plot', J);
    assignin('base', 'alpha_plot', alpha);
catch ME
    disp('Warning: Could not assign variables to base workspace for plotting.');
    disp(ME.message);
end

%% Save all data
folderName = fullfile('Saved Data', 'BEMT_DATA');
if ~exist(folderName, 'dir'), mkdir(folderName); end
fileName = sprintf('%s_%.2g_%ddeg.mat', Sel_Prop, J, alpha);
filePath = fullfile(folderName, fileName);
save(filePath);
disp(['Workspace variables have been saved to ' filePath]);
end
