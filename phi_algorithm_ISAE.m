function [H,VA_W,VT_W,CL_3D,CD_3D,phi,b_alpha,b_alpha_y] = phi_algorithm_ISAE(Ua,UT,UR,Ut_pr,beta,i,sigmal,Cla)
 global CL_P CD_P omega R c x_blade Use_3D_polar
%% ISAE BEMT phi algorihm
       phi_1 = 0; phi_2 = 150; % initial points
       itr = 1;
       lambda_y = 0;
       alpha_0 = fzero(CL_P{i},0);      
       while abs(phi_1-phi_2) > 1.0e-6
       phi = (phi_1+phi_2)/2;
       b_alpha =  beta - phi; % pitch-phi(bn)-p_alpha_ind(bn); % prop aoa
       b_alpha_y = b_alpha * cosd(lambda_y);
%% cl,cd 3D correction / New part from section E / Need some work
if Use_3D_polar == true && x_blade(i) < 0.8  %% according to snel 80% chord length
    J_loc = 2*pi*Ua/(omega*2*R+2*pi*UT);
    R0 = 1/(1+J_loc^2)*(x_blade(i)*R/c(i));
      fL = tanh(R0^-2)^3; fD = fL/2;
       CD_0 = CD_P{i}(alpha_0);
       dCl = Cla(i)*deg2rad((b_alpha-alpha_0)) - CL_P{i}(b_alpha); %0;
       dCd = CD_P{i}(b_alpha_y)-CD_0; %0;
       %
else, fL=0;fD=0;dCl=0;dCd=0;
end
       CL_3D = CL_P{i}(b_alpha) + fL*dCl;
       CD_3D = CD_P{i}(b_alpha_y) + fD*dCd;
        % Blade tip loss correction
%         F = 2/pi*acos(exp(B*(x_blade(i)*R-1)/(2*x_blade (i)*R*sind(phi))));
%         KT = (1-(1-F))*cosd(phi);
%         KP = (1-(1-F))*sind(phi);
%% Induced flow fractions       
       VA_W = sigmal/(4)*(CL_3D*cosd(phi) - CD_3D*sind(phi))*cscd(phi);
       VT_W = sigmal/(4)*(CL_3D*sind(phi) + CD_3D*cosd(phi))*cscd(phi);
       G = sind(phi) - VA_W;
       H = cosd(phi) + VT_W;
       I_new = sind(phi)*(Ut_pr*G-Ua*H);   
         if I_new < 0
             phi_1 = phi;
         else
             phi_2 = phi;
         end
                % radial flow correction
                V_rel = Ut_pr/H; 
                V_t = V_rel*cosd(phi); %tangential prop velocity (parallel)
                lambda_y = atand(UR/V_t);
                %
         itr = itr+1; % iterate as well
         if itr > 100
             warning('Inflow Iteration exceeds limit itr=%d',itr)
             break
         end
       end

%         disp(itr-1) 