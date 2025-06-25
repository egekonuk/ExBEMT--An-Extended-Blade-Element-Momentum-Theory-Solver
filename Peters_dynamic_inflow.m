function [va_avg,nu_coeff] = Peters_dynamic_inflow(dT,dn,dM,d2n,d2M,Ua,Uy,va_avg)
global omega R rho
       nu_0 = 0; nu_s = 0; nu_c = 0; nu_2s = 0; nu_2c = 0; % initialize
       nu_set{1} = zeros(5,1);
       err1 = 100; err2 = 100; err3 = 100;  err4 = 100;  err5 = 100;
       iter = 1;
%        Co_Mat{itr} = zeros(3,1);
       %initialize force/moment coefficients
       CT = dT/(rho*(omega/(2*pi))^2*(2*R)^4);
       Cn = dn/(rho*(omega/(2*pi))^2*(2*R)^5);
       CM = dM/(rho*(omega/(2*pi))^2*(2*R)^5);
       C2n = d2n/(rho*(omega/(2*pi))^2*(2*R)^6);
       C2M = d2M/(rho*(omega/(2*pi))^2*(2*R)^4);
       Co_Mat = [CT;Cn;CM;C2n;C2M];
        %% Non-dimensionalize 
        lambda = Ua/(omega*R);
        mu = Uy/(omega*R);
        nu = va_avg/(omega*R);
        lambda_T = sqrt(Ua^2+Uy^2)/(omega*R);

       while abs(err1(end))+abs(err2(end))+abs(err3(end))+abs(err4(end))+abs(err5(end)) > 1.0e-8 %  check convergence 
           
%           inflow ratio
          V_T = sqrt((lambda+nu)^2+mu^2);         
%           mass-flow parameter
          V = (mu^2+(2*nu+lambda)*(nu+lambda))/V_T;
%           V matrix
          V_mat = [V_T 0 0 0 0; 0 V 0 0 0; 0 0 V 0 0; 0 0 0 V 0; 0 0 0 0 V];
          
          % skew wake angle
          chi = atan(abs(lambda+nu)/mu);

%            V_T = sqrt((lambda_T+nu)^2*sin(chi)^2+((lambda_T+nu)*cos(chi)+nu)^2);
%            V = (lambda_T+nu)^2*sin(chi)^2+((lambda_T+nu)*cos(chi)+nu)*((lambda_T+nu)*cos(chi)+2*nu)/V_T;
%            V_mat = [V_T 0 0 0 0; 0 V 0 0 0; 0 0 V 0 0; 0 0 0 V 0; 0 0 0 0 V];

          % effective skew angle He[1989]
  %        chi = atan(((lambda_T+nu)*sin(chi))/((lambda_T+nu)*cos(chi)+nu)); 
% 
         % % Gain matrix L 
         %   LF = sqrt((1-sin(chi))/(1+sin(chi)));
         % 
         %   LF1 = 15*pi/64*LF;
         %   LF2 = 4/(1+sin(chi));
         % 
         %   LF44 = -sin(chi)*(11-5*sin(chi))/(1+sin(chi));
         %   LF24 = 105*pi/128*LF;
         %   LF42 = -2205*pi/2048*LF;
         % 
         %   LF55 = -6*(1+sin(chi)^2)/(1+sin(chi))^2;
         %   LF51 = -3/7*LF;
         %   LF53 = -3*pi/4*sin(chi)*(1-sin(chi));
         %   LF35 = 2*sin(chi)*(1-sin(chi));
         % 
         %   L = [0.5 0       -LF1         0     LF51;...
         %        0   LF2     0            LF24  0;...
         %        LF1 0       LF2*sin(chi) 0     LF35;...
         %        0   LF42   0            LF44     0;...
         %        0   0       LF53            0     LF55]; % L gain matrix  5x5    

           [L,~] = inflowgains(chi,2,'false');

           %induced flow matrix solution
           nu_set{iter+1} = L*inv(V_mat)*Co_Mat; 
           %coefficents
           nu_0(end+1) = nu_set{iter+1}(1); 
           nu_s(end+1) = nu_set{iter+1}(2);
           nu_c(end+1) = nu_set{iter+1}(3);
           nu_2s(end+1) = nu_set{iter+1}(4);
           nu_2c(end+1) = nu_set{iter+1}(5);
           %new va_avg( axial induced flow)
           nu = (1/2*[1 0 0 0 0]*inv(L)*nu_set{iter+1}); % peters1988- Eq28     %update the inflow velocity
           va_avg = nu*omega*R;  % this gives different solution (should be same??)

%             CT = 2*nu*sqrt(mu^2+(lambda+nu)^2);
            
          err1 = diff(nu_0); err2 = diff(nu_s); err3 = diff(nu_c); err4 = diff(nu_2s); err5 = diff(nu_2c);% error
            % new iteration 
          iter = iter+1;
       end


                  % axial induced flow distribution
           nu_coeff = nu_set{end};
%            nu = nu_0(end) + nu_s(end)*x_blade.*(sin(psi))+ nu_c(end)*x_blade.*(cos(psi));
        