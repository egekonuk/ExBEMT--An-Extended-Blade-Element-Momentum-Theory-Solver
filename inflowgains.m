function [L,M] = inflowgains(chi,m,alter)

M1 = 8/(3*pi); M2 = -16/(45*pi); M3 = M2; M4 = -256/(1575*pi); M5 = M4; % peters 1980 (-) sign converted to plus

X = tan(chi/2);

if alter == true

    L = [0.5 0   -15*pi/64*X;...
    0    2*(1+X^2) 0 ;...       
    15*pi/64*X  0   2*(1-X^2)]; % L gain matrix  3x3  

    M = [M1 0   0;...
        0   M2  0;...       
        0   0   M3]; % M apparent mass 3x3  
    return
end

if m == 1

LF = sqrt((1-sin(chi))/(1+sin(chi)));
LF1 = 15*pi/64*LF;
LF2 = 4/(1+sin(chi));

L = [0.5 0   -LF1;...
    0    LF2 0; ...       
    LF1  0   LF2*sin(chi)]; % L gain matrix  5x5  
M = [M1 0   0 ;...
     0   M2  0 ;...       
     0   0   M3]; % M apparent mass 3x3  
elseif m ==2

LF = sqrt((1-sin(chi))/(1+sin(chi)));
LF13 = -15*pi/64*LF;
LF31 = 3*pi/8*LF;
LF2 = -4/(1+sin(chi));

LF44 = -sin(chi)*(11-5*sin(chi))/(1+sin(chi));
LF24 = -105*pi/128*LF;
LF42 = 45*pi/32*LF;

LF55 = -6*(1+sin(chi)^2)/(1+sin(chi))^2;
LF51 = 3/5*LF;
LF53 = pi*sin(chi)*(1-sin(chi));
LF35 = -2*sin(chi)*(1-sin(chi));

L = [0.5   0      LF13           0     0;...
    0      LF2    0             LF24  0;...
    LF31    0      LF2*sin(chi)  0     LF35;...
    0      LF42   0             LF44  0;...
    LF51   0      LF53          0     LF55]; % L gain matrix  5x5
M = [M1 0   0  0   0; ...
    0   M2  0  0   0; ...       
    0   0   M3 0   0; ...
    0   0   0  M4  0; ...
    0   0   0  0   M5]; % M apparent mass 5x5  
else
    error('No L matrix formation exist for given harmonics')
end