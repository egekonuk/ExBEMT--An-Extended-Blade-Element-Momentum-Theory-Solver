function [CL_ext, CD_ext] = extend_polar(alpha_rad, alpha_stall_pos_rad, CL_stall_pos, CD_stall_pos, alpha_stall_neg_rad, CL_stall_neg, CD_stall_neg, CD_max)
% EXTEND_POLAR Extrapolates airfoil data to +/- 180 degrees using the Viterna(1982) method.
%   This method is more robust than pure curve fitting as it is anchored by
%   a physical value for the maximum drag coefficient (flat plate drag).
%
%   Inputs:
%   alpha_rad           - Current angle of attack in radians.
%   alpha_stall_pos_rad - Positive stall angle from original data (radians).
%   CL_stall_pos        - CL at the positive stall angle.
%   CD_stall_pos        - CD at the positive stall angle.
%   alpha_stall_neg_rad - Negative stall angle from original data (radians).
%   CL_stall_neg        - CL at the negative stall angle.
%   CD_stall_neg        - CD at the negative stall angle.
%
%   Outputs:
%   CL_ext              - Extended Lift Coefficient.
%   CD_ext              - Extended Drag Coefficient.

    % Assume a flat plate drag coefficient at 90 degrees. 
    % CD_max = 2;

    % The Viterna equations are applied differently for positive and negative stall.
    if alpha_rad > 0
        % Use positive stall point for blending
        alpha_s = alpha_stall_pos_rad;
        CL_s = CL_stall_pos;
        CD_s = CD_stall_pos;
    else
        % Use negative stall point for blending
        alpha_s = alpha_stall_neg_rad;
        CL_s = CL_stall_neg;
        CD_s = CD_stall_neg;
    end

    % --- Calculate Viterna Blending Coefficients ---
    % These coefficients ensure a smooth C1 continuous blend from the stall
    % point to the high-angle sinusoidal model.
    
    % Avoid division by zero if stall angle is exactly 90 degrees
    if abs(abs(alpha_s) - pi/2) < 1e-6
        cos_alpha_s = 1e-6;
    else
        cos_alpha_s = cos(alpha_s);
    end

    B2 = (CD_s - CD_max * (sin(alpha_s)^2)) / cos_alpha_s;
    
    % The lift coefficient equation has a singularity at alpha = 0 and 180.
    % The model is only valid post-stall, so we check if sin(alpha_s) is near zero.
    if abs(sin(alpha_s)) < 1e-6
        A2 = 0; % No lift contribution from this term if stalled at 0 AoA
    else
        A2 = (CL_s - (CD_max / 2) * sin(2 * alpha_s)) * (sin(alpha_s) / (cos_alpha_s^2));
    end
    
    % --- Apply the Viterna Equations ---
    
    % Add a small epsilon to sin(alpha_rad) to avoid division by zero at alpha = +/-180 deg
    sin_alpha_rad_safe = sin(alpha_rad);
    if abs(sin_alpha_rad_safe) < 1e-9
       sin_alpha_rad_safe = 1e-9 * sign(alpha_rad);
    end
    
    CL_ext = (CD_max / 2) * sin(2 * alpha_rad) + A2 * (cos(alpha_rad)^2) / sin_alpha_rad_safe;
    CD_ext = CD_max * (sin(alpha_rad)^2) + B2 * cos(alpha_rad);

end