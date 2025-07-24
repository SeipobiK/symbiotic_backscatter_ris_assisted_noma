function [w_k] = mrt_beamforming(para, H_k)
    % Inputs:
    %   - para: struct containing system parameters
    %   - H_k:  Channel matrix for user k (size [M, 1], where M = BS antennas)
    % Output:
    %   - w_k:  MRT beamforming vector (normalized to meet power constraints)
    
    M = para.M;              % Number of BS antennas
    P_max = para.P_max;      % Max transmit power
    
    % MRT beamforming (maximizes signal power)
    w_k = H_k / norm(H_k);   % Unit-norm direction
    
    % Scale to meet power constraint
    w_k = sqrt(P_max) * w_k; % w_k^H * w_k = P_max
end