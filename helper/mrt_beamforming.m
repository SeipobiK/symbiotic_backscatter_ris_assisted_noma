function [w_k] = mrt_beamforming(para, H_n, H_f)
    % Modified MRT Beamforming for 2-Cluster NOMA System
    % Inputs:
    %   - para: struct containing system parameters
    %   - H_n: Channel matrix for near users [M x K] (M = BS antennas, K = users)
    %   - H_f: Channel matrix for far users [M x K]
    % Output:
    %   - w_k: Beamforming vectors [M x 2] (one per cluster)
    
    M = para.M;              % Number of BS antennas
    K = para.K;              % Number of clusters (should be 2)
    P_max = para.P_max;      % Max transmit power
    
    assert(K == 2, 'This function is designed for 2 clusters');
    
    % Initialize beamforming matrix
    w_k = zeros(M, K);
    
    % Power allocation factors (can be customized)
    alpha_n = 0.55;  % More power to near users
    alpha_f = 0.45;  % Less power to far users
    
    % Cluster 1 (Near users) - MRT
    H_n_combined = mean(H_n, 2);  % Average channel for near users
    w_k(:,1) = H_n_combined / norm(H_n_combined);
    w_k(:,1) = sqrt(P_max * alpha_n) * w_k(:,1);
    
    % Cluster 2 (Far users) - MRT
    H_f_combined = mean(H_f, 2);  % Average channel for far users
    w_k(:,2) = H_f_combined / norm(H_f_combined);
    w_k(:,2) = sqrt(P_max * alpha_f) * w_k(:,2);
    
    % Verify power constraints
    total_power = sum(diag(w_k' * w_k));
    if total_power > P_max * 1.01
        warning('Power constraint violated. Normalizing...');
        w_k = w_k * sqrt(P_max / total_power);
    end
end