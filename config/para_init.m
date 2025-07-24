function [values] = para_init()
    % Author: Xidong Mu (base) | Modified by Kgomotjo Seipobi
    
    values.noise_dB = -90;% noise power in dB
    values.noise = 10^(values.noise_dB/10); 
    values.alpha_k_n=0.1;
    values.alpha_k_f=0.9;
    values.weights_n= 0.5; % weight for near user
    values.weights_f= 0.3; % weight for far user
    values.weights_c= 0.2; % weight for backscatter
    
    
    values.M = 8; % overall antennas
    values.RIS_size = [2,15]; % reflecting elements at RIS
    values.N = values.RIS_size(1)*values.RIS_size(2); % reflecting elements at RIS
    
    values.n = 1; % equivalent noise power
    values.K = 3; % user number
    values.P_max = 50; % maximum transmit power in W
    values.eta = 0.5; % backscatter coefficient
    values.R_min_f = 0.01; % minimum rate requirement in bps/Hz
    values.R_min_n = 0.01; % minimum backscatter rate for near user in bps/Hz
    values.R_c_min = 0.001; % minimum backscatter rate for far user in bps/Hz
    values.nu_n = 1; % weight for near user
    
    values.nu_f = 1; % weight for far user
    values.nu_c = 1; % weight for backscatter
    values.max_iter = 10; % maximum iterations
    values.tol = 0.00001; % convergence tolerance
    
    values.pathloss= @(d) 30 + 22*log10(d); % path loss with d in m\
    
    values.rician = 10^(3/10); % rician factor
    
    values.RIS_loc = [0,0,0]; %Location of RIS-RIS
    values.BS_loc = [3, 90, 0]; % Location of BS
    
    % user locations
    range_f = [4, 6]; % range of user locations in meters
    range_n = [1, 3]; % range of user locations in meters
    values.user_loc = zeros(values.K, 3);
    for i = 1:values.K
        if i==2
            values.user_loc(i,:) = [ (range_f(2)-range_f(1))*rand(1) + range_f(1), 180*rand(1), 180*rand(1)-90 ];
        else
            values.user_loc(i,:) = [ (range_n(2)-range_n(1))*rand(1) + range_n(1), 180*rand(1), 180*rand(1)-90 ];
        end
    end

    end
    