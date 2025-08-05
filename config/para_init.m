function [values] = para_init()
    % Author: Xidong Mu (base) | Modified by Kgomotjo Seipobi
    
    values.noise_dB = -120;% noise power in dB
    values.noise = 10^(values.noise_dB/10); 
    values.alpha_k_n=0.1;
    values.alpha_k_f=0.9;
    values.weights_n= 0.5; % weight for near user
    values.weights_f= 0.3; % weight for far user
    values.weights_c= 0.2; % weight for backscatter
    values.K_u=3; % number of users in each cluster
    
    
    values.M = 4; % overall antennas
    values.RIS_size = [2,15]; % reflecting elements at RIS
    values.N = values.RIS_size(1)*values.RIS_size(2); % reflecting elements at RIS
    
    values.n = 1; % equivalent noise power
    values.K = 2; % user number
    values.P_max = 100; % maximum transmit power in W
    values.eta = 0.5; % backscatter coefficient
    values.R_min_f = 0.01; % minimum rate requirement in bps/Hz
    values.R_min_n = 0.01; % minimum backscatter rate for near user in bps/Hz
    values.R_c_min = 0.001; % minimum backscatter rate for far user in bps/Hz
    values.nu_n = 1; % weight for near user
    
    values.nu_f = 1; % weight for far user
    values.nu_c = 1; % weight for backscatter
    values.max_iter = 10; % maximum iterations
    values.tol = 0.00001; % convergence tolerance
    
    values.pathloss= @(d) 30 + 22*log10(d); % path loss with d in m
    
    values.rician = 10^(3/10); % rician factor
    
    values.RIS_loc = [0,0,0]; %Location of RIS-RIS
    values.BS_loc = [3, 90, 0]; % Location of BS

    angles = calculate_angles(); % Calculate angles for users

    values.BS_loc = [angles.BS_AoD.distance, ...       % D
                     angles.BS_AoD.elevation, ...     % ELV
                     angles.BS_AoD.azimuth];          % AZ
    
    % RIS Parameters (AoA - reciprocal of BS AoD)
    values.RIS_loc = [angles.RIS_AoA.distance, ...     % D (same as BS-RIS distance)
                      angles.RIS_AoA.elevation, ...   % ELV (negated BS elevation)
                      angles.RIS_AoA.azimuth];        % AZ
        % Cluster 1 users [User1, User2, BD]
    values.users_cluster1 = [27, 4, 0; 
                      32, 0, 0; 
                      30, 5, 0];
    
    % Cluster 2 users [User1, User2, BD]
    values.users_cluster2 = [43, 4, 0; 
                      37, 0, 0; 
                      40, 5, 0];

    % Cluster 1 distances
    values.BD_cluster1 = [30, 5, 0];
    values.User1_cluster1 = [27, 4, 0];
    values.User2_cluster1 = [32, 0, 0];

    % Cluster 2 distances
    values.BD_cluster2 = [40, 5, 0];
    values.User1_cluster2 = [43, 4, 0];
    values.User2_cluster2 = [37, 0, 0];


    values.userloc = zeros(3, 2, 3);  % 3 users, 2 clusters, 3 parameters

    
    % Cluster 1 (now cluster index 1 in 2nd dimension)
    for i = 1:3
        values.userloc(i, 1, 1) = angles.RIS_AoD.cluster1(i).distance;  % d
        values.userloc(i, 1, 2) = angles.RIS_AoD.cluster1(i).elevation; % elev
        values.userloc(i, 1, 3) = angles.RIS_AoD.cluster1(i).azimuth;   % az
    end
    
    % Cluster 2 (now cluster index 2 in 2nd dimension)
    for i = 1:3
        values.userloc(i, 2, 1) = angles.RIS_AoD.cluster2(i).distance;  % d
        values.userloc(i, 2, 2) = angles.RIS_AoD.cluster2(i).elevation; % elev
        values.userloc(i, 2, 3) = angles.RIS_AoD.cluster2(i).azimuth;   % az
    end

    % clusters = {'cluster1', 'cluster2'};
    % numUsers = 3;

    % for c = 1:length(clusters)
    %     clusterName = clusters{c};
    %     for u = 1:numUsers
    %         values.userloc(u, c, 1) = angles.RIS_AoD.(clusterName)(u).distance;  % d
    %         values.userloc(u, c, 2) = angles.RIS_AoD.(clusterName)(u).elevation; % elev
    %         values.userloc(u, c, 3) = angles.RIS_AoD.(clusterName)(u).azimuth;   % az
    %     end
    % end



     disp(values.userloc(1, 1, :)); % Display first user's location in cluster 1
     disp(values.userloc(1, 2, :)); % Display first user's location in cluster 2
     disp(values.userloc(2, 1, :)); % Display second user's location in cluster 1
     disp(values.userloc(2, 2, :)); % Display second user's location in cluster 2
     disp(values.userloc(3, 1, :)); % Display third user's location in cluster 1
     disp(values.userloc(3, 2, :)); % Display third user's location in cluster 2
     disp('BS Location:');
    disp(values.RIS_loc);
    disp(values.BS_loc);
   



    

    