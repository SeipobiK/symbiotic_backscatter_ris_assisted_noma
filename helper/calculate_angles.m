function [angles] = calculate_angles()
% Define coordinates
BS = [0,0,2];          % Base Station (x, y, z)
RIS = [5, 2, 2];   % Ze-RIS (x, y, z)
%s           % Cluster radii (5m each)

    % Cluster 2 users [User1, User2, BD]
    users_cluster1 = [8,2,2; 
                      8,1,2; 
                      9,3,0];

    
    % Cluster 2 users [User1, User2, BD]
    users_cluster2 = [8,0,2; 
                      9,0,2; 
                      7,0,2];
    
    % Initialize output structure
    angles = struct();
    
    % ===================================================================
    % 1. AoD at BS (BS -> RIS)
    % ===================================================================
    delta_BS_RIS = RIS - BS;
    angles.BS_AoD.azimuth = atan2d(delta_BS_RIS(2), delta_BS_RIS(1)); % Azimuth (φ)
    angles.BS_AoD.elevation = atan2d(delta_BS_RIS(3), sqrt(delta_BS_RIS(1)^2 + delta_BS_RIS(2)^2)); % Elevation (θ)
    angles.BS_AoD.distance = norm(delta_BS_RIS); % Distance (d)
    
    % ===================================================================
    % 2. AoA at RIS (BS -> RIS) - Reciprocal of BS AoD
    % ===================================================================
    angles.RIS_AoA.azimuth = mod(angles.BS_AoD.azimuth + 180, 360); % Reciprocal azimuth
    angles.RIS_AoA.elevation = -angles.BS_AoD.elevation; % Negated elevation
    angles.RIS_AoA.distance = angles.BS_AoD.distance;
    
    % ===================================================================
    % 3. AoD at RIS (RIS -> Users)
    % ===================================================================
    % Cluster 1 Users
    for i = 1:size(users_cluster1, 1)
        user = users_cluster1(i, :);
        delta_RIS_user = user - RIS;
        
        azimuth = atan2d(delta_RIS_user(2), delta_RIS_user(1));
        elevation = atan2d(delta_RIS_user(3), sqrt(delta_RIS_user(1)^2 + delta_RIS_user(2)^2));
        
        angles.RIS_AoD.cluster1(i).user = sprintf('User%d', i);
        angles.RIS_AoD.cluster1(i).azimuth = azimuth;
        angles.RIS_AoD.cluster1(i).elevation = elevation;
        angles.RIS_AoD.cluster1(i).distance = norm(delta_RIS_user);
    end
    
    % Cluster 2 Users
    for i = 1:size(users_cluster2, 1)
        user = users_cluster2(i, :);
        delta_RIS_user = user - RIS;
        
        azimuth = atan2d(delta_RIS_user(2), delta_RIS_user(1));
        elevation = atan2d(delta_RIS_user(3), sqrt(delta_RIS_user(1)^2 + delta_RIS_user(2)^2));
        
        angles.RIS_AoD.cluster2(i).user = sprintf('User%d', i);
        angles.RIS_AoD.cluster2(i).azimuth = azimuth;
        angles.RIS_AoD.cluster2(i).elevation = elevation;
        angles.RIS_AoD.cluster2(i).distance = norm(delta_RIS_user);
    end
end