function [H, g_1,g_2,g_b,f_1_3,f_2_3] = generate_channel(para, BS_array, RIS_array)
    % Author: Xidong Mu (base) | Modified by Kgomotjo Seipobi

    epsilon = para.rician; % Rician factor
    
    BS_loc = para.BS_loc;
    user_loc = para.user_loc;
    
    %% BS to RIS-RIS channel
    % NLOS
    H_NLOS = 1/sqrt(2) .* ( randn(para.N,para.M) + 1i*randn(para.N,para.M) );
    
    % LOS
    a_BR = steering_vector(BS_array, -BS_loc(2), -BS_loc(3));
    a_RB = steering_vector(RIS_array, BS_loc(2), BS_loc(3));
    H_LOS = a_RB*a_BR.';
    
    % pathloss
    path_loss = para.pathloss(BS_loc(1))';

    path_loss = sqrt(10.^(- path_loss/10));
    H = path_loss .* (sqrt(epsilon/(epsilon+1)) * H_LOS + sqrt(1/(epsilon+1)) * H_NLOS);


    
    
    g_NLOS = 1/sqrt(2) .* (randn(para.N,para.K) + 1i*randn(para.N,para.K));
    
    % LOS
    g_LOS = zeros(para.N,para.K);
    g = zeros(para.N,para.K);
    path_loss_g = zeros(para.K,1);  

    for k = 1:para.K
        g_LOS(:,k) = steering_vector(RIS_array, user_loc(k,2), user_loc(k,3));
        path_loss_g(k) = para.pathloss(user_loc(k,1));  % Changed to scalar for each user
        path_loss_g(k) = sqrt(10.^((-path_loss_g(k))/10));  % Fixed variable name
        g(:,k) = path_loss_g(k) .* (sqrt(epsilon/(epsilon+1)) * g_LOS(:,k) + sqrt(1/(epsilon+1)) * g_NLOS(:,k));  % Fixed indexing and parentheses
    
    end 
  
    % pathloss
    g_1= g(:,1);
    g_2= g(:,2);
    g_b= g(:,3);




    %channel from BD to user 1 and 2
    cartesian_coords = zeros(para.K, 3);
    for i = 1:para.K
        d = para.user_loc(i,1);
        theta = deg2rad(para.user_loc(i,2)); % elevation to radians
        phi = deg2rad(para.user_loc(i,3));   % azimuth to radians
    
        x = d * cos(theta) * cos(phi);
        y = d * cos(theta) * sin(phi);
        z = d * sin(theta);
        cartesian_coords(i,:) = [x, y, z];
    end
    
    % Compute distances
    dist_1_3 = norm(cartesian_coords(1,:) - cartesian_coords(3,:));
    dist_2_3 = norm(cartesian_coords(2,:) - cartesian_coords(3,:));
    path_loss_1_3 = para.pathloss(dist_1_3);
    path_loss_1_3=sqrt(10.^(-path_loss_1_3/10));
    path_loss_2_3 = para.pathloss(dist_2_3);
    path_loss_2_3=sqrt(10.^(-path_loss_2_3/10));    
    % channel gains
    f_1_3 = path_loss_1_3 * (randn + 1i*randn)/sqrt(2);
    f_2_3 = path_loss_2_3 * (randn + 1i*randn)/sqrt(2);

    end