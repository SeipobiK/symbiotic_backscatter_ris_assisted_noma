function [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = feasible(para,H_n, H_f, H_n_c, H_f_c,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n)
    % Extract configuration parameters

    numClusters = 2; % Number of clusters
    M = para.M; % Number of BS antennas
    alpha_n = para.alpha_k_n; % Near user path loss factor
    alpha_f = para.alpha_k_f; % Far user path loss factor
    R_n_min = para.R_min_n; % Minimum rate for near user
    R_f_min = para.R_min_f; % Minimum rate for far user
    R_c_min = para.R_c_min; % Minimum rate for backscatter user
    eta_k = para.eta; % Backscatter coefficient
    P_max = para.P_max; % Maximum transmit power
    para.noise = para.noise * (para.scal)^2;  % Noise scales with power
    para.P_max = para.P_max;
    % Initialize CVX variables
    cvx_begin quiet sdp
        cvx_solver mosek
        % cvx_precision high
        % cvx_precision medium
        % cvx_expert true
        cvx_solver_settings( ...
        'MSK_DPAR_INTPNT_TOL_PFEAS', 1e-10, ...
        'MSK_DPAR_INTPNT_TOL_DFEAS', 1e-10, ...
        'MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-10 ...
       );
        cvx_solver_settings('MSK_DPAR_DATA_TOL_X', 1e-10); % Tighter tolerance
        cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 60); % More time

        variable W(M,M,numClusters) Hermitian semidefinite
        variable A_n(numClusters) nonnegative % Slack variable for near users
        variable B_n(numClusters) nonnegative % Slack variable for near users
        variable A_f(numClusters) nonnegative % Slack variable for far users
        variable B_f(numClusters) nonnegative % Slack variable for far users
        variable A_c_f(numClusters) nonnegative % Slack variable for backscatter devices at far user
        variable B_c_f(numClusters) nonnegative % Slack variable for backscatter devices at far user
        variable A_c_n(numClusters) nonnegative % Slack variable for backscatter devices at near user
        variable B_c_n(numClusters) nonnegative % Slack variable for backscatter devices at near user
        variable delta_g nonnegative % Slack variable for backscatter devices at near user
        variable R_n(numClusters)  nonnegative% Slack variable for backscatter devices at near user
        variable R_f(numClusters)  nonnegative% Slack variable for backscatter devices at near user
        variable R_c_n(numClusters)  nonnegative% Slack variable for backscatter devices at near user
    
        minimise(delta_g);

        subject to
    
                sum_power = 0;
                for k = 1:numClusters
                    sum_power = sum_power + real(trace(W(:,:,k)));
                end
                sum_power-P_max <=delta_g;

            for c = 1:numClusters      
                % delta_g>=1e-6-A_n(c);
                % delta_g>=1e-6-B_n(c);
                % delta_g>=1e-6-A_f(c);
                % delta_g>=1e-6-B_f(c);
                % delta_g>=1e-6-A_c_n(c);
                % delta_g>=1e-6-B_c_n(c);

                % Ensure that the slack variables are positive
                A_n(c) + delta_g >=1e-4; 
                B_n(c) + delta_g >= 1e-4;
                A_f(c) + delta_g>= 1e-4;
                B_f(c) + delta_g>= 1e-4;
                A_c_n(c) + delta_g>= 1e-4;
                B_c_n(c) + delta_g>= 1e-4;  
                
                A_n(c) + delta_g <=1e+3; 
                B_n(c) + delta_g <= 1e+3;
                A_f(c) + delta_g <= 1e+3;
                B_f(c) + delta_g <= 1e+3;
                A_c_n(c) + delta_g <= 1e+3;
                B_c_n(c) + delta_g <= 1e+3;  

                taylor_approx_far(c) = log2(1 + inv_pos(A_f_prev(c) * B_f_prev(c))) -  ...
                                    (log2(exp(1)) * inv_pos(A_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (A_f(c) - A_f_prev(c)) - ...
                                    (log2(exp(1)) * inv_pos(B_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (B_f(c) - B_f_prev(c));

                    
                taylor_approx_n(c)= log2(1 + inv_pos(A_n_prev(c) * B_n_prev(c))) -  ...
                                   (log2(exp(1)) * inv_pos(A_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (A_n(c) - A_n_prev(c)) - ...
                                   (log2(exp(1)) * inv_pos(B_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (B_n(c) - B_n_prev(c));


                    % Backscatter device constraints (Taylor approximation) at near user
                taylor_approx_backscatter_n(c) = log2(1 + 1 ./ (A_c_prev_n(c) * B_c_prev_n(c))) - ...
                                                (log2(exp(1)) / (A_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (A_c_n(c) - A_c_prev_n(c)) - ...
                                                (log2(exp(1)) / (B_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (B_c_n(c) - B_c_prev_n(c));



                    R_f(c)- taylor_approx_far(c)<=delta_g;
                    R_n(c) - taylor_approx_n(c) <=delta_g;
                    R_c_n(c) - taylor_approx_backscatter_n(c)<=delta_g;

                    R_f_min- R_f(c)<=delta_g;
                    R_n_min - R_n(c) <= delta_g;
                    R_c_min - R_c_n(c) <= delta_g;

                    
                    inter_cluster_interference_near = 0;
                    inter_cluster_interference_far = 0;
                    inter_cluster_interference_near_b=0;
                    for j = 1:numClusters
                        if j ~= c
                            % Near user inter cluster interference
                            inter_cluster_interference_near = inter_cluster_interference_near + ...
                                real(trace(W(:,:,j)  * H_n{c}' * H_n{c}));

                            inter_cluster_interference_far= inter_cluster_interference_far + ...
                                real(trace(W(:,:,j) * H_f{c}' * H_f{c})); 

                            inter_cluster_interference_near_b= inter_cluster_interference_near_b + ...
                                real(trace(W(:,:,j)  * H_n_c{c}' * H_n_c{c})); 

                            
                        end

                    end

                        % Define slack variables based on cascaded channel

                       inv_pos(A_n(c))- real(trace(W(:,:,c)  * H_n{c}' * H_n{c})) * alpha_n<=delta_g; % Near user

                        % inter cluster interference  + backscatter interference + noise power
                       delta_g>=inter_cluster_interference_near + ...
                                real(trace(W(:,:,c) * H_n_c{c}' * H_n_c{c})) * eta_k + para.noise - B_n(c);

                        
                        inv_pos(A_f(c)) - real(trace(W(:,:,c) * H_f{c}' * H_f{c})) * alpha_f <=delta_g; 


                       delta_g>= inter_cluster_interference_far + ...
                                real(trace(W(:,:,c)  * H_f{c}' * H_f{c}))  * alpha_n + ...
                                real(trace(W(:,:,c)  * H_f_c{c}' * H_f_c{c}))   * eta_k + para.noise - B_f(c) ;


                        inv_pos(A_c_n(c)) - real(trace(W(:,:,c)  * H_n_c{c}' * H_n_c{c})) * eta_k<=delta_g;

                       delta_g>= inter_cluster_interference_near_b + para.noise - B_c_n(c) ;


            end

            % Positive semidefinite constraint
    
            for k = 1:numClusters
                W(:,:,k) + delta_g * eye(M) == hermitian_semidefinite(M);
            end
            
    cvx_end

    % Return optimized results
    obj_prev = cvx_optval; 
    A_n_opt = A_n;
    B_n_opt = B_n;
    A_f_opt = A_f;
    B_f_opt = B_f;
    A_c_n_opt = A_c_n;
    B_c_n_opt = B_c_n;
    W_opt = W;
    status = cvx_status;
end