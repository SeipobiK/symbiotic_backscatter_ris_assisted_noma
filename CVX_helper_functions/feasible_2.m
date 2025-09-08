function [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, status] = feasible_2(para, H_n, H_f, H_n_c, H_f_c, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n)
    % Extract configuration parameters
    numClusters = 2;
    M = para.M;
    alpha_n = para.alpha_k_n;
    alpha_f = para.alpha_k_f;
    R_n_min = para.R_min_n;
    R_f_min = para.R_min_f;
    R_c_min = para.R_c_min;
    eta_k = para.eta;
    P_max = para.P_max;
    para.noise = para.noise * (para.scal)^2;

    cvx_begin quiet sdp
        cvx_solver mosek
        cvx_expert true
        cvx_precision high
        % cvx_solver_settings( ...
        % 'MSK_IPAR_PRESOLVE_USE', 1, ...   % tighter than default 1e-8
        % 'MSK_IPAR_NUM_THREADS', 1);
        % cvx_solver_settings('MSK_WRITE_DATA_FILE', 'debug_task.opf');
        cvx_solver_settings('MSK_IPAR_MIO_NUMERICAL_EMPHASIS_LEVEL', 2);

        % cvx_solver_settings( 'MSK_IPAR_LOG', 3, 'MSK_DPAR_PRESOLVE_TOL_X', 1e-20);


        
        % Variable declarations
        variable W(M, M, numClusters) Hermitian semidefinite
        variable A_n(numClusters) nonnegative
        variable B_n(numClusters) nonnegative
        variable A_f(numClusters) nonnegative
        variable B_f(numClusters) nonnegative
        variable A_c_f(numClusters) nonnegative
        variable B_c_f(numClusters) nonnegative
        variable A_c_n(numClusters) nonnegative
        variable B_c_n(numClusters) nonnegative
        variable delta_g nonnegative
        variable R_n(numClusters) nonnegative
        variable R_f(numClusters) nonnegative
        variable R_c_n(numClusters) nonnegative
        
        % Declare dual variables
        % dual variable dual_constraint1
        % dual variable dual_constraint3
        % dual variable dual_constraint5
        dual variables dual y_1
        
        minimize(delta_g)
        
        subject to
            % Power constraint
            sum_power = 0;
            for k = 1:numClusters
                sum_power = sum_power + real(trace(W(:, :, k)));
            end
            sum_power - P_max <= delta_g;
            
            for c = 1:numClusters
                % Variable bounds
                % A_n(c) + delta_g >= 1e-4;
                % B_n(c) + delta_g >= 1e-4;
                % A_f(c) + delta_g >= 1e-4;
                % B_f(c) + delta_g >= 1e-4;
                % A_c_n(c) + delta_g >= 1e-4;
                % B_c_n(c) + delta_g >= 1e-4;
                
                % A_n(c) + delta_g <= 1e+1;
                % B_n(c) + delta_g <= 1e+1;
                % A_f(c) + delta_g <= 1e+1;
                % B_f(c) + delta_g <= 1e+1;
                % A_c_n(c) + delta_g <= 1e+1;
                % B_c_n(c) + delta_g <= 1e+1;
                
                % Taylor approximations
                taylor_approx_far(c) = log2(1 + inv_pos(A_f_prev(c) * B_f_prev(c))) - ...
                    (log2(exp(1)) * inv_pos(A_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (A_f(c) - A_f_prev(c)) - ...
                    (log2(exp(1)) * inv_pos(B_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (B_f(c) - B_f_prev(c));
                
                taylor_approx_n(c) = log2(1 + inv_pos(A_n_prev(c) * B_n_prev(c))) - ...
                    (log2(exp(1)) * inv_pos(A_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (A_n(c) - A_n_prev(c)) - ...
                    (log2(exp(1)) * inv_pos(B_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (B_n(c) - B_n_prev(c));
                
                taylor_approx_backscatter_n(c) = log2(1 + 1 ./ (A_c_prev_n(c) * B_c_prev_n(c))) - ...
                    (log2(exp(1)) / (A_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (A_c_n(c) - A_c_prev_n(c)) - ...
                    (log2(exp(1)) / (B_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (B_c_n(c) - B_c_prev_n(c));
                
                % Rate constraints
                R_f(c) - taylor_approx_far(c) <= delta_g;
                R_n(c) - taylor_approx_n(c) <= delta_g;
                R_c_n(c) - taylor_approx_backscatter_n(c) <= delta_g;
                
                R_f_min - R_f(c) <= delta_g;
                R_n_min - R_n(c) <= delta_g;
                R_c_min - R_c_n(c) <= delta_g;
                
                % Inter-cluster interference calculations
                inter_cluster_interference_near = 0;
                inter_cluster_interference_far = 0;
                inter_cluster_interference_near_b = 0;
                test_int=0;
                testa_b=0;
                
                for j = 1:numClusters
                    if j ~= c
                        inter_cluster_interference_near = inter_cluster_interference_near + ...
                            real(trace(W(:, :, j) * H_n{c}' * H_n{c}));
                        test_int = test_int + trace(H_n{c}' * H_n{c});
                            
                        
                        inter_cluster_interference_far = inter_cluster_interference_far + ...
                            real(trace(W(:, :, j) * H_f{c}' * H_f{c}));
                        
                        inter_cluster_interference_near_b = inter_cluster_interference_near_b + ...
                            real(trace(W(:, :, j) * H_f_c{c}' * H_f_c{c}));

                            testa_b = testa_b + trace(H_n_c{c}' * H_n_c{c});
                    end
                end
                
                % Create separate constraints for problematic ones
                % Constraint 1

            disp(['Near Interference :', num2str(test_int)]);
            disp(['far BST Interference :', num2str(testa_b)]);

            disp(['Near channel :',num2str(trace(H_n{c}' * H_n{c}))]);
            disp(['Far channel :',num2str(trace(H_f{c}' * H_f{c}))]);
            disp(['BST Near channel :',num2str(trace(H_n_c{c}' * H_n_c{c}))]);
            disp(['BST Far channel :',num2str(trace(H_f_c{c}' * H_f_c{c}))]);
            disp(['Noise :',num2str(para.noise)]);
            disp(['alpha_f :',num2str(alpha_f)]);
            disp(['alpha_f :',num2str(alpha_n)]);
            
               inv_pos(A_n(c)) - real(trace(W(:, :, c) * H_n{c}' * H_n{c})) * alpha_n <= delta_g;
                
                % Constraint 2
                inter_cluster_interference_near + ...
                    real(trace(W(:, :, c) * H_n_c{c}' * H_n_c{c})) * eta_k + para.noise - B_n(c) <= delta_g;
                
                % Constraint 3
                inv_pos(A_f(c)) - real(trace(W(:, :, c) * H_f{c}' * H_f{c})) * alpha_f <= delta_g;
                
                % Constraint 4
                inter_cluster_interference_far + ...
                    real(trace(W(:, :, c) * H_f{c}' * H_f{c})) * alpha_n + ...
                    real(trace(W(:, :, c) * H_f_c{c}' * H_f_c{c})) * eta_k + para.noise - B_f(c) <= delta_g;
                
                % Constraint 5
                inv_pos(A_c_n(c)) - real(trace(W(:, :, c) * H_n_c{c}' * H_n_c{c})) * eta_k <= delta_g;
                
                % Constraint 6
                inter_cluster_interference_near_b + para.noise - B_c_n(c) <= delta_g;
            end

            
            
            % PSD constraints
            for k = 1:numClusters
                W(:, :, k) + delta_g * eye(M) == hermitian_semidefinite(M);
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