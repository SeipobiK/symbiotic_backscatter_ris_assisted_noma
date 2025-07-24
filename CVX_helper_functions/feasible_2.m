function [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, status] = feasible_2(para, H_n, H_f, H_n_c, H_f_c, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n)
    % Extract parameters
    numClusters = 3;
    M = para.M;
    alpha_n = para.alpha_k_n;
    alpha_f = para.alpha_k_f;
    R_n_min = para.R_min_n;
    R_f_min = para.R_min_f;
    R_c_min = para.R_c_min;
    eta_k = para.eta;
    P_max = para.P_max;
    noise = para.noise * (1e+4^2);
    
    % Precompute frequently used matrices
    HnH = cell(numClusters,1);
    HfH = cell(numClusters,1);
    HncH = cell(numClusters,1);
    HfcH = cell(numClusters,1);
    for c = 1:numClusters
        HnH{c} = H_n{c}'*H_n{c};
        HfH{c} = H_f{c}'*H_f{c};
        HncH{c} = H_n_c{c}'*H_n_c{c};
        HfcH{c} = H_f_c{c}'*H_f_c{c};
    end
    
    % Precompute Taylor approximation terms
    log2e = log2(exp(1));
    denom_n = 1 + A_n_prev.*B_n_prev;
    denom_f = 1 + A_f_prev.*B_f_prev;
    denom_c = 1 + A_c_prev_n.*B_c_prev_n;
    
    cvx_begin quiet
        cvx_solver mosek
        cvx_precision high
        cvx_solver_settings( ...
        'MSK_DPAR_INTPNT_TOL_PFEAS', 1e-14, ...
        'MSK_DPAR_INTPNT_TOL_DFEAS', 1e-14, ...
        'MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-14 ...
       );
        
        variable W(M,M,numClusters) Hermitian semidefinite
        variable A_n(numClusters) nonnegative
        variable B_n(numClusters) nonnegative
        variable A_f(numClusters) nonnegative
        variable B_f(numClusters) nonnegative
        variable A_c_n(numClusters) nonnegative
        variable B_c_n(numClusters) nonnegative
        variable delta_g nonnegative
        variable R_n(numClusters) nonnegative
        variable R_f(numClusters) nonnegative
        variable R_c_n(numClusters) nonnegative
        
        minimize(delta_g)
        
        subject to
            % Power constraint with slack
            sum_power = 0;
            for k = 1:numClusters
                sum_power = sum_power + trace(W(:,:,k));
            end
            sum_power - P_max <= delta_g;
            
            % Non-negativity constraints
            A_n >= 1e-4;
            B_n >= 1e-4;
            A_f >= 1e-4;
            B_f >= 1e-4;
            A_c_n >= 1e-4;
            B_c_n >= 1e-4;
            
            for c = 1:numClusters
                % Taylor approximation constraints with slack
                R_n(c) <= log2(1 + 1/(A_f_prev(c)*B_f_prev(c))) - ...
                         (log2e/(A_f_prev(c)*denom_f(c))) * (A_f(c) - A_f_prev(c)) - ...
                         (log2e/(B_f_prev(c)*denom_f(c))) * (B_f(c) - B_f_prev(c)) + delta_g;
                
                R_f(c) <= log2(1 + 1/(A_n_prev(c)*B_n_prev(c))) - ...
                         (log2e/(A_n_prev(c)*denom_n(c))) * (A_n(c) - A_n_prev(c)) - ...
                         (log2e/(B_n_prev(c)*denom_n(c))) * (B_n(c) - B_n_prev(c)) + delta_g;
                
                R_c_n(c) <= log2(1 + 1/(A_c_prev_n(c)*B_c_prev_n(c))) - ...
                           (log2e/(A_c_prev_n(c)*denom_c(c))) * (A_c_n(c) - A_c_prev_n(c)) - ...
                           (log2e/(B_c_prev_n(c)*denom_c(c))) * (B_c_n(c) - B_c_prev_n(c)) + delta_g;
                
                % Minimum rate constraints with slack
                R_f_min - R_f(c) <= delta_g;
                R_n_min - R_n(c) <= delta_g;
                R_c_min - R_c_n(c) <= delta_g;
                
                % Interference calculations
                inter_near = 0; inter_far = 0; inter_near_b = 0; inter_far_b = 0;
                for j = [1:c-1, c+1:numClusters]  % Skip c
                    inter_near = inter_near + trace(W(:,:,j) * HnH{c});
                    inter_far = inter_far + trace(W(:,:,j) * HfH{c});
                    inter_near_b = inter_near_b + trace(W(:,:,j) * HncH{c});
                    inter_far_b = inter_far_b + trace(W(:,:,j) * HfcH{c});
                end
                
                % Slack variable constraints with delta_g
                inv_pos(A_n(c)) - trace(W(:,:,c) * HnH{c}) * alpha_n <= delta_g;
                delta_g >= inter_near + trace(W(:,:,c) * HncH{c}) * eta_k + noise - B_n(c);
                
                inv_pos(A_f(c)) - trace(W(:,:,c) * HfH{c}) * alpha_f <= delta_g;
                delta_g >= inter_far + trace(W(:,:,c) * HfH{c}) * alpha_n + ...
                          trace(W(:,:,c) * HfcH{c}) * eta_k + noise - B_f(c);
                
                inv_pos(A_c_n(c)) - trace(W(:,:,c) * HncH{c}) * eta_k <= delta_g;
                delta_g >= inter_near_b + noise - B_c_n(c);
            end
            
            % PSD constraints with slack
            for k = 1:numClusters
                W(:,:,k) + delta_g * eye(M) == hermitian_semidefinite(M);
            end
    cvx_end
    
    % Output results
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