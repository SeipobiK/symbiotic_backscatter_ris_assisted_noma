function [V_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, status] = relaxed_passive(para, w_k, G_all, g_1_all, ...
    g_2_all, g_b_all, f1_all, f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n)

    % Extract parameters
    numClusters = para.K; % Number of clusters
    N = para.N; % Number of BS antennas
    alpha_n = para.alpha_k_n; % Near user path loss factor
    alpha_f = para.alpha_k_f; % Far user path loss factor
    R_n_min = para.R_min_n; % Minimum rate for near user
    R_f_min = para.R_min_f; % Minimum rate for far user
    R_c_min = para.R_c_min; % Minimum rate for backscatter user
    eta_k = para.eta; % Backscatter coefficient
    noise_power = para.noise * (1e+4)^2; % Noise power
    
    % Precompute identity matrix
    eye_N = eye(N);
    
    % =============================================
    % PRE-COMPUTATION OF CHANNEL PRODUCTS
    % =============================================
    
    % Initialize cell arrays for channel products
    H_n_H_n = cell(numClusters, 1);  % H_n * H_n'
    H_f_H_f = cell(numClusters, 1);  % H_f * H_f'
    H_n_c_H_n_c = cell(numClusters, 1); % H_n_c * H_n_c'
    H_f_c_H_f_c = cell(numClusters, 1); % H_f_c * H_f_c'
    
    % Initialize cell arrays for interference terms
    inter_terms_near = cell(numClusters, numClusters);
    inter_terms_far = cell(numClusters, numClusters);
    inter_terms_near_b = cell(numClusters, numClusters);
    inter_terms_far_b = cell(numClusters, numClusters);
    
    % Precompute all channel products and interference terms
    for c = 1:numClusters
        % Compute channel matrices
        H_n = diag(g_1_all{c}' * eye_N) * eye_N * G_all * f1_all{c} * w_k(:, c);
        H_f = diag(g_2_all{c}' * eye_N) * eye_N * G_all * f2_all{c} * w_k(:, c);
        H_n_c = diag(g_b_all{c}' * eye_N) * eye_N * G_all * f1_all{c} * w_k(:, c);
        H_f_c = diag(g_b_all{c}' * eye_N) * eye_N * G_all * f2_all{c} * w_k(:, c);
        
        % Store channel products
        H_n_H_n{c} = H_n * H_n';
        H_f_H_f{c} = H_f * H_f';
        H_n_c_H_n_c{c} = H_n_c * H_n_c';
        H_f_c_H_f_c{c} = H_f_c * H_f_c';
        
        % Precompute interference terms for each cluster pair
        for j = 1:numClusters
            if j ~= c
                term_near = diag(g_1_all{c}' * eye_N) * eye_N * G_all * f1_all{c} * w_k(:, j);
                term_far = diag(g_2_all{c}' * eye_N) * eye_N * G_all * f2_all{c} * w_k(:, j);
                term_near_b = diag(g_b_all{c}' * eye_N) * eye_N * G_all * f1_all{c} * w_k(:, j);
                term_far_b = diag(g_b_all{c}' * eye_N) * eye_N * G_all * f2_all{c} * w_k(:, j);
                
                inter_terms_near{c,j} = term_near * term_near';
                inter_terms_far{c,j} = term_far * term_far';
                inter_terms_near_b{c,j} = term_near_b * term_near_b';
                inter_terms_far_b{c,j} = term_far_b * term_far_b';
            end
        end
    end
    
    % =============================================
    % VECTORIZED TAYLOR APPROXIMATION TERMS
    % =============================================
    
    % Precompute Taylor approximation coefficients
    log2e = log2(exp(1));
    
    % Near user terms
    taylor_coeff_A_n = log2e ./ (A_n_prev .* (1 + A_n_prev .* B_n_prev));
    taylor_coeff_B_n = log2e ./ (B_n_prev .* (1 + A_n_prev .* B_n_prev));
    taylor_const_n = log2(1 + 1./(A_n_prev .* B_n_prev)) + ...
                     taylor_coeff_A_n .* A_n_prev + ...
                     taylor_coeff_B_n .* B_n_prev;
    
    % Far user terms
    taylor_coeff_A_f = log2e ./ (A_f_prev .* (1 + A_f_prev .* B_f_prev));
    taylor_coeff_B_f = log2e ./ (B_f_prev .* (1 + A_f_prev .* B_f_prev));
    taylor_const_f = log2(1 + 1./(A_f_prev .* B_f_prev)) + ...
                     taylor_coeff_A_f .* A_f_prev + ...
                     taylor_coeff_B_f .* B_f_prev;
    
    % Backscatter terms
    taylor_coeff_A_c_n = log2e ./ (A_c_prev_n .* (1 + A_c_prev_n .* B_c_prev_n));
    taylor_coeff_B_c_n = log2e ./ (B_c_prev_n .* (1 + A_c_prev_n .* B_c_prev_n));
    taylor_const_c_n = log2(1 + 1./(A_c_prev_n .* B_c_prev_n)) + ...
                       taylor_coeff_A_c_n .* A_c_prev_n + ...
                       taylor_coeff_B_c_n .* B_c_prev_n;
    
    % =============================================
    % CVX OPTIMIZATION
    % =============================================
    
    cvx_begin quiet sdp
        cvx_solver mosek
        % cvx_solver_settings( ...
        %     'MSK_DPAR_INTPNT_TOL_PFEAS', 1e-14, ...
        %     'MSK_DPAR_INTPNT_TOL_DFEAS', 1e-14, ...
        %     'MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-14 ...
        % );

        % Variables
        variable V(N, N) Hermitian semidefinite
        variable A_n(numClusters) nonnegative
        variable B_n(numClusters) nonnegative
        variable A_f(numClusters) nonnegative
        variable B_f(numClusters) nonnegative
        variable A_c_n(numClusters) nonnegative
        variable B_c_n(numClusters) nonnegative
        variable R_n(numClusters) nonnegative
        variable R_f(numClusters) nonnegative
        variable R_c_n(numClusters) nonnegative

        % Objective function: Maximize weighted sum rate
        maximize(sum(para.weights_n*R_n + para.weights_f*R_f + para.weights_c*R_c_n))

        subject to
            % Vectorized Taylor approximation constraints
            R_n <= taylor_const_n - taylor_coeff_A_n .* A_n - taylor_coeff_B_n .* B_n;
            R_f <= taylor_const_f - taylor_coeff_A_f .* A_f - taylor_coeff_B_f .* B_f;
            R_c_n <= taylor_const_c_n - taylor_coeff_A_c_n .* A_c_n - taylor_coeff_B_c_n .* B_c_n;
            
            % Minimum rate constraints
            R_n >= R_n_min;
            R_f >= R_f_min;
            R_c_n >= R_c_min;

            for c = 1:numClusters
                % Sum interference terms from precomputed matrices
                inter_cluster_near = 0;
                inter_cluster_far = 0;
                inter_cluster_near_b = 0;
                inter_cluster_far_b = 0;
                
                for j = 1:numClusters
                    if j ~= c
                        inter_cluster_near = inter_cluster_near + real(trace(V * inter_terms_near{c,j}));
                        inter_cluster_far = inter_cluster_far + real(trace(V * inter_terms_far{c,j}));
                        inter_cluster_near_b = inter_cluster_near_b + real(trace(V * inter_terms_near_b{c,j}));
                        inter_cluster_far_b = inter_cluster_far_b + real(trace(V * inter_terms_far_b{c,j}));
                    end
                end

                % Channel quality constraints using precomputed products
                inv_pos(A_n(c)) <= real(trace(V * H_n_H_n{c})) * alpha_n;
                B_n(c) >= inter_cluster_near + real(trace(V * H_n_c_H_n_c{c})) * eta_k + noise_power;
                
                inv_pos(A_f(c)) <= real(trace(V * H_f_H_f{c})) * alpha_f;
                B_f(c) >= inter_cluster_far + real(trace(V * H_f_H_f{c})) * alpha_n + ...
                         real(trace(V * H_f_c_H_f_c{c})) * eta_k + noise_power;
                
                % Backscatter constraints
                inv_pos(A_c_n(c)) <= real(trace(V * H_n_c_H_n_c{c})) * eta_k;
                B_c_n(c) >= inter_cluster_near_b + noise_power;
            end
            
            % Unit modulus constraints for RIS elements
            for i = 1:N
                real(V(i, i)) == 1;
                imag(V(i, i)) == 0;
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
    V_opt = V;
    status = cvx_status;
end