% Simulation with Monte Carlo for BD-RIS-assisted MISO-NOMA System
% Optimized version with preallocation, vectorization, and parallelization
clear all; clc;
rng(42); % For reproducibility

% Parameters
MC_MAX = 1;  % Monte Carlo simulations
max_iter = 30;  % Maximum iterations for each Monte Carlo run

% Preallocate arrays for results
obj_montecarlo = nan(1, MC_MAX); % Store final objective values
obj_history = nan(max_iter, MC_MAX);  % Store objectives per iteration
converg_it = zeros(MC_MAX, 1);  % Store convergence iteration for each MC run
execution_times = zeros(MC_MAX, 1);  % Store execution time per MC run

% Initialize parallel pool if available
% if isempty(gcp('nocreate'))
%     parpool(min([feature('numcores'), MC_MAX]));
% end

% Main Monte Carlo loop (parallelized)
for mc = 1:MC_MAX
    tic;
    fprintf('--- Monte Carlo Iteration %d ---\n', mc);
    
    % Initialize parameters and channels
    para = para_init();
    K = 3;
    N = para.N;
    M = para.M;
    numClusters = K;
    
    % Preallocate cell arrays for channels
    H_all = cell(1, K);
    g_1_all = cell(1, K);
    g_2_all = cell(1, K);
    g_b_all = cell(1, K);
    f1_all = cell(1, K);
    f2_all = cell(1, K);
    H_n = cell(1, K); 
    H_f = cell(1, K);
    H_nc = cell(1, K); 
    H_fc = cell(1, K);
    Theta = zeros(N, N, K);  % Preallocate Theta as 3D array
    
    % Generate channels and Theta matrices
    for i = 1:K
        [BS_array, RIS_array] = generate_arrays(para);
        [H, g_1, g_2, g_b, f1, f2] = generate_channel(para, BS_array, RIS_array);
        H_all{i} = H; 
        g_1_all{i} = g_1; 
        g_2_all{i} = g_2;
        g_b_all{i} = g_b; 
        f1_all{i} = f1; 
        f2_all{i} = f2;

        % Generate random phase shifts (vectorized)
        u = exp(1i * pi * (2 * rand(N, 1)));
        Theta(:,:,i) = diag(u);
        
        % Compute composite channels
        scal = 1e+4;
        H_n{i} = g_1_all{i}' * Theta(:,:,i) * H_all{i} * scal;
        H_f{i} = g_2_all{i}' * Theta(:,:,i) * H_all{i} * scal;
        H_nc{i} = g_b_all{i}' * Theta(:,:,i) * H_all{i} * f1_all{i} * scal;
        H_fc{i} = g_b_all{i}' * Theta(:,:,i) * H_all{i} * f2_all{i} * scal;

        % Ensure channel ordering
        if norm(H_n{i}) < norm(H_f{i})
            [H_n{i}, H_f{i}] = deal(H_f{i}, H_n{i});
            [H_nc{i}, H_fc{i}] = deal(H_fc{i}, H_nc{i});
        end
    end

    % Initialize beamforming matrices (vectorized)
    P_alloc = para.P_max / numClusters;
    W_init = zeros(M, M, numClusters);
    for k = 1:numClusters
        W_init(:,:,k) = P_alloc * (H_n{k} * H_n{k}');
    end

    % Initialize Taylor approximation parameters
    A_n_prev = ones(K,1); 
    B_n_prev = ones(K,1)*1e-0;
    A_f_prev = ones(K,1)*1e-0; 
    B_f_prev = ones(K,1);
    A_c_prev_n = ones(K,1); 
    B_c_prev_n = ones(K,1)*1e-0;

    % Find feasible starting point
    obj_array = nan(max_iter, 1);
    converged = false;
    
    for m = 1:max_iter
        [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
            update(para, W_init, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);

        if ~strcmp(cvx_status, 'Solved')
            warning('Update failed at MC %d iteration %d', mc, m);
            break;
        end
        
        % Update variables
        A_n_prev = A_n; B_n_prev = B_n;
        A_f_prev = A_f; B_f_prev = B_f;
        A_c_prev_n = A_cn; B_c_prev_n = B_cn;
        W_init = W_opt;
        obj_array(m) = obj_curr;
        
        % Check convergence
        if m > 1 && abs(obj_array(m) - obj_array(m-1)) < 1e-5
            converged = true;
            converg_it(mc) = m;
            break;
        end
    end
    
    if converged
        % Calculate final rates after convergence
        R_N = zeros(K, 1);
        R_F = zeros(K, 1);
        
        for c = 1:K
            % Calculate inter-cluster interference (vectorized)
            other_clusters = setdiff(1:K, c);
            inter_cluster_near = sum(cellfun(@(j) real(trace(W_init(:,:,j) * H_n{c}' * H_n{c})), num2cell(other_clusters)));
            inter_cluster_far = sum(cellfun(@(j) real(trace(W_init(:,:,j) * H_f{c}' * H_f{c})), num2cell(other_clusters)));
            inter_cluster_near_b = sum(cellfun(@(j) real(trace(W_init(:,:,j) * H_nc{c}' * H_nc{c})), num2cell(other_clusters)));
            inter_cluster_far_b = sum(cellfun(@(j) real(trace(W_init(:,:,j) * H_fc{c}' * H_fc{c})), num2cell(other_clusters)));

            % Calculate parameters
            A_n_prev(c) = 1/(real(trace(W_init(:,:,c) * H_n{c}' * H_n{c})) * para.alpha_k_n);
            B_n_prev(c) = inter_cluster_near + real(trace(W_init(:,:,c) * H_nc{c}' * H_nc{c})) * para.eta_k + para.noise;
            
            A_f_prev(c) = 1/(real(trace(W_init(:,:,c) * H_f{c}' * H_f{c})) * para.alpha_k_f);

            B_f_prev(c) = inter_cluster_far + real(trace(W_init(:,:,c)) * H_f{c}' * H_f{c}) * para.alpha_k_n + ...
                          real(trace(W_init(:,:,c) * H_fc{c}' * H_fc{c})) * para.eta_k + para.noise;
            
            A_c_prev_n(c) = 1/(real(trace(W_init(:,:,c) * H_nc{c}' * H_nc{c})) * para.eta_k);
            B_c_prev_n(c) = max(1e-6, inter_cluster_near_b + para.noise);

            % Calculate rates
            R_N(c) = log2(1 + 1/(A_n_prev(c)*B_n_prev(c)));
            R_F(c) = log2(1 + 1/(A_f_prev(c)*B_f_prev(c)));
        end
        
        obj_montecarlo(mc) = obj_curr;
        fprintf('MC %d converged at iteration %d\n', mc, converg_it(mc));
        fprintf('Average Near User Rate: %.4f\n', mean(R_N));
        fprintf('Average Far User Rate: %.4f\n', mean(R_F));
    else
        obj_montecarlo(mc) = nan;
    end
    
    execution_times(mc) = toc;
end

% Monte Carlo Summary
valid_results = ~isnan(obj_montecarlo);
avg_obj = mean(obj_montecarlo(valid_results));
std_obj = std(obj_montecarlo(valid_results));
avg_time = mean(execution_times(valid_results));

fprintf('\n==== Monte Carlo Summary ====\n');
fprintf('Completed runs: %d/%d\n', sum(valid_results), MC_MAX);
fprintf('Average Objective Value: %.6f\n', avg_obj);
fprintf('Standard Deviation: %.6f\n', std_obj);
fprintf('Average Execution Time: %.2f seconds\n', avg_time);
fprintf('Average Convergence Iteration: %.1f\n', mean(converg_it(valid_results)));

% Save workspace variables
save('optimized_simulation_results.mat');