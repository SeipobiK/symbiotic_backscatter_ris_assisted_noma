clear; clc;
addpath(genpath('/home/morolong/Documents/Msc/symbiotic_backscatter_ris_assisted_noma'));
rng('shuffle');  % different random numbers on each run
para=para_init();
[BS_array, RIS_array] = generate_arrays(para);
[H,g,f] = generate_channel(para, BS_array, RIS_array);
max_inner_iters=30;
max_inner_iter=2;
MC_MAX=1;
max_iter=50;
V_max=0;
obj_history = zeros(max_iter, 1);
obj_history_mc = zeros(max_iter, 1);





tic
% obj_inner = zeros(10);
W_test = cell(1,3);
v_test= cell(1,3);
w_k_inner_iters = zeros(para.M, para.K,max_inner_iter);
Theta_2=zeros(para.N,para.N,max_inner_iters);
Theta_inner = zeros(para.N, para.N,max_inner_iter);
Theta=zeros(para.N, para.N);
obj_inner_iteration=zeros(max_inner_iter,1)

for mc = 1:MC_MAX
    
    disp(['--- Monte Carlo Iteration ', num2str(mc), ' ---']);
    
    para = para_init();
    K  = para.K;
    N = para.N;
    M = para.M;
    numClusters = para.K; % Number of clusters;
    u = exp(1i * pi * (2 * rand(N, 1)));  % Single random phase shift vector (N x 1)
    Theta_1 = diag(u); 
    % norm(Theta);
    % disp(norm(Theta));
    
    H_all = cell(1, K);
    g_1_all = cell(1, K);
    g_2_all = cell(1, K);
    g_b_all = cell(1, K);
    f1_all = cell(1, K);
    f2_all = cell(1, K);

    
    H_n = cell(1, K); H_f = cell(1, K);
    H_nc = cell(1, K); H_fc = cell(1, K);
    u = zeros(N, K);
    [BS_array, RIS_array] = generate_arrays(para);
    [H, g,f] = generate_channel(para, BS_array, RIS_array);
    g_1_all{1} = g(:, 1, 1); 
    g_1_all{2} = g(:, 2, 1);
    H_n{1}  = g_1_all{1}' * Theta_1 * H*para.scal;
    H_n{2}  = g_1_all{2}' * Theta_1 * H*para.scal;


    [w_mrt] = mrt_beamforming(para, H_n{1},H_n{2});
    w_k=w_mrt;
    

        Theta_inner(:, :,1)=Theta_1;
        w_k_inner_iters(:,:,1)=w_mrt;
        % Innerlayer iteration
        for tau_1 = 1:max_inner_iter
                    
                    for i = 2:numClusters
                        g_1_all{i} = g(:, i, 1); 
                        g_2_all{i} = g(:, i, 2);
                        g_b_all{i} = g(:, i, 3); 
                        f1_all{i} = f(i,1);
                        f2_all{i} = f(i,2);
                        G_all=H*(para.scal); % Scale the channel matrix
                        H_all=G_all;
                        H_n{i}  = g_1_all{i}' * Theta_inner(:, :,tau_1-1) * G_all;
                        H_f{i}  = g_2_all{i}' * Theta_inner(:, :,tau_1-1)* G_all;
                        H_nc{i} = g_b_all{i}' * Theta_inner(:, :,tau_1-1)* G_all * f1_all{i};
                        H_fc{i} = g_b_all{i}' * Theta_inner(:, :,tau_1-1) * G_all * f2_all{i};
                    end

                    % Initialize Taylor points
                    A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
                    A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
                    A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;


                    % %Algorithm 2: Find a feasible solution
                    % % Find feasible starting point
                    [W_feasible, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, cvx_status, converged] = ...
                    find_feasible_solution(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n,...
                    B_c_prev_n,max_inner_iters);

                    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 


                    % Algorithm 1: Active beamforming optimization
                    for m = 1:max_iter
                        % Update beamforming and Taylor parameters
                        [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
                            update(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);
                        
                        if ~strcmp(cvx_status, 'Solved')
                            disp(cvx_status);
                            warning('Update failed at MC %d iteration %d', mc, m);
                            
                            break;
                        end

                        for k = 1:para.K
                            [W_max,max_eigenvalue_w]=max_eigVect(W_opt(:,:,k));
                            w_k(:, k) = sqrt(max_eigenvalue_w)*W_max;
                            % W(:,:,k) = w_k(:,k) * w_k(:,k)'; 
                        end
                        
                       [WSR, R_n, R_f, R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all,para.alpha_k_n, para.alpha_k_f, Theta_1);

                        % Update variables
                        A_n_prev = A_n; B_n_prev = B_n;
                        A_f_prev = A_f; B_f_prev = B_f;
                        A_c_prev_n = A_cn; B_c_prev_n = B_cn;
                        % W_init = W_opt;
                        obj_history(m) = obj_curr;
                        obj_history_mc(m) = WSR;  % Store WSR for this iteration

                        % Display progress
                        disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_curr)]);
                        disp(['    WSR: ', num2str(WSR)]);
                        if m > 1
                            disp(['    Change: ', sprintf('%.10f', abs(obj_history(m) - obj_history(m-1)))]);
                            disp(['    Change(Calculated): ', sprintf('%.10f', obj_history_mc(m) - obj_history_mc(m-1))]);
                        end
                        for k = 1:size(W_opt, 3)
                            disp(['    Rank(W_', num2str(k), '): ', num2str(rank(W_opt(:,:,k)))]);
                                    current_eig = eig(W_opt(:,:,k));
                                    sorted_eig = sort(current_eig, 'descend');  
                                    disp(['    Eigenvalues of W_', num2str(k), ': ', num2str(sorted_eig')]);
                        end

                        % Check convergence
                        if m > 1 && abs(obj_history(m) - obj_history(m-1)) < 1e-3
                            disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
                            converged = true;
                            obj_history = obj_history(1:m);  % Trim unused entries
                            disp(cvx_status);
                            break;
                        end
                    end

                    % Extract the optimal beamforming vectors
                    W= zeros(M, M, numClusters);
                    for k = 1:K
                        [W_max,max_eigenvalue_w]=max_eigVect(W_opt(:,:,k));
                        w_k(:, k) = sqrt(max_eigenvalue_w)*W_max;
                        W(:,:,k) = w_k(:,k) * w_k(:,k)'; 
                        x=norm(w_k(:, k))^2;
                        fprintf('%.10f', x);
                    end
                    w_k_inner_iters(:,:,tau_1) = w_k;  % Store the beamforming vectors for this iteration
                    disp(['Iteration ', num2str(tau_1), ' completed.']);


                    % Passive BF
                    % Find taylor initial points

                    % %Algorithm 2: Find a feasible solution
                    % % Find feasible starting point
                    % % % Algo 5: Find a feasible solution for the passive BF
                    [V_feasible,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = feasible_passive(para,w_k,G_all, g_1_all,...
                    g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);

                    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 


                    % % % ===================Sequential constraint relaxation Algorithm=============================

                    % Solve Relaxed Passive BF without Rank relaxation constrain
                            % Solve the relaxed (P2.2b)
                    [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_history,converged] =passive_BF_opt(para,w_k,G_all, g_1_all,...
                    g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc);

                    % % % Extract V
                    [V_mo,max_egvl]=max_eigVect(V_opt);
                    v_opt=sqrt(max_egvl)*V_mo;
                    Theta_2(:,:,tau_1)=diag(v_opt);
                    Theta_1=diag(v_opt);

                    % calculate Sum rate
                    [WSR, R_n, R_f, R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all,para.alpha_k_n, para.alpha_k_f, Theta_1);
                    obj_inner_iteration(tau_1)=WSR;

                    

        end
        disp('%.10f',obj_inner_iteration);

end