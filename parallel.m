% Simulation with Monte Carlo for BD-RIS-assisted MISO-NOMA System
% Author: Xidong Mu (base) | Modified by ChatGPT
% Date: 2025-07-02

para = para_init();
max_inner_iters=20;
addpath(genpath('/home/morolong/Documents/Msc'))
alpha_n = para.alpha_k_n; % Near user path loss factor
alpha_f = para.alpha_k_f; % Far user path loss factor
% close all; 
clear all; clc;
rng(42); % For reproducibility
MC_MAX = 1;  % Monte Carlo simulations
max_iter = 30;  % Maximum iterations for each Monte Carlo run
obj_montecarlo = zeros(1, MC_MAX); % Store objective values
obj_history = zeros(max_iter, MC_MAX);  % Store objectives per iteration
converg_it = zeros(max_iter, MC_MAX);  % Store objectives per iteration

tic
for mc = 1:MC_MAX
    disp(['--- Monte Carlo Iteration ', num2str(mc), ' ---']);
    
    para = para_init();
    K = 3;
    N = para.N;
    M = para.M;
    numClusters = 3;

    H_all = cell(1, K);
    g_1_all = cell(1, K);
    g_2_all = cell(1, K);
    g_b_all = cell(1, K);
    f1_all = cell(1, K);
    f2_all = cell(1, K);
    
    H_n = cell(1, K); H_f = cell(1, K);
    H_nc = cell(1, K); H_fc = cell(1, K);
    u = zeros(N, K);

    for i = 1:K
        [BS_array, RIS_array] = generate_arrays(para);
        [H, g_1, g_2, g_b, f1, f2] = generate_channel(para, BS_array, RIS_array);
        H_all{i} = H; g_1_all{i} = g_1; g_2_all{i} = g_2;
        g_b_all{i} = g_b; f1_all{i} = f1; f2_all{i} = f2;

        % Theta matrix
        for m = 1:N
            u(m, i) = exp(1i * pi * (2 * rand(1)));
        end
        Theta(:,:,i) = diag(u(:,i));


        % u = exp(1i*pi*(2*rand(N,1)));  % Random phase shifts
        % Theta = diag(u);  % N x N diagonal reflection matrix
        
        scal = 1e+4;
        H_n{i}  = g_1_all{i}' * Theta(:,:,i) * H_all{i} * scal;
        H_f{i}  = g_2_all{i}' * Theta(:,:,i) * H_all{i} * scal;
        H_nc{i} = g_b_all{i}' * Theta(:,:,i) * H_all{i} * f1_all{i} * scal;
        H_fc{i} = g_b_all{i}' * Theta(:,:,i) * H_all{i} * f2_all{i} * scal;

        E_h_sq = mean(abs(H_all{i}(:)).^2);
        disp(['E[|h|^2] for H_all{', num2str(i), '} = ', num2str(E_h_sq)]);

        % Ensure ordering
        if norm(H_n{i}) < norm(H_f{i})
            [H_n{i}, H_f{i}] = deal(H_f{i}, H_n{i});
        end
    end

    % Initial beamforming matrix
    W_init = zeros(M, M, numClusters);
    P_alloc = para.P_max / numClusters;
    for k = 1:numClusters
        W_init(:,:,k) = P_alloc * H_n{k} * H_n{k}';
    end
    max_inner_iters=20;

    % Initialize Taylor points
    A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
    A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
    A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;
    V_eignemax=zeros(N, N);


    %Algorithm 2: Find a feasible solution
    % Find feasible starting point
    [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, cvx_status, converged] = ...
    find_feasible_solution(para, W_init, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n,...
     B_c_prev_n,max_inner_iters);

    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
    A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 

    if ~converged
        warning('Feasible point not found in MC iteration %d.', mc);
        obj_montecarlo(mc) = NaN;
        break;
    else
        disp(['Feasible solution found with objective value: ', num2str(obj_prev)]);
    end

    % Algorithm 3: Active beamforming optimization
    [W_opt, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, obj_history_mc, converged] = ...
    active_Bf_opt(para, W_init, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc);

    if ~converged
        warning('MC %d did not converge.', mc);
    end
    obj_history(1:length(obj_history_mc), mc) = obj_history_mc;  % Store results

end
toc
