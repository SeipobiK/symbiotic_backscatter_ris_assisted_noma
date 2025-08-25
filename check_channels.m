clear all; clc;
addpath(genpath('/home/morolong/Documents/Msc/symbiotic_backscatter_ris_assisted_noma'));
rng('shuffle');  % different random numbers on each run
para=para_init();
[BS_array, RIS_array] = generate_arrays(para);
[H,g,f] = generate_channel(para, BS_array, RIS_array);
max_inner_iters=30;
MC_MAX=1;
max_iter=50;
V_max=0;


tic
% obj_inner = zeros(10);
W_test = cell(1,3);
v_test= cell(1,3);
w_k = zeros(para.M, para.K);
Theta_2 = zeros(para.N, para.N);
Theta=zeros(para.N, para.N);
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
   
        % Innerlayer iteration
        for tau_1 = 1:5
                    for i = 1:numClusters
                        g_1_all{i} = g(:, i, 1); 
                        g_2_all{i} = g(:, i, 2);
                        g_b_all{i} = g(:, i, 3); 
                        f1_all{i} = f(i,1);
                        f2_all{i} = f(i,2);
                        G_all=H*(para.scal); % Scale the channel matrix
                        H_all=G_all;

                        % H_n{i}  = g_1_all{i}' * Theta_1 * G_all;
                        % H_f{i}  = g_2_all{i}' * Theta_1* G_all;
                        % H_nc{i} = g_b_all{i}' * Theta_1* G_all * f1_all{i};
                        % H_fc{i} = g_b_all{i}' * Theta_1 * G_all * f2_all{i};
                        % disp('-------------------------');
                        % disp(['H near = ','Cluster ', num2str(tau_1),': ',num2str(norm(H_n{i})^2)]);
                        % disp(['H far = ','Cluster ', num2str(tau_1),': ',num2str(norm(H_f{i})^2)]);
                        % disp(['H near BST = ','Cluster ', num2str(tau_1),': ',num2str(norm(H_nc{i})^2)]);
                        % disp(['H far BST = ','Cluster ', num2str(tau_1),': ',num2str(norm(H_fc{i})^2)]);
  
                        disp('-------------------------');
                    end

                    A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
                    A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
                    A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;
                    V_eignemax=zeros(N, N);
                    % Algo 3: Find a feasible solution for the passive BF
                    [V_opt_feasible,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = feasible_passive(para,w_k,G_all, g_1_all,...
                    g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);

                    % Update Taylor points
                    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 



                    % Sequential constraint relaxation Algorithm for  obtaining the optimal passive BF

                    % Solve the relaxed (P2.2b)
                    [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_history,converged] =passive_BF_opt(para,w_k,G_all, g_1_all,...
                    g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc)

                    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
                    [V_mo,max_egvl]=max_eigVect(V_opt);
                    v_opt=sqrt(max_egvl)*V_mo;
                    Theta_2=diag(v_opt);
                    Theta=diag(v_opt);


                    for i = 1:numClusters
                        H_n{i}  = g_1_all{i}' * Theta_2 * G_all;
                        H_f{i}  = g_2_all{i}' * Theta_2* G_all;
                        H_nc{i} = g_b_all{i}' * Theta_2* G_all * f1_all{i};
                        H_fc{i} = g_b_all{i}' * Theta_2 * G_all * f2_all{i};
                    end
                    
                    % Initialize Taylor points
                    A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
                    A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
                    A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;


                    % %Algorithm 2: Find a feasible solution
                    % % Find feasible starting point
                    [V_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, cvx_status, converged] = ...
                    find_feasible_solution(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n,...
                    B_c_prev_n,max_inner_iters);

                    

                    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 

                    % % % Algorithm 1: Active beamforming optimization
                    % [W_opt, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, obj_history,obj_history_mc, converged,status] = ...
                    % active_Bf_opt(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc);
                    %  disp('cvx_optval :');
                    % disp(obj_history);
                    % disp('Calculated :');
                    % disp(obj_history_mc);

                    for m = 1:max_iter
                        % Update beamforming and Taylor parameters
                        [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
                            update(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);

                        for k = 1:2
                            [W_max,max_eigenvalue_w]=max_eigVect(W_opt(:,:,k));
                            w_k(:, k) = sqrt(max_eigenvalue_w)*W_max;
                            W(:,:,k) = w_k(:,k) * w_k(:,k)'; 
                        end
                        
                       [WSR, R_n, R_f, R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all,para.alpha_k_n, para.alpha_k_f, Theta);


                        if ~strcmp(cvx_status, 'Solved')
                            disp(cvx_status);
                            warning('Update failed at MC %d iteration %d', mc, m);
                            
                            break;
                        end

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
                    
                    % Recover w_k from W_opt using eigenvalue decomposition
                    W= zeros(M, M, numClusters);
                    for k = 1:K
                        [W_max,max_eigenvalue_w]=max_eigVect(W_opt(:,:,k));
                        w_k(:, k) = sqrt(max_eigenvalue_w)*W_max;
                        W(:,:,k) = w_k(:,k) * w_k(:,k)'; 
                        % disp(['w_k(:,', num2str(k), ') = ', num2str(w_k(:, k)')]);
                        % disp(norm(w_k(:,k))^2);
                        x=norm(w_k(:, k))^2;
                        fprintf('%.10f', x);
                    end
                    disp('-------------------------');
                    disp(['Total power in MC ', num2str(mc), ': ', num2str(norm(w_k(:,1))^2 + norm(w_k(:,2))^2)]);
                    disp('-------------------------');

                    disp('-------------------------');
                    disp(['Total power in MC ', num2str(mc), ': ', num2str(norm(W_opt(:,:,1)) + norm(W_opt(:,:,2)))]);
                    disp('-------------------------');

                    disp('2-------------------------');
                    disp(norm(G_all));
                    disp('-------------------------');
                    % [WSR,R_n,R_f,R_c_n] = Compute_WSR(para,w_k,G_all, g_1_all,...
                    %         g_2_all,g_b_all,f1_all,f2_all, para.alpha_k_n, para.alpha_k_n, Theta);

                    [WSR, R_n, R_f, R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all,para.alpha_k_n, para.alpha_k_f, Theta);

                    disp(['WSR: ', num2str(WSR)]);

                    disp('cvx_optval :');
                    disp(obj_history');
                    disp(size(obj_history));
                    disp('Calculated :');
                    disp(size(obj_history_mc))
                    disp(obj_history_mc');
                    disp(obj_history_mc(1));
                    

                    % Find initial point for the passive BF
                    % Initialize Taylor points

                    [WSR_, R_n, R_f, R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all,para.alpha_k_n, para.alpha_k_f, Theta_2);
                    obj_inner_(tau_1) = WSR_; 



            %         % extract the maximum eigenvalue and corresponding eigenvector
            %         [V_max_init,max_eigenvalue_v_init]=max_eigVect(V_opt_init);
            %         obj_history = zeros(max_iter, 1);
            %         U_opt = zeros(N,N,max_iter);
            %         U_opt(:,:,1) = V_opt_init;
            %         Delta = zeros(max_iter,1); 
            %         e_new = zeros(max_iter,1);
            %         max_eigenvalue=zeros(max_iter, 1);
            %         V_max_eigVect = zeros(N,max_iter);
            %         V_max_eigVect(:,1)=V_max_init;
            %         Delta(1)= (1-max_eigenvalue_v_init/trace(V_opt_init))/15;
            %         disp(Delta(1));
            %         e_new(1) = min(1, max_eigenvalue_v_init/trace(U_opt(:,:,1)) + Delta(1));
            %         disp(['e_new(1) = ', num2str(e_new(1))]);


            %         for m=2:max_iter

            %             [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = Passive_BF(para,w_k,G_all, g_1_all,...
            %                 g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,V_max_eigVect(:,m-1),e_new(m-1));

            %                 if strcmp(cvx_status, 'Solved')
            %                     % Update variables
            %                     A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
            %                     A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
            %                     A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
            %                     [V_max, lambda_max] = max_eigVect(V_opt);
            %                     U_opt(:,:,m) = V_opt;
            %                     max_eigenvalue(m)=lambda_max;
            %                     V_max_eigVect(:,m)=V_max;
            %                     obj_history(m) = obj_prev;
            %                     current_eig_ = eig(V_opt);
            %                     sorted_eig_ = sort(current_eig_, 'descend');
            %                     fprintf('%.2e\n', sorted_eig_');
            %                     disp('max_eigenvalue');
            %                     disp(e_new(m));
            %                     Delta(m)=Delta(1);
            %                     Theta = diag(V_max);   
            %                     Theta_2 = diag(V_max); %  Update Theta with the new maximum eigenvector     
            %                     for i = 1:K
            %                         g_1_all{i} = g(:, i, 1); 
            %                         g_2_all{i} = g(:, i, 2);
            %                         g_b_all{i} = g(:, i, 3); 
            %                         f1_all{i} = f(i,1);
            %                         f2_all{i} = f(i,2);
            %                         G_all=H*(para.scal); % Scale the channel matrix
            %                         H_all=G_all;
                                    
            %                         H_n{i}  = g_1_all{i}' * Theta * G_all;
            %                         H_f{i}  = g_2_all{i}' * Theta* G_all;
            %                         H_nc{i} = g_b_all{i}' * Theta* G_all * f1_all{i};
            %                         H_fc{i} = g_b_all{i}' * Theta * G_all * f2_all{i};
            %                     end
            %                     [SINR_n, SINR_f, SINR_b,WSR] = calculate_WSR(para, W, H_n, H_f, H_nc, H_fc);
            
            %                 else
            %                     A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
            %                     A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
            %                     A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt;                                 
            %                     U_opt(:,:,m) = U_opt(:,:,m-1); % Keep the previous solution
            %                     max_eigenvalue(m) = max_eigenvalue(m-1); % Keep the previous maximum eigenvalue
            %                     V_max_eigVect(:,m) = V_max_eigVect(:,m-1); % Keep the previous eigenvector
            %                     obj_history(m) = obj_history(m-1); % Keep the previous objective value
            %                     Delta(m)=Delta(1)/2;
            %                     disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_prev)]);
            %                     disp(['Iteration: ', num2str(m), ' | Calculated WSR Objective: ', sprintf('%.10f', WSR)]);

            %                     % break; % Exit the loop if CVX did not solve the problem
                                
            %                 end
            %                 disp(trace(U_opt(:,:,m)));

            %                 e_new(m) = min(1, max_eigenvalue(m)/trace(U_opt(:,:,m)) + Delta(m));
            %                 disp(e_new(m));

            %                 % Display progress
            %                 disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_prev)]);
            %                 disp(['Iteration: ', num2str(m), ' | Calculated WSR Objective: ', sprintf('%.10f', WSR)]);

            %                 if m > 1
            %                     disp(['    Change: ', sprintf('%.10f', abs(obj_history(m) - obj_history(m-1)))]);
            %                     disp(['    Relaxation parameter: ', sprintf('%.10f', e_new(m))]);
                                
            %                 end

            %                 if  abs(obj_history(m) - obj_history(m-1)) < 1e-2  && strcmp(cvx_status, 'Solved') && abs(1-e_new(m)) < 1e-4
            %                         disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
            %                         disp(['MC ', num2str(mc), ' objective ', num2str(obj_prev)]);                   
            %                         converged = true;
            %                         obj_history = obj_history(1:m);  % Trim unused entries
            %                         [V_max,max_eigenvalue]=max_eigVect(U_opt(:,:,m) );
            %                         V_max= V_max * sqrt(max_eigenvalue); % Scale the maximum eigenvector
            %                         Theta_2 = diag(V_max); 
            %                         Theta = diag(V_max);
            %                         break;
            %                 end
            %         end
            % disp('------------------------'); 


            % [WSR_, R_n, R_f, R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all,para.alpha_k_n, para.alpha_k_f, Theta_2);
            % obj_inner_(tau_1) = WSR_; 


            % disp(['WSR for tau_1 = ', num2str(tau_1), ': ', num2str(WSR_)]);        
            % disp('------------------------');
            %             [V_max,lam]=max_eigVect(V_opt);
            %              V_max=sqrt(lam) * V_max; % Scale the eigenvector
                    
            %             % Compute checks
            %             V_t   = norm(V_max)^2;
            %             V_tra = trace(V_max * V_max');
            %             V_p   = trace(V_opt_init);

            %             disp(['Cluster ', num2str(k)]);
            %             disp(['  Max Eigenvalue: ', num2str(lam)]);
            %             disp(['  Rank Estimate: ', num2str(rank(W_opt(:,:,k),1e-6))]);
            %             disp(['  V_p   = ', num2str(V_p, '%.10f')]);
            %             disp(['  V_t   = ', num2str(V_t, '%.10f')]);
            %             disp(['  V_tra = ', num2str(V_tra, '%.10f')]);
            %             disp('------------------------');	

        end 
        disp('------------------------');
        disp(obj_inner_');
        disp('------------------------');


end



