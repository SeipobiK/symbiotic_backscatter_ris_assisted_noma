clear all; clc;
addpath(genpath('/home/morolong/Documents/Msc/symbiotic_backscatter_ris_assisted_noma'));
rng('shuffle');  % different random numbers on each run
para=para_init();
[BS_array, RIS_array] = generate_arrays(para);
[H,g,f] = generate_channel(para, BS_array, RIS_array);
max_feasible=20;
max_iter=5;
V_max=0;
objectives=zeros(para.outer_iter,1);
K=para.K;


tic
% obj_inner = zeros(10);
W_test = cell(1,3);
v_test= cell(1,3);
w_k = zeros(para.M, para.K);
Theta_2 = zeros(para.N, para.N);
Theta=zeros(para.N, para.N);
% obj_history_all = zeros(max_iter, para.MC_MAX);
obj_history_activeBF=zeros(max_iter, 1);
obj_history=zeros(max_iter, 1);


rng_seeds = randi(1e6, para.MC_MAX);


Num= 2;  % Don't exceed MC runs or available cores

% if isempty(gcp('nocreate'))  % Check if pool exists
%     pool = parpool('local', Num);  % Create pool with N workers
%     fprintf('Using %d workers\n', pool.NumWorkers);
% else
%     pool = gcp;  % Get current pool
%     fprintf('Existing pool with %d workers found\n', pool.NumWorkers);
% end

for mc = 1:para.MC_MAX
    disp(['--- Monte Carlo Iteration ', num2str(mc), ' ---']);
    
    para = para_init();
    K  = para.K;
    N = para.N;
    M = para.M;
    numClusters = para.K; % Number of clusters;
    u = exp(1i * pi * (2 * rand(N, 1)));  % Single random phase shift vector (N x 1)
    Theta = diag(u); 
    Theta_1 = diag(u); 
    % norm(Theta);
    % disp(norm(Theta));
    
    H_all = cell(1, 1);
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

    x=norm(w_k(:,1))^2 +norm(w_k(:,2))^2;
    disp(x);
        
        for i = 1:numClusters
            g_1_all{i} = g(:, i, 1); 
            g_2_all{i} = g(:, i, 2);
            g_b_all{i} = g(:, i, 3); 
            f1_all{i} = f(i,1);
            f2_all{i} = f(i,2);
            G_all=H*(para.scal); % Scale the channel matrix
            H_all=G_all;
        end

        A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
        A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
        A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;

        % Algo 3: Find a feasible solution for the passive BF
        [V_opt_feasible,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status,converged] = feasible_passive(para,w_k,G_all, g_1_all,...
        g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n); 
        % A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
        % A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
        % A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt;       
        
        if ~converged
            break;
        end
                    
        % Solve the relaxed problem
        [V_opt_init,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_curr_p,cvx_status] = relaxed_passive(para,w_k,G_all, g_1_all,...
        g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);

        A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
        A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
        A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt;  

        % Extract V_0 and relaxation parameter_0
        step_size = zeros(para.outer_iter,1);
        relax_parameter=zeros(para.outer_iter,1);
        max_eigenvalue_V_opt=zeros(para.outer_iter,1);
        max_eigVector_V_opt = zeros(N,para.outer_iter);
        [V_max, lambda_max] = max_eigVect(V_opt_init);
        max_eigenvalue_V_opt(1)=lambda_max;
        max_eigVector_V_opt(:,1)=V_max;
        U_opt = zeros(N,N,para.outer_iter);
        U_opt(:,:,1) = V_opt_init;
        obj_history(1) = obj_curr_p;
        initial_ratio = max_eigenvalue_V_opt(1)/trace(U_opt(:,:,1));
        step_size(1) = 0.1*(1 - initial_ratio); % More conservative
        relax_parameter(1)=min(1, max_eigenvalue_V_opt(1)/trace(U_opt(:,:,1)) + step_size(1));
        disp(relax_parameter(1));

        A_n_prev_p = ones(K,1); B_n_prev_p = ones(K,1)*1e-0;
        A_f_prev_p = ones(K,1)*1e-0; B_f_prev_p = ones(K,1);
        A_c_prev_n_p = ones(K,1); B_c_prev_n_p = ones(K,1)*1e-0;

        [~,A_n_opt_p, B_n_opt_p, A_f_opt_p, B_f_opt_p, A_c_n_opt_p, B_c_n_opt_p,~,~,converged_p] = ...
            find_feasible_solution(para,Theta,G_all, g_1_all,...
            g_2_all,g_b_all,f1_all,f2_all, A_n_prev_p, B_n_prev_p, A_f_prev_p, B_f_prev_p, A_c_prev_n_p, B_c_prev_n_p, max_feasible, mc);

        if ~converged_p
            break;
        end

        A_n_prev_p = A_n_opt_p; B_n_prev_p = B_n_opt_p; 
        A_f_prev_p = A_f_opt_p; B_f_prev_p = B_f_opt_p; 
        A_c_prev_n_p = A_c_n_opt_p;  B_c_prev_n_p = B_c_n_opt_p; 

        % for i = 1:numClusters
        %     H_n{i}  = g_1_all{i}' * Theta * G_all;
        %     H_f{i}  = g_2_all{i}' * Theta * G_all;
        %     H_nc{i} = g_b_all{i}' * Theta * G_all * f1_all{i};
        %     H_fc{i} = g_b_all{i}' * Theta * G_all * f2_all{i};
        % end

        % % % Update beamforming and Taylor parameters
        % [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
        %     update(para, H_n, H_f, H_nc, H_fc, A_n_prev_p, B_n_prev_p, A_f_prev_p, B_f_prev_p, A_c_prev_n_p, B_c_prev_n_p);

        % A_n_prev_p = A_n; B_n_prev_p = B_n; 
        % A_f_prev_p = A_f; B_f_prev_p = B_f;
        % A_c_prev_n_p = A_cn;  B_c_prev_n_p = B_cn;
        % obj_history_activeBF(1) = obj_curr;
        % for k = 1:2
        %     [W_max,max_eigenvalue_w]=max_eigVect(W_opt(:,:,k));
        %     w_k(:, k) = sqrt(max_eigenvalue_w)*W_max;
        %     % W(:,:,k) = w_k(:,k) * w_k(:,k)'; 
        % end   

        

        %  [WSR, R_n, R_f, R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, para.alpha_k_n, para.alpha_k_f, Theta);
        %  objectives(1)=WSR;




        % Innerlayer iteration
        for m = 2:para.outer_iter

                for i = 1:numClusters
                    g_1_all{i} = g(:, i, 1); 
                    g_2_all{i} = g(:, i, 2);
                    g_b_all{i} = g(:, i, 3); 
                    f1_all{i} = f(i,1);
                    f2_all{i} = f(i,2);
                    G_all=H*(para.scal); % Scale the channel matrix
                    H_all=G_all;
                end

                % ===================Passive BF OPT==================================
                    % Update beamforming and Taylor parameters
                    % Algo 3: Find a feasible solution for the passive BF 

                    [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = Passive_BF(para,w_k,G_all, g_1_all,...
                        g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,max_eigVector_V_opt(:,m-1),relax_parameter(m-1));

                    if strcmp(cvx_status, 'Solved')
                        % Update variables
                        A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                        A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                        A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
                        obj_history(m) = obj_prev;
                        % Extract Max eigen Value and Max eigen Vector 
                        [V_max_, eig_max] = max_eigVect(V_opt);
                        U_opt(:,:,m)=V_opt;
                        max_eigenvalue_V_opt(m)=eig_max;
                        max_eigVector_V_opt(:,m)=V_max_;
                        current_eig_ = eig(V_opt);
                        sorted_eig_ = sort(current_eig_, 'descend');
                        fprintf('%.2e\n', sorted_eig_'); 
                        step_size(m)=step_size(1);
                        % % 4. Calculate current ratio
                        current_ratio = eig_max/trace(V_opt);
                        
                        % 5. Enforce monotonic improvement
                        if m > 1 && current_ratio < relax_parameter(m-1) - 1e-6
                            step_size(m) = step_size(m-1)/2;
                            relax_parameter(m) = relax_parameter(m-1)*0.9;
                            disp('LARGEEEEEE');
                            break;
                        end
                    

                        % Display progress
                        disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_history(m))]);
                        if m > 1
                            disp(['    Change: ', sprintf('%.10f', (obj_history(m) - obj_history(m-1)))]);
                        end
                        disp(['    Rank(V_', num2str(m), '): ', num2str(rank(V_opt))]);

                        % disp(1-relax_parasmeter(m-1));


                    else
                        U_opt(:,:,m) = U_opt(:,:,m-1);
                        max_eigenvalue_V_opt(m) = max_eigenvalue_V_opt(m-1);
                        obj_history(m) = obj_history(m-1);
                        step_size(m)=step_size(m-1)/2;
                        disp(step_size(m));
                        obj_history = obj_history(1:m);
                        % relax_parameter = relax_parameter(1:m);
                        if step_size(m)<1e-3
                            break;
                        end
                        
                    end
                    % disp(max_eigenvalue_V_opt(1)/trace(U_opt(:,:,m)));
                    current_ratio = max_eigenvalue_V_opt(m)/trace(U_opt(:,:,m));
                    relax_parameter(m) = min(0.9999, current_ratio + 0.5*(1-current_ratio)); % 10% step

                    % relax_parameter(m) = min(0.99999,max_eigenvalue_V_opt(m)/trace(U_opt(:,:,m)) + step_size(m));

                    
                    

                    % % Check convergence
                    % if m > 1 && abs(obj_history(m) - obj_history(m-1)) < 1e-3 && strcmp(cvx_status, 'Solved') && abs(relax_parameter(m)) >= 0.99999
                    %     disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
                    %     converged = true;
                    %     obj_history = obj_history(1:m); 
                    %     % relax_parameter = relax_parameter(1:m); % Trim unused entries
                    %     break;
                    % end
                    
                                
                % ===========================END Passive BF OPT=====================================================4


                % extract Theta
                [V_mo,max_egvl]=max_eigVect(V_opt);
                v_opt=sqrt(max_egvl)*V_mo;
                Theta=diag(v_opt);

                % % % =============================Active BF=========================================================
                %     % Update beamforming and Taylor parameters
                %              for i = 1:numClusters
                %                 H_n{i}  = g_1_all{i}' * Theta * G_all;
                %                 H_f{i}  = g_2_all{i}' * Theta * G_all;
                %                 H_nc{i} = g_b_all{i}' * Theta * G_all * f1_all{i};
                %                 H_fc{i} = g_b_all{i}' * Theta * G_all * f2_all{i};
                %              end

                %             [~,A_n_opt_p, B_n_opt_p, A_f_opt_p, B_f_opt_p, A_c_n_opt_p, B_c_n_opt_p,~,~,converged_p] = ...
                %                 find_feasible_solution(para,Theta,G_all, g_1_all,...
                %                 g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_feasible, mc);

                %             if ~converged_p
                %                 break;
                %             end

                %             A_n_prev_p = A_n_opt_p; B_n_prev_p = B_n_opt_p; 
                %             A_f_prev_p = A_f_opt_p; B_f_prev_p = B_f_opt_p; 
                %             A_c_prev_n_p = A_c_n_opt_p;  B_c_prev_n_p = B_c_n_opt_p; 
                %             % Update beamforming and Taylor parameters
                %             [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
                %                 update(para, H_n, H_f, H_nc, H_fc, A_n_prev_p, B_n_prev_p, A_f_prev_p, B_f_prev_p, A_c_prev_n_p, B_c_prev_n_p);



                %             for k = 1:2
                %                 [W_max,max_eigenvalue_w]=max_eigVect(W_opt(:,:,k));
                %                 w_k(:, k) = sqrt(max_eigenvalue_w)*W_max;
                %                 % W(:,:,k) = w_k(:,k) * w_k(:,k)'; 
                %             end
                            
                %             [WSR, R_n, R_f, R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, para.alpha_k_n, para.alpha_k_f, Theta);


                %             if ~strcmp(cvx_status, 'Solved')
                %                 disp(cvx_status);
                %                 disp('-------------------------');
                %                 testt=trace(H_n{2}'*H_n{2});
                %                 testt1=trace(H_nc{2}'*H_nc{2});
                %                 disp(testt);
                %                 disp('-------------------------');
                %                 disp(testt1);
                %                 disp('-------------------------');
                %                 warning('Update failed at MC %d iteration %d', mc, m);
                %                 break;
                %             end

                %             A_n_prev_p = A_n; B_n_prev_p = B_n; 
                %             A_f_prev_p = A_f; B_f_prev_p = B_f;
                %             A_c_prev_n_p = A_cn;  B_c_prev_n_p = B_cn;
                %             obj_history_activeBF(1) = obj_curr;

                            
                %             % disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_prev_p')]);
                %             % disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_prev_p')]);
                %             % disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_prev_p')]);
                %             % disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_prev_p')]);
                %             % disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_prev_n_p')]);
                %             % disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_prev_n_p')]);

                %             % Display progress
                %             disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_curr)]);
                %             disp(['    WSR: ', num2str(WSR)]);
                %             % if m > 1
                %             %     disp(['    Change: ', sprintf('%.10f', abs(obj_history_activeBF(m) - obj_history_activeBF(m-1)))]);
                %             %     disp(['    Change(Calculated): ', sprintf('%.10f', obj_history_activeBF(m) - obj_history_activeBF(m-1))]);
                %             % end
                %             for k = 1:size(W_opt, 3)
                %                 disp(['    Rank(W_', num2str(k), '): ', num2str(rank(W_opt(:,:,k)))]);
                %                         current_eig = eig(W_opt(:,:,k));
                %                         sorted_eig = sort(current_eig, 'descend');  
                %                         disp(['    Eigenvalues of W_', num2str(k), ': ', num2str(sorted_eig')]);
                %             end

                            % Check convergence
                            % if m > 1 && abs(obj_history_activeBF(m) - obj_history_activeBF(m-1)) < 1e-3
                            %     disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
                            %     converged = true;
                            %     obj_history = obj_history(1:m);  % Trim unused entries
                            %     disp(cvx_status);
                            %     break;
                            % end
                % ====================End ACtive BF====================================


                % ====================Calculate SR====================================
                [WSR, R_n, R_f, R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, para.alpha_k_n, para.alpha_k_f, Theta);
                objectives(m)=WSR;
        end 
end
disp('fiest   :');
disp(objectives(1));
disp(objectives');

disp('-------------------------------');
disp(para.noise * (1e+4)^2);
disp('-------------------------------');

    for i = 1:numClusters
        H_n{i}  = g_1_all{i}' * Theta * G_all;
        H_f{i}  = g_2_all{i}' * Theta * G_all;
        H_nc{i} = g_b_all{i}' * Theta * G_all * f1_all{i};
        H_fc{i} = g_b_all{i}' * Theta * G_all * f2_all{i};
         disp('startt');
        disp('============================');
        disp(trace(H_n{i}'*H_n{i}));
        disp(norm(H_f{i})^2);
        disp(norm(H_nc{i})^2);
        disp(trace(H_fc{i}'*H_fc{i}));
        disp('============================');
    end    


toc