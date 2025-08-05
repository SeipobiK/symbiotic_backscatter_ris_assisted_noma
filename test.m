clear all; clc;
addpath(genpath('/home/morolong/Documents/Msc/symbiotic_backscatter_ris_assisted_noma'));
rng('shuffle');  % different random numbers on each run
para=para_init();
[BS_array, RIS_array] = generate_arrays(para);
[H,g,f] = generate_channel(para, BS_array, RIS_array);
max_inner_iters=30;
MC_MAX=1;
max_iter=50;



tic
for mc = 1:MC_MAX
    disp(['--- Monte Carlo Iteration ', num2str(mc), ' ---']);
    
    para = para_init();
    K  = para.K;
    N = para.N;
    M = para.M;
    numClusters = para.K; % Number of clusters;
    u = exp(1i * pi * (2 * rand(N, 1)));  % Single random phase shift vector (N x 1)
    Theta = diag(u); 

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

    % AO starts here
    for m = 1:10
            for i = 1:K
                g_1_all{i} = g(:, i, 1); 
                
                g_2_all{i} = g(:, i, 2);

                g_b_all{i} = g(:, i, 3); 
                
                f1_all{i} = f(i,1);
                
                f2_all{i} = f(i,2);

                G_all=H*(1e+4); % Scale the channel matrix
                H_all=G_all;

                scal = 1e+4;
                H_n{i}  = g_1_all{i}' * Theta * G_all;
                H_f{i}  = g_2_all{i}' * Theta* G_all;
                H_nc{i} = g_b_all{i}' * Theta* G_all * f1_all{i};
                H_fc{i} = g_b_all{i}' * Theta * G_all * f2_all{i};

                E_h_sq = mean(abs(H_fc{i} (:)).^2);
                disp(['E[|h|^2] for H_n{', num2str(i), '} = ', num2str(E_h_sq)]);

                % Ensure ordering
                if norm(H_n{i}) < norm(H_f{i})
                    [H_n{i}, H_f{i}] = deal(H_f{i}, H_n{i});
                end
            end


            % Initialize Taylor points
            A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
            A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
            A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;
            V_eignemax=zeros(N, N);


            % %Algorithm 2: Find a feasible solution
            % % Find feasible starting point
            [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, cvx_status, converged] = ...
            find_feasible_solution(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n,...
            B_c_prev_n,max_inner_iters);

            A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
            A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
            A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
            disp(cvx_status);

            % obj_history(1:length(obj_history_mc), mc) = obj_history_mc;  % Store results

            w_k = zeros(M, numClusters);
            for k=1:numClusters
                current_eig = eig(W_opt(:,:,k));
                sorted_eig = sort(current_eig, 'descend');    
                % Correct dominance ratio
                sum_except_max = sum(sorted_eig(2:end));
                [V,D] = eig(W_opt(:,:,k));
                disp(max(diag(D)));

                [max_eigenvalue, index] = max(diag(D)); % The diagonal of D contains the eigenvalues.
                max_eigenvector =sqrt(max_eigenvalue) * V(:, index); % The corresponding column in V is the associated eigenvector.W
                disp(max_eigenvector);
                w_k(:,k)=max_eigenvector;
        
                
            end
            epsln_1=0;


   
        disp(['--- Iteration ', num2str(m), ' ---']);
            % % Algorithm 1: Active beamforming optimization
            [W_opt, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, obj_history_mc, converged] = ...
            active_Bf_opt(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc);
            disp(cvx_status);
                if ~converged
                    warning('MC %d did not converge.', mc);
                end

                            % Initialize Taylor points
            A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
            A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
            A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;
            V_eignemax=zeros(N, N);

            [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = feasible_passive(para,w_k,G_all, g_1_all,...
            g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,epsln_1,V_eignemax);

            A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
            A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
            A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 


            [V_opt_, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, status] = relaxed_passive(para, w_k, G_all, g_1_all, ...
            g_2_all, g_b_all, f1_all, f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);
            disp(status);

                current_eig_ = eig(V_opt_);
            sorted_eig_ = sort(current_eig_, 'descend');

            disp(sorted_eig_);

            obj_history = zeros(max_iter, 1);
            V_eignemaxx = zeros(N,max_iter);
            U_opt = zeros(N,N,max_iter);
            U_opt(:,:,1) = V_opt_;
            [V_max, lambda_max] = max_eigVect(V_opt);
            Delta = zeros(max_iter,1); 
            e_new = zeros(max_iter,1);
            Delta(1)=0.005;
            max_eigenvalue=zeros(max_iter, 1);
            V_eignemaxx(:,1)=V_max;
            e_new(1) = min(1, lambda_max/trace(U_opt(:,:,1)) + Delta(1));

        for m=2:max_iter

            [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = Passive_BF(para,w_k,G_all, g_1_all,...
                g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,V_eignemaxx(:,m-1),e_new(m-1));

                if strcmp(cvx_status, 'Solved')
                    % Update variables
                    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
                    [V_max, lambda_max] = max_eigVect(V_opt);
                    U_opt(:,:,m) = V_opt;
                    max_eigenvalue(m)=lambda_max;
                    V_eignemaxx(:,m)=V_max;
                    obj_history(m) = obj_prev;
                    current_eig_ = eig(V_opt);
                    sorted_eig_ = sort(current_eig_, 'descend');
                    fprintf('%.2e\n', sorted_eig_');
                    disp('max_eigenvalue');
                    disp(e_new(m));
                    
                end
                disp(trace(U_opt(:,:,m)));

                e_new(m) = min(1, max_eigenvalue(m)/trace(U_opt(:,:,m)) + Delta(1));
                disp(e_new(m));

                % Display progress
                disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_prev)]);
                if m > 1
                    disp(['    Change: ', sprintf('%.10f', abs(obj_history(m) - obj_history(m-1)))]);
                end

            if e_new(m)==1 && abs(obj_history(m) - obj_history(m-1)) < 1e-3 
                    disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
                    converged = true;
                    obj_history = obj_history(1:m);  % Trim unused entries
                    [V_max,max_eigenvalue]=max_eigVect(U_opt(:,:,m) );
                    Theta = diag(V_max); 


                    break;
            end
        end
            disp(e_new);
            disp((para.BD_cluster1));
            disp('Theta:');
            disp(size(Theta));
            disp(V_max);    
            % Asign  theta for next AO iteration
    end
    % end Ao

end

toc

