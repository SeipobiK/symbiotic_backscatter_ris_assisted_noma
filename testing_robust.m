clear all; clc;
addpath(genpath('/home/morolong/Documents/Msc/symbiotic_backscatter_ris_assisted_noma'));
rng('shuffle');  % different random numbers on each run
para=para_init();
[BS_array, RIS_array] = generate_arrays(para);
[H,g,f] = generate_channel(para, BS_array, RIS_array);
max_feasible=30;
max_iter=15;
V_max=0;
objectives=zeros(para.outer_iter,1);
objectives1=zeros(para.outer_iter,1);


tic
% obj_inner = zeros(10);
W_test = cell(1,3);
v_test= cell(1,3);
w_k = zeros(para.M, para.K);
Theta_2 = zeros(para.N, para.N);
Theta=zeros(para.N, para.N);
% obj_history_all = zeros(max_iter, para.MC_MAX);
obj_history_activeBF=zeros(max_iter, para.outer_iter);
obj_history_passiveBF=zeros(max_iter, para.outer_iter);

rng_seeds = randi(1e6, para.MC_MAX);


Num= 2;  % Don't exceed MC runs or available cores

if isempty(gcp('nocreate'))  % Check if pool exists
    pool = parpool('local', Num);  % Create pool with N workers
    fprintf('Using %d workers\n', pool.NumWorkers);
else
    pool = gcp;  % Get current pool
    fprintf('Existing pool with %d workers found\n', pool.NumWorkers);
end

parfor mc = 1:para.MC_MAX
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

    %  W(0) points
    [w_mrt] = mrt_beamforming(para, H_n{1},H_n{2});
    w_k=w_mrt;
                    for i = 1:numClusters
                        g_1_all{i} = g(:, i, 1); 
                        g_2_all{i} = g(:, i, 2);
                        g_b_all{i} = g(:, i, 3); 
                        f1_all{i} = f(i,1);
                        f2_all{i} = f(i,2);
                        G_all=H*(para.scal); % Scale the channel matrix
                        H_all=G_all;
                    end



    % ============================Feasible Taylor points===============================================================
    % % %   Active BF taylor points
    %     A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
    %     A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
    %     A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;

    %     [~,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status,converged] = ...
    %         find_feasible_solution(para,Theta,G_all, g_1_all,...
    %         g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_feasible, mc);

                    
    %     if ~converged
    %          break;
    %     end

    %     % Update Taylor points
    %     A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
    %     A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
    %     A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 




    %   Passive BF taylor points
    %   Passsive BF taylor points
        % A_n_prev_p = ones(K,1); B_n_prev_p = ones(K,1)*1e-0;
        % A_f_prev_p = ones(K,1)*1e-0; B_f_prev_p = ones(K,1);
        % A_c_prev_n_p = ones(K,1); B_c_prev_n_p = ones(K,1)*1e-0;

        % [V_opt_feasible,A_n_opt_p, B_n_opt_p, A_f_opt_p, B_f_opt_p, A_c_n_opt_p, B_c_n_opt_p,obj_prev,status] = feasible_passive(para,w_k,G_all, g_1_all,...
        % g_2_all,g_b_all,f1_all,f2_all,A_n_prev_p, B_n_prev_p, A_f_prev_p, B_f_prev_p,  A_c_prev_n_p, B_c_prev_n_p);

        % % Update Taylor points
        % A_n_prev_p = A_n_opt_p; B_n_prev_p = B_n_opt_p; 
        % A_f_prev_p = A_f_opt_p; B_f_prev_p = B_f_opt_p; 
        % A_c_prev_n_p = A_c_n_opt_p;  B_c_prev_n_p = B_c_n_opt_p;

        % for m=1:30
        %     disp(['Iteration: ', num2str(m), ' A_c_prev_n_p: ', num2str(A_c_prev_n_p')]);
        %     disp(['Iteration: ', num2str(m), ' B_c_prev_n_p: ', num2str(B_c_prev_n_p')]);

        %     [V_opt_feasible,A_n_opt_p, B_n_opt_p, A_f_opt_p, B_f_opt_p, A_c_n_opt_p, B_c_n_opt_p,obj_prev,status] = relaxed_passive(para,w_k,G_all, g_1_all,...
        %         g_2_all,g_b_all,f1_all,f2_all,A_n_prev_p, B_n_prev_p, A_f_prev_p, B_f_prev_p,  A_c_prev_n_p, B_c_prev_n_p);
        %     disp(obj_prev);

        %     disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_n_opt_p')]);
        %     disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_n_opt_p')]);

        %     % disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_prev')]);
        %     % disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_prev')]);
        %     % disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_prev')]);
        %     % disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_prev')]);


        %     % Update Taylor points
        %     A_n_prev_p = A_n_opt_p; B_n_prev_p = B_n_opt_p; 
        %     A_f_prev_p = A_f_opt_p; B_f_prev_p = B_f_opt_p; 
        %     A_c_prev_n_p = A_c_n_opt_p;  B_c_prev_n_p = B_c_n_opt_p;



        % end




    % ============================Feasible Taylor points=================================================================================================


   
        % Innerlayer iteration
        for tau_1 = 1:para.outer_iter

                    for i = 1:numClusters
                        g_1_all{i} = g(:, i, 1); 
                        g_2_all{i} = g(:, i, 2);
                        g_b_all{i} = g(:, i, 3); 
                        f1_all{i} = f(i,1);
                        f2_all{i} = f(i,2);
                        G_all=H*(para.scal); % Scale the channel matrix
                        H_all=G_all;
                    end

                    % A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
                    % A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
                    % A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;

                    % % Algo 3: Find a feasible solution for the passive BF
                    % [V_opt_feasible,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status,converged] = feasible_passive(para,w_k,G_all, g_1_all,...
                    % g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);

                    % if ~converged
                    %     break;
                    % end

                    % % Update Taylor points
                    % A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    % A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    % A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 

                    % % Sequential constraint relaxation Algorithm for  obtaining the optimal passive BF
                    % % Solve the relaxed problem(P2.2b)
                    % disp(tau_1);
                    % [V_opt,~, ~, ~, ~, ~, ~,obj_history,obj_history_mc,converged] =passive_BF_opt(para,w_k,G_all, g_1_all,...
                    % g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc);
                    % % obj_history_passiveBF(:, tau_1) = obj_history(1:max_iter);
                    % % Store the objective history for this MC run
                    % % Store the results properly

                    % if ~converged
                    %     break;
                    % end


                    % [V_mo,max_egvl]=max_eigVect(V_opt);
                    % v_opt=sqrt(max_egvl)*V_mo;
                    % Theta=diag(v_opt);





                    % ====================Active BF===================================================================================================

                    disp(['------------ Inner Iteration ', num2str(tau_1), ' -------------']);
                        % %   Active BF taylor points
                    A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
                    A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
                    A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;

                    [~,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status,converged] = ...
                        find_feasible_solution(para,Theta,G_all, g_1_all,...
                        g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_feasible, mc);

                                
                    if ~converged
                        break;
                    end

                    % Update Taylor points
                    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 

                    
                    [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_op, obj_history_,obj_history_mc, converged,cvx_status] = ...
                    active_Bf_opt(para,Theta,G_all, g_1_all,...
                    g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc);
                    % obj_history_activeBF(:, tau_1) = obj_history_(1:max_iter);

                    % A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    % A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    % A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 

                    if converged
                        obj_history_all(:,mc)=obj_history_;
                    end

                    % Extract Soln
                    for k = 1:K
                        [W_max,max_eigenvalue_w]=max_eigVect(W_opt(:,:,k));
                        w_k(:, k) = sqrt(max_eigenvalue_w)*W_max;
                    end

                    % ====================Active  BF====================================================================================================



    

                %     % ====================Passive BF====================================================================================================

                %     % Sequential constraint relaxation Algorithm for  obtaining the optimal passive BF
                %     % Solve the relaxed problem(P2.2b)
                %     % disp(tau_1);
                %     A_n_prev_p = ones(K,1); B_n_prev_p = ones(K,1)*1e-0;
                %     A_f_prev_p = ones(K,1)*1e-0; B_f_prev_p = ones(K,1);
                %     A_c_prev_n_p = ones(K,1); B_c_prev_n_p = ones(K,1)*1e-0;

                %     [V_opt_feasible,A_n_opt_p, B_n_opt_p, A_f_opt_p, B_f_opt_p, A_c_n_opt_p, B_c_n_opt_p,obj_prev,status] = feasible_passive(para,w_k,G_all, g_1_all,...
                %     g_2_all,g_b_all,f1_all,f2_all,A_n_prev_p, B_n_prev_p, A_f_prev_p, B_f_prev_p,  A_c_prev_n_p, B_c_prev_n_p);

                %     % Update Taylor points
                %     A_n_prev_p = A_n_opt_p; B_n_prev_p = B_n_opt_p; 
                %     A_f_prev_p = A_f_opt_p; B_f_prev_p = B_f_opt_p; 
                %     A_c_prev_n_p = A_c_n_opt_p;  B_c_prev_n_p = B_c_n_opt_p;
                %     [V_opt,A_n_opt_p, B_n_opt_p, A_f_opt_p, B_f_opt_p, A_c_n_opt_p, B_c_n_opt_p,obj_history,obj_history_mc,converged] =passive_BF_opt(para,w_k,G_all, g_1_all,...
                %     g_2_all,g_b_all,f1_all,f2_all, A_n_prev_p, B_n_prev_p, A_f_prev_p, B_f_prev_p,  A_c_prev_n_p, B_c_prev_n_p, max_iter, mc);

                %     % Update Taylor points
                %     % A_n_prev_p = A_n_opt_p; B_n_prev_p = B_n_opt_p; 
                %     % A_f_prev_p = A_f_opt_p; B_f_prev_p = B_f_opt_p; 
                %     % A_c_prev_n_p = A_c_n_opt_p;  B_c_prev_n_p = B_c_n_opt_p;
                %     % obj_history_passiveBF(:, tau_1) = obj_history(1:max_iter);
                %     % Store the objective history for this MC run
                %     % Store the results properly


                %     % % Extract Soln
                %     [V_mo,max_egvl]=max_eigVect(V_opt);
                    
                %     v_opt=sqrt(max_egvl)*V_mo;
                %     theta = exp(1j * angle(v_opt));
                %     % Theta=diag(theta);

                % [U,S] = eig((V_opt+V_opt')/2);
                % [lambda1, idx] = max(real(diag(S)));
                % v = sqrt(lambda1)*U(:,idx);

                % cand{1} = exp( 1j*angle(v));      % +phase
                % cand{2} = exp(-1j*angle(v));      % -phase   <-- often the right one
                % cand{3} = conj(cand{1});          % conjugate
                % cand{4} = conj(cand{2});          % conjugate of -phase

                % bestWSR = -inf; bestTheta = [];
                % for t = 1:4
                %     theta = cand{t};
                %     Theta  = diag(theta);
                %     % Use ONE interface consistently (see §2 below)
                %     WSR = calculate_WSR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, para.alpha_k_n, para.alpha_k_f, Theta);
                %     if WSR > bestWSR
                %         bestWSR = WSR; bestTheta = Theta;
                %     end
                % end
                % Theta = bestTheta;
                   
                    % ====================Passive BF====================================================================================================
                    
                    

        
        
                    % [SR,R_n,R_f,R_c_n] = Compute_WSR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, para.alpha_k_n, para.alpha_k_f, v_opt);
                    % objectives1(tau_1)=SR;   
                    % [qq,R_n,R_f,R_c_n] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, para.alpha_k_n, para.alpha_k_f, Theta);
                    % % objectives1(tau_1)=SR; 
                    %  [WSR, SINR_n, SINR_f, SINR_b] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all,para.alpha_k_n, para.alpha_k_f,Theta);

                    % [WSR,R_n,R_f,R_c_n,~,~] = calculate_WSR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, para.alpha_k_n, para.alpha_k_f, Theta);
                    % objectives1(tau_1)=WSR;  
                    % disp('Correct  :')
                    % disp(rank(V_opt));
                    % disp(WSR);


                    

        end 
end
% disp('fiest   :');
% disp(objectives1);
% disp(objectives');
% disp(norm(Theta)^2);



% disp(norm(v_opt)^2);
% disp(sqrt(max_egvl));

% % disp(R_c_n);
% % disp(R_n);


% disp('-------------------------------');
% disp(para.noise * (1e+4)^2);
% disp('-------------------------------');

%     % for i = 1:numClusters
%     %     H_n{i}  = g_1_all{i}' * Theta * G_all;
%     %     H_f{i}  = g_2_all{i}' * Theta * G_all;
%     %     H_nc{i} = g_b_all{i}' * Theta * G_all * f1_all{i};
%     %     H_fc{i} = g_b_all{i}' * Theta * G_all * f2_all{i};
%     %      disp('startt');
%     %     disp('============================');
%     %     disp(trace(H_n{i}'*H_n{i}));
%     %     disp(norm(H_f{i})^2);
%     %     disp(norm(H_nc{i})^2);
%     %     disp(trace(H_fc{i}'*H_fc{i}));
%     %     disp('============================');



%     % end    


% toc
% % % % % Calculate average
avg_obj = mean(obj_history_all, 2);
disp(avg_obj);

% % % % Plot
% % figure;


plot(1:max_iter, avg_obj, 'LineWidth', 2);
xlabel('Iteration'); ylabel('Weighted Sum Rate(bits/sec)'); title('Convergence');
grid on;
xlim([1, max_iter]);
dateStr = datestr(now,'yyyymmdd');
baseName_conv = sprintf('results/MC%d_%s_%2f_ConvergenceBehavior', para.MC_MAX, dateStr,para.P_max);
saveas(gcf, [baseName_conv '.png']);

% Create base line plot
figure;
h = plot(1:max_iter, avg_obj, 'b-', 'LineWidth', 1.5);
hold on;

% Add markers at specific points
marker_points = 1:5:max_iter;  % Show markers every 5 iterations
plot(marker_points, avg_obj(marker_points), 'bo', ...
    'MarkerSize', 6, ...
    'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', 'w');




% % % Add error bars if needed
% % % Final touches
% % hold off;
% % xlabel('Iteration Number');
% % ylabel('Objective Value');
% % title('Convergence with Marked Points');
% % legend('Objective', 'Markers', 'Std Dev');
% % grid on;

% % disp(MC_MAX);



% % % Cleanup parallel pool
delete(gcp('nocreate'));

% % %% Visualization
% % % Find maximum number of iterations across all MC runs
% % max_iters = max(cellfun(@length, obj_history_all));

% % % Pad shorter runs with NaNs for uniform matrix
% % obj_matrix = nan(MC_MAX, max_iters);
% % for mc = 1:MC_MAX
% %     obj_matrix(mc, 1:length(obj_history_all{mc})) = obj_history_all{mc};
% % end

% % % Calculate statistics
% % mean_obj = mean(obj_matrix, 1, 'omitnan');
% % std_obj = std(obj_matrix, 0, 1, 'omitnan');

% % % Create figure
% % figure('Position', [100 100 1200 500], 'Color', 'w');

% % % Plot 1: All MC runs with mean
% % subplot(1,2,1);
% % hold on;
% % for mc = 1:MC_MAX
% %     plot(obj_history_all{mc}, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
% % end
% % h_mean = plot(mean_obj, 'LineWidth', 2, 'Color', [0 0.447 0.741]);
% % xlabel('Iteration', 'FontWeight', 'bold');
% % ylabel('Objective Value', 'FontWeight', 'bold');
% % title('Monte Carlo Simulations with Mean', 'FontSize', 12);
% % grid on;
% % legend(h_mean, 'Mean', 'Location', 'best');

% % % % Plot 2: Mean ± standard deviation
% % % subplot(1,2,2);
% % % shadedErrorBar(1:max_iters, mean_obj, std_obj, ...
% % %     'lineprops', {'Color', [0 0.447 0.741], 'LineWidth', 2});
% % % xlabel('Iteration', 'FontWeight', 'bold');
% % % ylabel('Objective Value', 'FontWeight', 'bold');
% % % title('Mean ± 1 Standard Deviation', 'FontSize', 12);
% % % grid on;

% %     % % Add shadedErrorBar function if needed
% %     % if ~exist('shadedErrorBar', 'file')
% %     %     warning('shadedErrorBar not found - installing from MATLAB Central');
% %     %     websave('shadedErrorBar.m', ...
% %     %         'https://raw.githubusercontent.com/raacampbell/shadedErrorBar/master/shadedErrorBar.m');
% %     %     rehash;
% %     % end
