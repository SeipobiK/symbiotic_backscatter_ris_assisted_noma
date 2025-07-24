addpath(genpath('/home/morolong/Documents/Msc'))
para = para_init();
antenna_configs = [8];
num_configs = length(antenna_configs);
MC_MAX = 1;  % Monte Carlo runs
max_iter = 20;
numClusters = 3;
K = 3;
N = para.N;



% Storage for results (using temporary variables for parallel execution)
eig_history_cell = cell(numClusters,num_configs, MC_MAX);
dominance_ratios_cell = cell(num_configs, MC_MAX);
obj_history_cell = cell(num_configs, MC_MAX);
convergence_data_cell = cell(MC_MAX, 1);
W_opt_cell = cell(MC_MAX, num_configs);


% Pre-generate all random seeds for reproducibility
rng_seeds = randi(1e6, MC_MAX, 1);
tic;
% --- Main parallel loop ---
for mc = 1:MC_MAX
    % Set random seed for reproducibility
    rng(rng_seeds(mc));
    
    % Temporary variables for this MC run
    local_eig_history = cell(num_configs, 1);
    local_dominance_ratios = cell(num_configs, 1);
    local_obj_history = cell(num_configs, 1);
    
    % Generate arrays (RIS and BS)
    [BS_array, RIS_array] = generate_arrays(para);
    
    % Generate channels for max antenna configuration (M=8)
    H_all_max = cell(1, 1);  
    g_1_all_max = cell(1, K);  
    g_2_all_max = cell(1, K);  
    g_b_all_max = cell(1, K);  
    f1_all_max = cell(1, K);  
    f2_all_max = cell(1, K);  
    
    for i = 1:K
        [H, g_1, g_2, g_b, f1, f2] = generate_channel(para, BS_array, RIS_array);
        H_all_max = H;  
        g_1_all_max{i} = g_1;  
        g_2_all_max{i} = g_2;  
        g_b_all_max{i} = g_b;  
        f1_all_max{i} = f1;  
        f2_all_max{i} = f2;  
    end

    % Run simulations for each antenna config
    for a = 1:num_configs
        M_current = antenna_configs(a);
        local_para = para;  % Create local copy
        local_para.M = M_current;
        
        % Extract relevant channel dimensions for current M
        H_all = cell(1, 1);
        g_1_all = cell(1, K);
        g_2_all = cell(1, K);
        g_b_all = cell(1, K);
        f1_all = cell(1, K);
        f2_all = cell(1, K);
        
        for i = 1:K
            H_all = H_all_max(:, 1:M_current);
            g_1_all{i} = g_1_all_max{i};
            g_2_all{i} = g_2_all_max{i};
            g_b_all{i} = g_b_all_max{i};
            f1_all{i} = f1_all_max{i};
            f2_all{i} = f2_all_max{i};
            G_all = H_all_max * (1e+4);  % Scale the channel matrix
        end
        
        % Initialize variables
        u = zeros(N, K);
        Theta = zeros(N, N, K);
        H_n = cell(1, K); H_f = cell(1, K);
        H_nc = cell(1, K); H_fc = cell(1, K);
        
        for i = 1:K
            u(:, i) = exp(1i * pi * (2 * rand(N, 1)));
            Theta(:,:,i) = diag(u(:,i));
            
            scal = 1e+4;
            H_n{i} = g_1_all{i}' * Theta(:,:,i) * H_all * scal;
            H_f{i} = g_2_all{i}' * Theta(:,:,i) * H_all * scal;
            H_nc{i} = g_b_all{i}' * Theta(:,:,i) * H_all * f1_all{i} * scal;
            H_fc{i} = g_b_all{i}' * Theta(:,:,i) * H_all * f2_all{i} * scal;
        end
        
        % Initialize Taylor points
        A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
        A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
        A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;
        
        % Find feasible starting point
        obj_array = zeros(max_iter, 1);
        converged = false;
        
        for n = 1:1
            [W_optp, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_prev, cvx_status] = ...
                feasible_2(local_para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);
            
            if strcmp(cvx_status, 'Solved')
                for m = 1:13
                    [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, cvx_status] = ...
                        feasible_2(local_para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);
                    
                    A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
                    A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
                    A_c_prev_n = A_c_n_opt; B_c_prev_n = B_c_n_opt; 
                    
                    if ~strcmp(cvx_status, 'Solved')
                        break;
                    end
                    
                    if strcmp(cvx_status, 'Solved') && abs(obj_prev) < 1e-7
                        converged = true;
                        break;
                    end
                end
            end
        end
        
        % Main optimization loop
        for m = 1:max_iter
            [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
            update_2(local_para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);
            
            if strcmp(cvx_status, 'Solved')
                A_n_prev = A_n; B_n_prev = B_n;
                A_f_prev = A_f; B_f_prev = B_f;
                A_c_prev_n = A_cn; B_c_prev_n = B_cn;
                obj_array(m) = obj_curr;
                
                if m > 1 && abs(obj_array(m) - obj_array(m-1)) < 1e-5
                    % W_k=W_opt;

                    continue;
                end
            else
                break;
            end
        end
        
        % Store eigenvalues and dominance ratios
        temp_eig = cell(numClusters, 1);
        temp_ratios = zeros(numClusters, 1);
        W_opt_cell{mc, a} = W_opt;

        w_k = zeros(M_current, numClusters);
        for k = 1:numClusters
            current_eig = eig(W_opt(:,:,k));
            sorted_eig = sort(current_eig, 'descend');
            temp_eig{k} = sorted_eig;
            temp_ratios(k) = sorted_eig(1) / sum(sorted_eig(2:end));
            disp(sorted_eig);
            [W,D] = eig(W_opt(:,:,k));
            [max_eigenvalue, index] = max(diag(D)); % The diagonal of D contains the eigenvalues.
            max_eigenvector = sqrt(max_eigenvalue) * W(:, index); % The corresponding column in W is the associated eigenvector.
            w_k(:,k)=max_eigenvector;
            disp(max_eigenvector);
        end


        A_n_prev = ones(K,1); B_n_prev = ones(K,1)*1e-0;
        A_f_prev = ones(K,1)*1e-0; B_f_prev = ones(K,1);
        A_c_prev_n = ones(K,1); B_c_prev_n = ones(K,1)*1e-0;
        max_Veigenvector=zeros(N, N);

        [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = feasible_passive(para,w_k,G_all, g_1_all,...
        g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,epsln_1,max_Veigenvector);

        A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
        A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
        A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 

        epsln_1=0;

    for m=1:15

        [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = Passive_BF(para,w_k,G_all, g_1_all,...
          g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,epsln_1,max_Veigenvector);

        if ~strcmp(cvx_status, 'Solved')
            warning('Update failed at MC %d iteration %d', mc, m);
            disp(cvx_status);
            break;
        end

        % Update variables
        A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
        A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
        A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
        obj_history(m) = obj_prev;disp(D);

        % Display progress
        disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_prev)]);
        if m > 1
            disp(['    Change: ', sprintf('%.10f', abs(obj_history(m) - obj_history(m-1)))]);
        end

       if m > 1 && abs(obj_history(m) - obj_history(m-1)) < 1e-4
            disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
            converged = true;
            obj_history = obj_history(1:m);  % Trim unused entries
            break;
        end
    end
       
        

        [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = Passive_BF(para,w_k,G_all, g_1_all,...
          g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,max_Veigenvector,epsln_1);
    
        % --- Passive Beamforming Optimization Loop ---
% Initialize variables for passive beamforming
U_opt = zeros(N,N,max_iter);  % Stores V_opt at each iteration
obj_history = zeros(max_iter, 1);  % Objective value history
Delta = zeros(max_iter,1);  % Delta values for convergence control
e_new = zeros(50,1);  % Epsilon values for rank-1 constraint

% Initial values
Delta(1) = 0.01;  % Initial delta value
U_opt(:,:,1) = V_opt;  % Initial solution from feasible_passive
[V_max] = max_eigVect(V_opt); 
disp(trace(V_opt)); % Get dominant eigenvector
% Delta(1) = (1-min(1, V_max/trace(V_opt) + Delta(1)))/20;
% disp(Delta(1));
e_new(1) = min(1, V_max/trace(V_opt) + Delta(1));

disp(e_new(1));

% khsvclsiug

% Main optimization loop
for l = 2:40
    disp(V_max);
    % Call Passive_BF with current parameters
    [V_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, status] = ...
        Passive_BF(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, ...
                 A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, ...
                 V_max, e_new(l-1));
                 disp('whaaq');
                 disp(status);
                 disp(trace(V_opt)); 
    
    % Update variables based on optimization status
    if strcmp(status, 'Solved')

        % Store successful iteration results
        U_opt(:,:,l) = V_opt;
        obj_history(l) = obj_prev;
        
        % Update parameters for next iteration
        Delta(l) = Delta_0;
        [V_max] = max_eigVect(U_opt(:,:,l));
        e_new(l) = min(1, V_max/trace(V_opt) + Delta(l));
        
        % % Update Taylor approximation points
        A_n_prev = A_n_opt; 
        B_n_prev = B_n_opt;
        A_f_prev = A_f_opt; 
        B_f_prev = B_f_opt;
        A_c_prev_n = A_c_n_opt;  
        B_c_prev_n = B_c_n_opt;

        if m > 1
            disp(['    Change: ', sprintf('%.10f', abs(obj_history(m) - obj_history(m-1)))]);
            disp(e_new(l));
            disp(obj_prev);
            disp(['Objective: ', num2str(obj_prev), ' Parameter : ',' Iteration: ',num2str(l), ' ',num2str(e_new(l)')]);
        end
        
        % Check convergence
        if l > 1 && abs(obj_history(l) - obj_history(l-1)) < 1e-5
            disp(['Converged at iteration ', num2str(l)]);
            obj_history = obj_history(1:l);  % Trim unused entries
            U_opt = U_opt(:,:,1:l);  % Trim unused entries
            break;
        end
    else
        % If optimization failed, reduce delta and reuse previous solution
        U_opt(:,:,l) = U_opt(:,:,l-1);
        Delta(l) = Delta(l-1)/2; 
        [V_max] = max_eigVect(U_opt(:,:,l));
        e_new(l) = min(1, V_max/trace(U_opt(:,:,l)) + Delta(l));
        obj_history(l) = obj_history(l-1);
        
        % If delta becomes too small, exit
        if Delta(l) < 1e-6
            disp('Delta too small - stopping optimization');
            break;
        end
    end
end

% % Final solution is the last successful iteration
% V_opt_final = U_opt(:,:,find(obj_history ~= 0, 1, 'last'));


        disp(status);
        disp(e_new(1));
        disp(e_new(40));
        disp('Final Objective Value:');
        disp(obj_prev);
        mkvdlin

       


        







        
        local_eig_history{a} = temp_eig;
        local_dominance_ratios{a} = temp_ratios;
        local_obj_history{a} = obj_array;


        dominance_ratios_cell{mc, a} = local_dominance_ratios{a};
        obj_history_cell{mc, a} = local_obj_history{a};
    end
    

end
toc;
% --- Post-process results ---
% Initialize properly sized arrays
dominance_ratios = zeros(numClusters, MC_MAX, num_configs);
obj_history = zeros(max_iter, MC_MAX, num_configs);
dominant_vals = zeros(num_configs, numClusters, MC_MAX);
rest_sum_vals = zeros(num_configs, numClusters, MC_MAX);

% Process all MC runs
for mc = 1:MC_MAX
    for a = 1:num_configs
        % Store dominance ratios and objective history
        if ~isempty(dominance_ratios_cell{mc,a})
            dominance_ratios(:, mc, a) = dominance_ratios_cell{mc,a};
        end
        if ~isempty(obj_history_cell{mc,a})
            obj_history(:, mc, a) = obj_history_cell{mc,a};
        end
        
        % Process eigenvalues from W_opt_cell
        if ~isempty(W_opt_cell{mc,a})
            for k = 1:numClusters
                eigs_k = eig(W_opt_cell{mc,a}(:,:,k));
                eigs_k = sort(eigs_k, 'descend');
                disp(eigs_k);
                dominant_vals(a, k, mc) = eigs_k(1);
                rest_sum_vals(a, k, mc) = sum(eigs_k(2:end));
            end
        end
    end
end



% --- Clean up parallel pool ---
delete(gcp('nocreate'));

% --- Clear temporary variables ---
clear local_* temp_* H_all* g_* f_* 

% --- Post-process results ---
% Initialize properly sized arrays
dominance_ratios = zeros(numClusters, MC_MAX, num_configs);
obj_history = zeros(max_iter, MC_MAX, num_configs);
dominant_vals = zeros(num_configs, numClusters, MC_MAX);
rest_sum_vals = zeros(num_configs, numClusters, MC_MAX);

% Process all MC runs
for mc = 1:MC_MAX
    for a = 1:num_configs
        % Store dominance ratios and objective history
        if ~isempty(dominance_ratios_cell{mc,a}) && all(size(dominance_ratios_cell{mc,a}) == [numClusters, 1])
            dominance_ratios(:, mc, a) = dominance_ratios_cell{mc,a};
        end
        
        if ~isempty(obj_history_cell{mc,a}) && (length(obj_history_cell{mc,a}) == max_iter)
            obj_history(1:length(obj_history_cell{mc,a}), mc, a) = obj_history_cell{mc,a};
        end
        
        % Process eigenvalues from W_opt_cell
        if ~isempty(W_opt_cell{mc,a})
            try
                for k = 1:numClusters
                    if size(W_opt_cell{mc,a},3) >= k  % Check if cluster exists
                        eigs_k = eig(W_opt_cell{mc,a}(:,:,k));
                        eigs_k = sort(eigs_k, 'descend');
                        disp(eigs_k);
                        dominant_vals(a, k, mc) = eigs_k(1);
                        rest_sum_vals(a, k, mc) = sum(eigs_k(2:end));
                    end
                end
            catch ME
                warning('Error processing eigenvalues for mc=%d, a=%d: %s', mc, a, ME.message);
            end
        end
    end
end


delete(gcp('nocreate'));

% Calculate averages across MC runs
avg_dominant = squeeze(mean(dominant_vals, 3, 'omitnan'));
avg_rest_sum = squeeze(mean(rest_sum_vals, 3, 'omitnan'));
avg_dominance = squeeze(mean(dominance_ratios, 2, 'omitnan'));
avg_obj_history = squeeze(mean(obj_history, 2, 'omitnan'));
dateStr = datestr(now,'yyyymmdd');

% --- Plotting ---
% 1. Dominant vs Rest Eigenvalues
figure;
for k = 1:numClusters
    subplot(1, numClusters, k);
    bar_data = [avg_dominant(:,k), avg_rest_sum(:,k)];
    bar(antenna_configs, bar_data, 'grouped');
    xlabel('Number of Antennas (M)');
    ylabel('Eigenvalue Magnitude');
    title(sprintf('Cluster %d: Dominant vs Sum of Rest', k));
    legend('Dominant Eigenvalue', 'Sum of Rest Eigenvalues');
    grid on;
end


baseName_c = sprintf('results/MC%d_%s_%2f_DominantVsRestEigenvalues', MC_MAX, dateStr,para.P_max);


saveas(gcf, [baseName_c '.png']);

% 2. Dominance Ratios
figure;
bar(antenna_configs, avg_dominance', 'grouped');
xlabel('Number of Antennas (M)');
ylabel('Average Dominance Ratio');
title('Dominance Ratio vs Number of Antennas');
legend(arrayfun(@(k) ['Cluster ', num2str(k)], 1:numClusters, 'UniformOutput', false));
grid on;
saveas(gcf, 'DominanceRatios.png');

baseName_dom = sprintf('results/MC%d_%s_%2f_DominanceRatios', MC_MAX, dateStr,para.P_max);

saveas(gcf, [baseName_dom '.png']);

% 3. Convergence Behavior
figure;
hold on;
colors = {'b', 'r', 'g'};
for a = 1:num_configs
    plot(1:max_iter, avg_obj_history(:,a), 'Color', colors{a}, 'LineWidth', 2,...
        'DisplayName', [num2str(antenna_configs(a)), ' antennas'],'Visible', 'on');
end
xlim([1, max_iter]);
xlabel('Iteration');
ylabel('Weighted sum rate (bits/s/Hz)');
title('Average Convergence Behavior');
legend('show');
grid on;
% saveas(gcf, 'ConvergenceBehavior.png');
print([baseName_c '.png'],'-dpng','-r300');

baseName_conv = sprintf('results/MC%d_%s_%2f_ConvergenceBehavior', MC_MAX, dateStr,para.P_max);


saveas(gcf, [baseName_conv '.png']);

% Save all results
save(sprintf('results/_MC%d_%s.mat', MC_MAX, datestr(now, 'yyyymmdd_HHMM')), ...
    'dominance_ratios', 'avg_dominance', 'obj_history', 'avg_obj_history', ...
    'W_opt_cell', 'dominant_vals', 'rest_sum_vals', 'avg_dominant', 'avg_rest_sum');

disp('All results saved successfully.');
disp('Parallel processing completed.');
% Test plotting and saving
toc;

