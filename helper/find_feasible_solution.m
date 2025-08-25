function [W_opt, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, obj_history,converged,cvx_status] = ...
    find_feasible_solution(para,Theta,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc)
    K=para.K;
    numClusters=para.K;


    % Initialize
    obj_history = zeros(max_iter, 1);
    converged = false;

    H_n = cell(1, K); H_f = cell(1, K);
    H_nc = cell(1, K); H_fc = cell(1, K);
    for i = 1:numClusters
        H_n{i}  = g_1_all{i}' * Theta * G_all;
        H_f{i}  = g_2_all{i}' * Theta * G_all;
        H_nc{i} = g_b_all{i}' * Theta * G_all * f1_all{i};
        H_fc{i} = g_b_all{i}' * Theta * G_all * f2_all{i};
    end    

    for m = 1:max_iter
        % Update beamforming and Taylor parameters
        [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
            feasible(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);

        if ~strcmp(cvx_status, 'Solved')
            warning('Update failed at MC %d iteration %d', mc, m);
            break;
        end

        % Update variables
        A_n_prev = A_n; B_n_prev = B_n;
        A_f_prev = A_f; B_f_prev = B_f;
        A_c_prev_n = A_cn; B_c_prev_n = B_cn;
        % W_init = W_opt;
        obj_history(m) = obj_curr;
        % disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_prev')]);
        % disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_prev')]);
        % disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_prev')]);
        % disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_prev')]);
        % disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_prev_n')]);
        % disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_prev_n')]);

        % Display progress
        disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.14f', obj_curr)]);
        if m > 1
            disp(['    Change: ', sprintf('%.14f', abs(obj_history(m)))]);
        end

        % Check convergence
        if m > 1 && abs(obj_history(m)) < 1e-8
            disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
            converged = true;
            obj_history = obj_history(1:m);  % Trim unused entries
            disp(cvx_status);
            break;
        end
    end
end