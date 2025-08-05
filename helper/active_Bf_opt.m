function [W_opt, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, obj_history, converged] = ...
    active_Bf_opt(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc)

    % Initialize
    obj_history = zeros(max_iter, 1);
    converged = false;

    for m = 1:max_iter
        % Update beamforming and Taylor parameters
        [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_curr, cvx_status] = ...
            update(para, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);

        if ~strcmp(cvx_status, 'Solved')
            warning('Update failed at MC %d iteration %d', mc, m);
            disp(cvx_status);
            break;
        end

        % Update variables
        A_n_prev = A_n; B_n_prev = B_n;
        A_f_prev = A_f; B_f_prev = B_f;
        A_c_prev_n = A_cn; B_c_prev_n = B_cn;
        % W_init = W_opt;
        obj_history(m) = obj_curr;

        % Display progress
        disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_curr)]);
        if m > 1
            disp(['    Change: ', sprintf('%.10f', abs(obj_history(m) - obj_history(m-1)))]);
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
end