function [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, cvx_status, converged] = ...
    find_feasible_solution(para, W_init, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_inner_iters)

    converged = false;

    % First feasibility check
    [W_opt, A_n, B_n, A_f, B_f, A_cn, B_cn, obj_prev, cvx_status] = ...
        feasible(para, W_init, H_n, H_f, H_nc, H_fc, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n);

    if ~strcmp(cvx_status, 'Solved')
        warning('Feasible point not found.');
        return;
    end

    % Inner optimization loop
    for m = 1:max_inner_iters
        [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, cvx_status] = ...
            feasible(para, W_opt, H_n, H_f, H_nc, H_fc, A_n, B_n, A_f, B_f, A_cn, B_cn);

        % Update Taylor approximation parameters
        A_n = A_n_opt; B_n = B_n_opt;
        A_f = A_f_opt; B_f = B_f_opt;
        A_cn = A_c_n_opt; B_cn = B_c_n_opt;

        % Display progress (optional)
        disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_opt')]);
        disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_opt')]);
        disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_opt')]);
        disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_opt')]);
        disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_n_opt')]);
        disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_n_opt')]);
        disp(['Iteration: ', num2str(m), ', Objective Value: ', sprintf('%.15f', obj_prev)]);
        disp(['Iteration: ', num2str(m), ', Status: ', cvx_status]);

        % Check convergence
        if strcmp(cvx_status, 'Solved') && abs(obj_prev) < 1e-13
            disp('Convergence achieved.');
            disp(['Objective Value: ', num2str(obj_prev)]);
            converged = true;
            break;
        end
    end
end