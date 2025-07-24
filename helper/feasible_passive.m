function [V_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt, obj_prev, cvx_status, converged] = ...
    feasible_passive(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,epsln_1,V_eignemax)
    converged = false;

    % First feasibility check
   [V_opt,epsln A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = initial_pointSearch(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,epsln_1,V_eignemax)
   
    if ~strcmp(cvx_status, 'Solved')
        warning('Feasible point not found.');
        return;
    end

    % Inner optimization loop
    for m = 1:25
        [V_opt,epsln A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = initial_pointSearch(para,w_k,G_all, g_1_all,...
          g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,epsln_1,V_eignemax)
        % Update Taylor approximation parameters
        A_n_prev = A_n_opt; B_n_prev = B_n_opt;
        A_f_prev = A_f_opt; B_f_prev = B_f_opt;
        A_c_prev_n = A_c_n_opt; B_c_prev_n = B_c_n_opt;

        % Display progress (optional)
        disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_opt')]);
        disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_opt')]);
        disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_opt')]);
        disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_opt')]);
        disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_n_opt')]);
        disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_n_opt')]);
        disp(['Iteration: ', num2str(m), ', Objective Value: ', sprintf('%.10f', obj_prev)]);
        disp(['Iteration: ', num2str(m), ', Status: ', cvx_status]);

        % Check convergence
        if strcmp(cvx_status, 'Solved') && abs(obj_prev) < 1e-5
            disp('Convergence achieved.');
            disp(['Objective Value: ', num2str(obj_prev)]);
            converged = true;
            break;
        end
    end
    % disp(size(B_c_prev_n));
end