function [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_history,converged] =passive_BF_opt(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc)

    % Initialize
    obj_history = zeros(max_iter, 1);
    converged = false;

    for m = 1:max_iter
        % Update beamforming and Taylor parameters
       [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_curr,cvx_status] = Passive_BF(para,w_k,G_all, g_1_all,...
         g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n)

        if ~strcmp(cvx_status, 'Solved')
            warning('Update failed at MC %d iteration %d', mc, m);
            break;
        end

        % Update variables
        A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
        A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
        A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
        obj_history(m) = obj_curr;

        % Display progress
        disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_curr)]);
        if m > 1
            disp(['    Change: ', sprintf('%.10f', abs(obj_history(m) - obj_history(m-1)))]);
        end
        
        disp(['    Rank(W_', num2str(m), '): ', num2str(rank(V_opt))]);
        

        % Check convergence
        if m > 1 && abs(obj_history(m) - obj_history(m-1)) < 1e-5
            disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
            converged = true;
            obj_history = obj_history(1:m);  % Trim unused entries
            break;
        end
    end
end