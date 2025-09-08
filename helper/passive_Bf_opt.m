function [V_opt, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, obj_history,obj_history_mc, converged,cvx_status] = ...
    passive_Bf_opt(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter, mc)

    K=para.K;
    numClusters=para.K;

    % Initialize
    obj_history = zeros(max_iter, 1);
    obj_history_mc = zeros(max_iter, 1);
    converged = false;    

    for m = 1:max_iter
        % Update beamforming and Taylor parameters
      [V_opt,A_n, B_n, A_f, B_f, A_cn, B_cn,obj_curr,cvx_status] = relaxed_passive(para,w_k,G_all, g_1_all,...
                g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);

      [V_max_, eig_max] = max_eigVect(V_opt);
      current_ratio = eig_max/trace(V_opt);
      obj_history_mc(m)=current_ratio;

        if ~strcmp(cvx_status, 'Solved')
            warning('Update failed at MC %d iteration %d', mc, m); 
            break;
        end

        % Update variables
        A_n_prev = A_n; B_n_prev = B_n;
        A_f_prev = A_f; B_f_prev = B_f;
        A_c_prev_n = A_cn; B_c_prev_n = B_cn;
        obj_history(m) = obj_curr;
        % obj_history_mc(m) = obj_curr;  % Store WSR for this iteration
        
        % disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_prev')]);
        % disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_prev')]);
        % disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_prev')]);
        % disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_prev')]);
        % disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_prev_n')]);
        % disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_prev_n')]);

        % Display progress
        disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_curr)]);
        % disp(['    WSR: ', num2str(WSR)]);
        if m > 1
            disp(['    Change: ', sprintf('%.10f', abs(obj_history(m) - obj_history(m-1)))]);
            disp(['    Change(Calculated): ', sprintf('%.10f', obj_history_mc(m) - obj_history_mc(m-1))]);
            disp(['Parameter :' ,num2str(current_ratio)]);
        end

        % Check convergence
        if m > 1 && abs(obj_history(m) - obj_history(m-1)) < 1e-3
            disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
            converged = true;
            obj_history = obj_history(1:m);  % Trim unused entries
            disp(cvx_status);
            continue;
        end
    end
end