function [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_history,obj_history_mc,converged] =passive_BF_opt(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, A_n_prev, B_n_prev, A_f_prev, B_f_prev, A_c_prev_n, B_c_prev_n, max_iter,mc)



    % Initialize
    obj_history = zeros(max_iter, 1);
    obj_history_mc = zeros(max_iter, para.MC_MAX);
    converged = false;
    N=para.N;

     % Solve the relaxed problem
    [V_opt_init,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_curr,cvx_status] = relaxed_passive(para,w_k,G_all, g_1_all,...
     g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);
            A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
            A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
            A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 
    % disp(cvx_status);
    disp(obj_curr);


    % % Extract V_0 and relaxation parameter_0
    step_size = zeros(max_iter,1);
    relax_parameter=zeros(max_iter,1);
    max_eigenvalue_V_opt=zeros(max_iter,1);
    max_eigVector_V_opt = zeros(N,max_iter);
    [V_max, lambda_max] = max_eigVect(V_opt_init);
    max_eigenvalue_V_opt(1)=lambda_max;
    max_eigVector_V_opt(:,1)=V_max;
    U_opt = zeros(N,N,max_iter);
    U_opt(:,:,1) = V_opt_init;
    obj_history(1) = obj_curr;
    initial_ratio = max_eigenvalue_V_opt(1)/trace(U_opt(:,:,1));
    step_size(1) = 0.5*(1 - initial_ratio); % More conservative

    relax_parameter(1)=min(1, max_eigenvalue_V_opt(1)/trace(U_opt(:,:,1)) + step_size(1));

    % disp(initial_ratio);

    disp(relax_parameter(1));

    disp(step_size(1) );
    

    % disp(step_size(1));
  
    
    for m = 2:max_iter
            %  disp(A_f_prev(1));
            %  disp(A_f_prev(2));
            % disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_prev')]);
            % disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_prev')]);
            % disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_prev')]);
            % disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_prev')]);
            % disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_prev_n')]);
            % disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_prev_n')]);
            


         [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = Passive_BF(para,w_k,G_all, g_1_all,...
            g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,max_eigVector_V_opt(:,m-1),relax_parameter(m-1));


        if strcmp(cvx_status, 'Solved')
            % Update variables
            A_n_prev = A_n_opt; B_n_prev = B_n_opt; 
            A_f_prev = A_f_opt; B_f_prev = B_f_opt; 
            A_c_prev_n = A_c_n_opt;  B_c_prev_n = B_c_n_opt; 

            obj_history(m) = obj_prev;
            obj_history_mc(m,mc)=obj_prev;

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
            
            % % 5. Enforce monotonic improvement
            % if m > 1 && current_ratio < relax_parameter(m-1) - 1e-6
            %     step_size(m) = step_size(m-1)/2;
            %     relax_parameter(m) = relax_parameter(m-1)*0.9;
            %     disp('LARGEEEEEE');
            %     break;
            % end
        

            % Display progress
            disp(['Iteration: ', num2str(m), ' | Objective: ', sprintf('%.10f', obj_history(m))]);
            if m > 1
                disp(['Change: ', sprintf('%.10f', (obj_history(m) - obj_history(m-1)))]);
            end
            disp(['    Rank(V_', num2str(m), '): ', num2str(rank(V_opt))]);

            disp(relax_parameter(m-1));


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

        relax_parameter(m) = min(1, current_ratio + 0.5*(1-current_ratio)); % 10% step

        disp(relax_parameter(m));

        % relax_parameter(m) = min(0.99999,max_eigenvalue_V_opt(m)/trace(U_opt(:,:,m)) + step_size(m));

        
        

        % Check convergence
        if m > 1 && abs(obj_history(m) - obj_history(m-1)) < 1e-4 && strcmp(cvx_status, 'Solved') && abs(1-relax_parameter(m)) <= 1e-5
            disp(['MC ', num2str(mc), ' converged at iteration ', num2str(m)]);
            converged = true;
            obj_history = obj_history(1:m); 
            % relax_parameter = relax_parameter(1:m); % Trim unused entries
            continue;
        end
    end
    % disp(obj_history);
    % disp(relax_parameter);
    % fprintf('%.10f',obj_history');
    % % fprintf('%.10f',relax_parameter');
    % disp(cvx_status);
    % disp(obj_prev);
    % fprintf('%.14f',relax_parameter(end-3));
    % disp(' + ');
    % fprintf('%.10f',relax_parameter);
end