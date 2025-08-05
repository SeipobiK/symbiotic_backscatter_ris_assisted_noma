function [V_opt,epsln,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = initial_pointSearch(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,epsln_1,V_eignemax)
   
   numClusters  = para.K;% Number of clusters
   N = para.N; % Number of BS antennas
   M = para.M; % Number of BS antennas
   alpha_n = para.alpha_k_n; % Near user path loss factor
   alpha_f = para.alpha_k_f; % Far user path loss factor
   R_n_min = para.R_min_n; % Minimum rate for near user
   R_f_min = para.R_min_f; % Minimum rate for far user
   R_c_min = para.R_c_min; % Minimum rate for backscatter user
   eta_k = para.eta; % Backscatter coefficient
   P_max = para.P_max; % Maximum transmit power
   para.noise = para.noise * (1e+4)^2;  % Noise scales with power
   para.P_max = para.P_max;
%    scal=1e+6;

   cvx_begin quiet   sdp
       % cvx_solver sedumi
       cvx_solver mosek
       cvx_precision high
       cvx_precision high
        % cvx_solver_settings( ...
        %     'MSK_DPAR_INTPNT_TOL_PFEAS', 1e-14, ...
        %     'MSK_DPAR_INTPNT_TOL_DFEAS', 1e-14, ...
        %     'MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-14 ...
        % );

       %  (59h)
       variable V(N,N) Hermitian semidefinite  
       variable A_n(numClusters) nonnegative % Slack variable for near users
       variable B_n(numClusters) nonnegative % Slack variable for near users
       variable A_f(numClusters) nonnegative % Slack variable for far users
       variable B_f(numClusters) nonnegative % Slack variable for far users
       variable A_c_f(numClusters) nonnegative % Slack variable for backscatter devices at far user
       variable B_c_f(numClusters) nonnegative % Slack variable for backscatter devices at far user
       variable A_c_n(numClusters) nonnegative % Slack variable for backscatter devices at near user
       variable B_c_n(numClusters) nonnegative % Slack variable for backscatter devices at near user
       variable R_n(numClusters)  nonnegative% Slack variable for backscatter devices at near user
       variable R_f(numClusters)  nonnegative% Slack variable for backscatter devices at near user
       variable R_c_n(numClusters)  nonnegative% Slack variable for backscatter devices at near user
    %    variable delta_p% Slack variable for backscatter devices at near user
       variable delta_p nonnegative  % This should be a variable, not an expression

       expressions taylor_approx_far(numClusters, 1) taylor_approx_n(numClusters, 1) taylor_approx_backscatter_n(numClusters, 1)
       

   


       % Objective function: Maximize weighted sum rate (59a)
       minimise delta_p

        subject to
  
           for c = 1:numClusters
               H_n_p=G_all*f1_all{c}*w_k(:, c);
               H_f_p=G_all*f2_all{c}*w_k(:, c);
               H_fc_p=G_all*f2_all{c}*w_k(:, c);
               H_nc_p=G_all*f1_all{c}*w_k(:, c);
       
           
       
               [J_t_n]=permut(H_n_p);
               [J_t_f]=permut(H_f_p);
               [J_t_nc]=permut(H_nc_p);
               [J_t_fc]=permut(H_fc_p);
       
               J_r_n = permut_JT(g_1_all{c}');
               J_r_f = permut_JT(g_2_all{c}');
               J_r_nc = permut_JT(g_b_all{c}');
               J_r_fc = permut_JT(g_b_all{c}');
       
               H_n{c}  = diag(g_1_all{c}'*J_r_n)*J_t_n*G_all*f1_all{c}*w_k(:, c);
            %    disp(H_n{c});
               H_f{c}  = diag(g_2_all{c}'*J_r_f)*J_t_f*G_all*f2_all{c}*w_k(:, c);
               H_n_c{c} = diag(g_b_all{c}'*J_r_nc)*J_t_nc*G_all*f1_all{c}*w_k(:, c);
               H_f_c{c} = diag(g_b_all{c}'*J_r_fc)*J_t_fc*G_all*f2_all{c}*w_k(:, c);  
    
                   A_n(c)+ delta_p >=1e-7;
                   B_n(c)+ delta_p >= 1e-7;
                   A_f(c)+ delta_p >= 1e-7;
                   B_f(c)+ delta_p >= 1e-7;
                   A_c_n(c)+ delta_p >= 1e-7;
                   B_c_n(c)+ delta_p >= 1e-7;               
                   
                   taylor_approx_far(c) = log2(1 + inv_pos(A_f_prev(c) * B_f_prev(c))) -  ...
                    (log2(exp(1)) * inv_pos(A_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (A_f(c) - A_f_prev(c)) - ...
                    (log2(exp(1)) * inv_pos(B_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (B_f(c) - B_f_prev(c));

    
                    taylor_approx_n(c)= log2(1 + inv_pos(A_n_prev(c) * B_n_prev(c))) -  ...
                                    (log2(exp(1)) * inv_pos(A_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (A_n(c) - A_n_prev(c)) - ...
                                    (log2(exp(1)) * inv_pos(B_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (B_n(c) - B_n_prev(c));


                        % Backscatter device constraints (Taylor approximation) at near user
                    taylor_approx_backscatter_n(c) = log2(1 + 1 ./ (A_c_prev_n(c) * B_c_prev_n(c))) - ...
                                                    (log2(exp(1)) / (A_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (A_c_n(c) - A_c_prev_n(c)) - ...
                                                    (log2(exp(1)) / (B_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (B_c_n(c) - B_c_prev_n(c));



                        R_f(c)- taylor_approx_far(c)<=delta_p;
                        R_n(c) - taylor_approx_n(c) <=delta_p;
                        R_c_n(c) - taylor_approx_backscatter_n(c)<=delta_p;

                       % (54b) and (54c)
                        delta_p>= R_f_min-R_f(c);  
                        delta_p>= R_n_min-R_n(c);
                        delta_p>= R_c_min- R_c_n(c);

                   inter_cluster_interference_near = 0;
                   inter_cluster_interference_far = 0;
                   inter_cluster_interference_near_b=0;
                   inter_cluster_interference_far_b=0;
                   for j = 1:numClusters
                       if j ~= c
                           % Near user inter cluster interference  
                        %    (size(V));
                        %    dispdisp(size((diag(g_1_all{c}'*J_r_n)*J_t_n*G_all*f1_all{c}*w_k(:, j))));
                           inter_cluster_interference_near = inter_cluster_interference_near + ...
                               real(trace(V * (diag(g_1_all{c}'*J_r_n)*J_t_n*G_all*f1_all{c}*w_k(:, j)) * (diag(g_1_all{c}'*J_r_n)*J_t_n*G_all*f1_all{c}*w_k(:, j))'));

                           inter_cluster_interference_far= inter_cluster_interference_far + ...
                               real(trace(V * (diag(g_2_all{c}'*J_r_f)*J_t_f*G_all*f2_all{c}*w_k(:, j)) * (diag(g_2_all{c}'*J_r_f)*J_t_f*G_all*f2_all{c}*w_k(:, j))')); 

                           inter_cluster_interference_near_b= inter_cluster_interference_near_b + ...
                               real(trace(V* (diag(g_b_all{c}'*J_r_nc)*J_t_nc*G_all*f1_all{c}*w_k(:, j)) * (diag(g_b_all{c}'*J_r_nc)*J_t_nc*G_all*f1_all{c}*w_k(:, j))')); 
                           
                           inter_cluster_interference_far_b= inter_cluster_interference_far_b + ...
                               real(trace(V * (diag(g_b_all{c}'*J_r_fc)*J_t_fc*G_all*f2_all{c}*w_k(:, j))*(diag(g_b_all{c}'*J_r_fc)*J_t_fc*G_all*f2_all{c}*w_k(:, j))')); 
               
                       end

                   end

                       % Define slack variables based on cascaded channel

                   delta_p >= inv_pos(A_n(c))-real(trace(V * H_n{c} * H_n{c}')) * alpha_n; % (59b)

                    % inter cluster interference  + backscatter interference + noise power
                   delta_p>=inter_cluster_interference_near + ...
                            real(trace(V * H_n_c{c} * H_n_c{c}')) * eta_k + para.noise- B_n(c) ;  %% (59c)

                       
                   delta_p>= inv_pos(A_f(c)) - real(trace(V * H_f{c} * H_f{c}')) * alpha_f; %% (59d)

                   % %% (59e)
                   delta_p>= inter_cluster_interference_far + ...
                           real(trace(V * H_f{c} * H_f{c}'))  * alpha_n + ...
                           real(trace(V * H_f_c{c} * H_f_c{c}'))   * eta_k + para.noise-B_f(c) ;

                   %% (59f)  
                   delta_p>=inv_pos(A_c_n(c)) - real(trace(V * H_n_c{c} * H_n_c{c}')) * eta_k;

                    %% (50g)  
                   delta_p >= inter_cluster_interference_near_b + para.noise- B_c_n(c) ;
                   %% (59j)

           end
        %    diag(V) == 1 + delta_p; 

          for m=1:N
            V(m,m) == 1 + delta_p;
          end

        %    V_eignemax'*V_eignemax >=epsln_1*trace(V); % (59g)


          V + delta_p * eye(N) == hermitian_semidefinite(N);
          
   cvx_end

       obj_prev = cvx_optval;
       A_n_opt = A_n;
       B_n_opt = B_n;
       A_f_opt = A_f;
       B_f_opt = B_f;
       A_c_n_opt = A_c_n;
       B_c_n_opt = B_c_n;
       V_opt = V;
       status = cvx_status;
       epsln = epsln_1; 

end