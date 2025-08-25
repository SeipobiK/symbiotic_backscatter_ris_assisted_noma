function [V_opt,A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,status] = Passive_BF(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n,V_max,epsln_1)



   
   numClusters = para.K; % Number of clusters
   N = para.N; % Number of BS antennas
   alpha_n = para.alpha_k_n; % Near user path loss factor
   alpha_f = para.alpha_k_f; % Far user path loss factor
   R_n_min = para.R_min_n; % Minimum rate for near user
   R_f_min = para.R_min_f; % Minimum rate for far user
   R_c_min = para.R_c_min; % Minimum rate for backscatter user
   eta_k = para.eta; % Backscatter coefficient
   noise = para.noise * (para.scal)^2;  % Noise scales with power
   para.P_max = para.P_max;

   cvx_begin quiet sdp
       % cvx_solver sedumi
       cvx_solver mosek
    %    cvx_precision high
    %    cvx_precision high


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

       expressions taylor_approx_far(numClusters, 1) taylor_approx_n(numClusters, 1) taylor_approx_backscatter_n(numClusters, 1)
       

   


       % Objective function: Maximize weighted sum rate (59a)
       maximize(sum(R_n + R_f + R_c_n)) 

        subject to

  
           for c = 1:numClusters
                A_n(c)  >=1e-4; 
                B_n(c)  >= 1e-4;
                A_f(c) >= 1e-4;
                B_f(c) >= 1e-4;
                A_c_n(c) >= 1e-4;
                B_c_n(c) >= 1e-4;  
                
                A_n(c)  <=1e+2; 
                B_n(c)  <= 1e+2;
                A_f(c)  <= 1e+2;
                B_f(c)  <= 1e+2;
                A_c_n(c)  <= 1e+2;
                B_c_n(c)  <= 1e+2;  
        
               H_n{c}  = diag(g_1_all{c}')*G_all*w_k(:, c);
               H_f{c}  = diag(g_2_all{c}')*G_all*w_k(:, c);
               H_n_c{c} = diag(g_b_all{c})*G_all*f1_all{c}*w_k(:, c);
               H_f_c{c} = diag(g_b_all{c}')*G_all*f2_all{c}*w_k(:, c);  
            
                   
                %    taylor_approx_far(c) = log2(1 + inv_pos(A_f_prev(c) * B_f_prev(c))) -  ...
                %     (log2(exp(1)) * inv_pos(A_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (A_f(c) - A_f_prev(c)) - ...
                %     (log2(exp(1)) * inv_pos(B_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (B_f(c) - B_f_prev(c));

    
                %     taylor_approx_n(c)= log2(1 + inv_pos(A_n_prev(c) * B_n_prev(c))) -  ...
                %                     (log2(exp(1)) * inv_pos(A_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (A_n(c) - A_n_prev(c)) - ...
                %                     (log2(exp(1)) * inv_pos(B_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (B_n(c) - B_n_prev(c));


                %         % Backscatter device constraints (Taylor approximation) at near user
                %     taylor_approx_backscatter_n(c) = log2(1 + 1 ./ (A_c_prev_n(c) * B_c_prev_n(c))) - ...
                %                                     (log2(exp(1)) / (A_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (A_c_n(c) - A_c_prev_n(c)) - ...
                %                                     (log2(exp(1)) / (B_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (B_c_n(c) - B_c_prev_n(c));


                    R_f(c) <= log2(1 + 1 ./ (A_f_prev(c) * B_f_prev(c))) -  ...
                    (log2(exp(1)) * 1 ./ (A_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (A_f(c) - A_f_prev(c)) - ...
                    (log2(exp(1)) * 1 ./ (B_f_prev(c) * (1 + A_f_prev(c) * B_f_prev(c)))) * (B_f(c) - B_f_prev(c));

                    
                    R_n(c) <= log2(1 + 1 ./ (A_n_prev(c) * B_n_prev(c))) -  ...
                    (log2(exp(1)) * 1 ./ (A_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (A_n(c) - A_n_prev(c)) - ...
                    (log2(exp(1)) * 1 ./ (B_n_prev(c) * (1 + A_n_prev(c) * B_n_prev(c)))) * (B_n(c) - B_n_prev(c));



                    R_c_n(c) <= log2(1 + 1 ./ (A_c_prev_n(c) * B_c_prev_n(c))) - ...
                    (log2(exp(1)) *1 ./  (A_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (A_c_n(c) - A_c_prev_n(c)) - ...
                    (log2(exp(1)) * 1 ./  (B_c_prev_n(c) * (1 + A_c_prev_n(c) * B_c_prev_n(c)))) * (B_c_n(c) - B_c_prev_n(c)); 

                       % (54b) and (54c)
                        R_f(c)>= R_f_min;  
                        R_n(c)>= R_n_min;
                        R_c_n(c)>= R_c_min;

                   inter_cluster_interference_near = 0;
                   inter_cluster_interference_far = 0;
                   inter_cluster_interference_near_b=0;
                   for j = 1:numClusters
                       if j ~= c
                           % Near user inter cluster interference  
                        %    (size(V));
                        %    dispdisp(size((diag(g_1_all{c}')*G_all*f1_all{c}*w_k(:, j))));
                           inter_cluster_interference_near= inter_cluster_interference_near + ...
                               real(trace(V * (diag(g_1_all{c}')*G_all*w_k(:, j)) * (diag(g_1_all{c}')*G_all*w_k(:, j))')); 

                           inter_cluster_interference_far= inter_cluster_interference_far + ...
                               real(trace(V * (diag(g_2_all{c}')*G_all*w_k(:, j)) * (diag(g_2_all{c}')*G_all*w_k(:, j))')); 

                           inter_cluster_interference_near_b= inter_cluster_interference_near_b + ...
                               real(trace(V* (diag(g_b_all{c}')*G_all*f1_all{c}*w_k(:, j)) * (diag(g_b_all{c})*G_all*f1_all{c}*w_k(:, j))')); 
 
              
                       end

                   end

                       inv_pos(A_n(c))<=real(trace(V  * H_n{c} * H_n{c}')) * alpha_n; % Near user

                        % inter cluster interference  + backscatter interference + noise power
                       B_n(c)>=inter_cluster_interference_near + ...
                                real(trace(V * H_n_c{c} * H_n_c{c}')) * eta_k + noise;

                        
                       inv_pos(A_f(c)) <=real(trace(V * H_f{c} * H_f{c}')) * alpha_f; 


                       B_f(c) >= inter_cluster_interference_far + ...
                                real(trace(V  * H_f{c} * H_f{c}'))  * alpha_n + ...
                                real(trace(V  * H_f_c{c} * H_f_c{c}'))   * eta_k + noise;


                       inv_pos(A_c_n(c)) <=real(trace(V  * H_n_c{c} * H_n_c{c}')) * eta_k;

                       B_c_n(c)>= inter_cluster_interference_near_b + noise ;
           end

          for m=1:N
            V(m,m) == 1;
          end 
         
          V_max'*V*V_max>=epsln_1*real(trace(V)); % (59g)

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
end