function [WSR,R_n,R_f,R_c_n,A_n,B_n] = Compute_WSR(para,w_k,G_all, g_1_all,...
    g_2_all,g_b_all,f1_all,f2_all, alpha_n, alpha_f, v)
         % power= abs(w_k(:,1))^2 + abs(w_k(:,2))^2;
         % disp(['Total power in Compute WSR: ', num2str(power)]);         
          numClusters = para.K; % Number of 
          N=para.N; % Number of reflecting elements at RIS
          alpha_f = para.alpha_k_f; % Far user weight
          alpha_n = para.alpha_k_n; % Near user weight
          WSR = zeros(numClusters, 1); % Initialize WSR vector
          A_n = zeros(numClusters, 1); % Initialize A_n vector
          B_n = zeros(numClusters, 1); % Initialize B_n vector
          A_f = zeros(numClusters, 1); % Initialize A_f vector
          B_f = zeros(numClusters, 1); % Initialize B_f vector
          A_c_n = zeros(numClusters, 1); % Initialize A_c_n vector
          B_c_n = zeros(numClusters, 1); % Initialize B_c_n vector
          R_n=zeros(numClusters, 1); % Initialize R_n vector
          R_f=zeros(numClusters, 1); % Initialize R_f vector    
          R_c_n=zeros(numClusters, 1); % Initialize R_c_n vector
          noise = para.noise* (para.scal)^2;  % Noise power
        %   G_all=G_all*para.scal; % Scale the channel matrix


           for c = 1:numClusters
                    inter_cluster_interference_near = 0;
                    inter_cluster_interference_far = 0;
                    inter_cluster_interference_near_b=0;
                    inter_cluster_interference_far_b=0;
                    for j = 1:numClusters
                        if j ~= c
                            % Near user inter cluster interference  
                            inter_cluster_interference_near = inter_cluster_interference_near + ...
                                abs(v'* (diag(g_1_all{c}')*G_all*w_k(:, j))).^2;

                            inter_cluster_interference_far= inter_cluster_interference_far + ...
                               abs(v'* (diag(g_2_all{c}')*G_all*w_k(:, j))).^2;

                            inter_cluster_interference_near_b= inter_cluster_interference_near_b + ...
                                abs(v'* (diag(g_b_all{c}')*G_all*f1_all{c}*w_k(:, j))).^2;
                                
                            inter_cluster_interference_far_b= inter_cluster_interference_far_b + ...
                                abs(v'* (diag(g_b_all{c}')*G_all*f2_all{c}*w_k(:, j))).^2;
                            
                        end

                    end
                    
                    A_n(c) = abs(v'* (diag(g_1_all{c}')*G_all*w_k(:, c))).^2 * alpha_n; % Near user

                    B_n(c) =inter_cluster_interference_near + ...
                            abs(v'* (diag(g_2_all{c}')*G_all*w_k(:, c))).^2 * para.eta + noise;

                        
                    A_f(c) = abs(v'* (diag(g_2_all{c}')*G_all*w_k(:, c))).^2 * alpha_f;


                    B_f(c) = inter_cluster_interference_far + ...
                            abs(v'* (diag(g_2_all{c}')*G_all*w_k(:, c))).^2 * alpha_n + ...
                            abs(v'* (diag(g_b_all{c}')*G_all*f2_all{c}*w_k(:, c))).^2 * para.eta + noise;


                    A_c_n(c) = abs(v'* (diag(g_b_all{c}')*G_all*f1_all{c}*w_k(:, c))).^2 * para.eta;

                    B_c_n(c) = inter_cluster_interference_near_b + noise;

                    % compute rates
                    R_f(c) = log2(1 + A_f(c) / B_f(c)); % Far user rate
                    R_n(c) = log2(1 + A_n(c) / B_n(c)); % Near user rate
                    R_c_n(c) = log2(1 + A_c_n(c) / B_c_n(c)); % Backscatter rate

            end

            % Compute WSR
            WSR=0;
            for c = 1:numClusters
                WSR =WSR+ R_f(c) + R_n(c) + R_c_n(c); % Weighted sum rate
            end
end
