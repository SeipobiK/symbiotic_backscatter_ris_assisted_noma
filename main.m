%Simulation codes for the unicast scenario
%Date: 01/10/2022
%Author: Xidong Mu

clear all;
tic;
%The number of users is K, the number of antennas at the BS is N, the number of reflection elements is M, the user is single-antenna
para = para_init();
K=3;
N=para.N;
M=para.M;
PL=12;
K_BS=2;%The Rician factor of the BS-RIS channel, 3dB
a_BS=2.2;%The path-loss exponent of the BS-RIS channel
K_SU=2;%The Rician factor of the RIS-user channel, 3dB
a_SU=2.2;%The path-loss exponent of the RIS-use channel
%Initialize the locations of the BS, the RIS, and K users
% para.noise = para.noise * (1e+4^2);  % Noise scales with power
% para.P_max = para.P_max / 1e+4;
location_BS=[-50;0;0];
location_RIS=[0;0;0];
%The users are uniformly distributed in a circle centered at the RIS with the radius of 3 m
radius(1)=3;
radius(2)=3;
maxOut=10;%Maximum number of outer iterations
maxIn=10;%Maximum number of inner iterations
c=10;
for k=1:K
    QoS(k)=10^(0/10);
end
numClusters=3; % Number of clusters

% para.RIS_size = [1,3]; % elements at RIS
% para.N = para.RIS_size(1)*para.RIS_size(2);

% [BS_array, RIS_array] = generate_arrays(para);
% [G, r, f_1, f_2] = generate_channel(para, BS_array, RIS_array);      


k = 3; % number of sets
H_all = cell(1, k);
g_1_all = cell(1, k);
g_2_all = cell(1, k);
g_b_all = cell(1, k);
f1_all = cell(1, k);
f2_all = cell(1, k);

H_n=cell(1, k); 
H_f=cell(1, k);
H_nc=cell(1, k);
H_fc=cell(1, k);
u = zeros(M, k); % Initialize u for phase shifts
% Generate channel realizations

        for i = 1:k
            [BS_array, RIS_array] = generate_arrays(para);
            [H, g_1,g_2,g_b,f_1_3,f_2_3] = generate_channel(para, BS_array, RIS_array);
            H_all{i} = H;
            g_1_all{i} = g_1;
            g_2_all{i} = g_2;
            g_b_all{i} = g_b;
            f1_all{i} = f_1_3;
            f2_all{i} = f_2_3;
            disp(['H_all size: ', num2str(size(H_all{i}))]);
  
          % Generate reflection matrix Theta
            Theta=[];
            for m = 1:N
                u(m,k) = exp(1i*pi*(2*rand(1,1)));  % Random phase shifts
            end
            Theta(:,:,k) = diag(u(:,k),0);  % M x M diagonal reflection matrix'
            % disp(['Theta size: ', num2str(size(Theta(:,:,k)))]);



            scal=1e+4;
            H_n{i} = g_1_all{i}' * Theta(:,:,k) * H_all{i}*scal;
            H_f{i} = g_2_all{i}' * Theta(:,:,k) * H_all{i}*scal;
            H_nc{i} = g_b_all{i}' * Theta(:,:,k) * H_all{i} * f1_all{i}*scal;
            H_fc{i} = g_b_all{i}' * Theta(:,:,k) * H_all{i} * f2_all{i}*scal;


            % H_n{i} = g_1_all{i}' * Theta(:,:,k) * H_all{i}/norm(H_n{i}, 'fro');
            % H_f{i} = g_2_all{i}' * Theta(:,:,k) * H_all{i}/norm(H_f{i}, 'fro');
            % H_nc{i} = g_b_all{i}' * Theta(:,:,k) * H_all{i} * f1_all{i}/norm(H_nc{i}, 'fro');
            % H_fc{i} = g_b_all{i}' * Theta(:,:,k) * H_all{i} * f2_all{i}/norm(H_fc{i}, 'fro');

            disp(['H_n NORM : ', num2str(norm(H_fc{i}, 'fro'))]);

            % disp(['H_n : ', num2str(H_fc{i})]);

            % disp(['H_f size: ', num2str(size(H_f{i}))]);
            % disp(['H_nc size: ', num2str(size(H_nc{i}))]);
            % disp(['H_fc size: ', num2str(size(H_fc{i}))]);
            % Ensure H_n and H_f are ordered correctly
        

            if norm(H_n{i}) < norm(H_f{i})  % Correct
                H_temp = H_n{i};
                H_n{i} = H_f{i};
                H_f{i} = H_temp;
            end
         
        end

        W_init = zeros(M, M, numClusters);
        P_alloc = para.P_max / numClusters;

        for k = 1:numClusters
            W_init(:,:,k) = P_alloc * H_n{k} * H_n{k}';  % M×M Hermitian matrix
        end



        inter_cluster_interference_near = 0;
        inter_cluster_interference_far = 0;
        inter_cluster_interference_near_b=0;
        inter_cluster_interference_far_b=0;

        numClusters = 3; % Number of clusters
        M = para.M; % Number of BS antennas
        alpha_n = para.alpha_k_n; % Near user path loss factor
        alpha_f = para.alpha_k_f; % Far user path loss factor
        R_n_min = para.R_min_n; % Minimum rate for near user
        R_f_min = para.R_min_f; % Minimum rate for far user
        R_c_min = para.R_c_min; % Minimum rate f    or backscatter user
        eta_k = para.eta; % Backscatter coefficient
        P_max = para.P_max; % Maximum transmit power

        A_n_prev = zeros(k,1);
        B_n_prev = zeros(k,1);
        A_f_prev = zeros(k,1);
        B_f_prev = zeros(k,1);
        A_c_prev_n = zeros(k,1);
        B_c_prev_n = zeros(k,1);





        for c = 1:k  
                        % disp(size(H_n{c}));
                        % disp(size(H_f{c}));
                        % disp(size(W_init(:,:,c)));
                        % gain_f = norm(H_f{c} * W_init(:,:,c) , 'fro')^2;
                        % gain_n = norm(H_n{c} * W_init{c}, 'fro')^2;

                        % if gain_n < gain_f
                        %     disp('noooo');
                        %     % Swap if ordering is incorrect
                        %     temp = H_n{c};
                        %     H_n{c} = H_f{c};
                        %     H_f{c} = temp;

                        %     temp = H_nc{c};
                        %     H_nc{c} = H_fc{c};
                        %     H_fc{c} = temp; 
                        % end

                        % gain_2 = norm(H_f{c} * W_init{c}, 'fro')^2;
                        % gain_1 = norm(H_n{c} * W_init{c}, 'fro')^2;

                        % disp(['Gain 1: ', num2str(gain_1)]);
                        % disp(['Gain 2: ', num2str(gain_2)]);

                        inter_cluster_interference_near = 0;
                        inter_cluster_interference_far = 0;
                        inter_cluster_interference_near_b=0;
                        inter_cluster_interference_far_b=0;
                        for j = 1:numClusters
                            if j ~= c
                                % Near user inter cluster interference
                                inter_cluster_interference_near = inter_cluster_interference_near + ...
                                    real(trace(W_init(:,:,j)  * H_n{c}' * H_n{c}));

                                inter_cluster_interference_far= inter_cluster_interference_far + ...
                                    real(trace(W_init(:,:,j)   * H_f{c}' * H_f{c})); 

                                inter_cluster_interference_near_b= inter_cluster_interference_near_b + ...
                                    real(trace(W_init(:,:,j)   * H_nc{c}' * H_nc{c})); 
                                
                                inter_cluster_interference_far_b= inter_cluster_interference_far_b + ...
                                    real(trace(W_init(:,:,j)   * H_fc{c}' * H_fc{c})); 
                                
                            end

                        end
                     

                        A_n_prev(c) = 1/(real(trace(W_init(:,:,c) * H_n{c}' * H_n{c})) * alpha_n); % Near user


                        % inter cluster interference  + backscatter interference + noise power
                        B_n_prev(c) =inter_cluster_interference_near + ...
                                real(trace(W_init(:,:,c)* H_nc{c}' * H_nc{c})) * eta_k + para.noise;

                        R_N=log2(1+1/(A_n_prev(c)*B_n_prev(c)));

                        % disp(['R_N: ', num2str(R_N)]);


                        A_f_prev(c) = 1/(real(trace(W_init(:,:,c)* H_f{c}' * H_f{c})) * alpha_f); %(c)%;


                        B_f_prev(c) = inter_cluster_interference_far + ...
                                real(trace(W_init(:,:,c)* H_f{c}' * H_f{c}))  * alpha_n + ...
                                real(trace(W_init(:,:,c) * H_fc{c}' * H_fc{c}))   * eta_k + para.noise;

                        A_c_prev_n(c) = 1/(real(trace(W_init(:,:,c) * H_nc{c}' * H_nc{c})) * eta_k);
                        B_c_prev_n(c) = inter_cluster_interference_near_b + para.noise;


                        B_c_prev_n(c) = max(1e-6,inter_cluster_interference_near_b + para.noise);

                        % R_F=log2(1+1/(A_f_prev(c)*B_f_prev(c)));

                        % disp(['R_F: ', num2str(R_F)]);

       end               



       for c = 1:numClusters
            %    % Initialize A_prev, B_prev to feasible values
            A_n_prev(c) = 1e+0;  
            B_n_prev(c) = 1e-0;  
            A_f_prev(c) =1e-0;  
            B_f_prev(c) = 1e+0;  
            A_c_prev_n(c) = 1e+0;
            B_c_prev_n(c) = 1e-0;
            B_c_prev_f(c) = 1e-0; % Initialize B_c_prev_f for backscatter devices at far user
            A_c_prev_f(c) = 1e+0; % Initialize A_c_prev_f for backscatter devices at far user
       end

    
       disp('Initial values:');
       disp(['A_n_prev: ', num2str(A_n_prev')]);
       disp(['B_n_prev: ', num2str(B_n_prev')]);
       disp(['A_f_prev: ', num2str(A_f_prev')]);
       disp(['B_f_prev: ', num2str(B_f_prev')]);
       disp(['A_c_prev_n: ', num2str(A_c_prev_n')]);
       disp(['B_c_prev_n: ', num2str(B_c_prev_n')]);
       
    % Feasible Initial Point Search
    for n=1:1

        [W_optp, A_n_optp, B_n_optp, A_f_optp, B_f_optp, A_c_n_optp, B_c_n_optp,obj_prev,cvx_status] = feasible(para,W_init,H_n, H_f, H_nc, H_fc,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);
        
      
        if(strcmp(cvx_status,'Infeasible')||strcmp(cvx_status,'Failed'))
            disp('Infeasible problem, try to change the parameters.');  
            disp(cvx_status);
            % disp(W_opt);
            break;


        elseif strcmp(cvx_status,'Solved')
            % disp(cvx_status);
            % disp(obj_prev);

            for m= 1:20
                [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = feasible(para,W_init,H_n, H_f, H_nc, H_fc,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);

                A_n_prev = A_n_opt; 
                B_n_prev = B_n_opt; 
                A_f_prev = A_f_opt; 
                B_f_prev = B_f_opt; 
                A_c_prev_n = A_c_n_opt; 
                B_c_prev_n = B_c_n_opt; 
                W_init=W_opt;
                % disp(W_init);
                                
                disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_opt')]);
                disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_opt')]);     
                disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_opt')]);
                disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_opt')]);
                disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_n_opt')]);
                disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_n_opt')]);
                disp(['Iteration: ', num2str(m), ', Objective Value: ', sprintf('%.10f', obj_prev)]);
                disp(['Iteration: ', num2str(m), ', Status: ', num2str(cvx_status)]);
                disp(['Iteration: ', num2str(m), ', Rank ', num2str(rank(W_opt(:,:,1)))]);  
                % disp(W_opt);

                if strcmp(cvx_status, 'Solved') && abs(obj_prev) < 1e-5
                    disp('Convergence achieved.');
                    break;
                end
            end
        else
            disp('Problem not solved, check the parameters.');
        end

    end



    [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_prev,cvx_status] = update(para,W_init,H_n, H_f, H_nc, H_fc,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);
    disp(cvx_status);


    obj_array = zeros(1, 100);   % Pre-allocate for speed
    for m= 1:100
        [W_opt, A_n_opt, B_n_opt, A_f_opt, B_f_opt, A_c_n_opt, B_c_n_opt,obj_curr,cvx_status] = update(para,W_init,H_n, H_f, H_nc, H_fc,A_n_prev, B_n_prev, A_f_prev, B_f_prev,  A_c_prev_n, B_c_prev_n);

        A_n_prev = A_n_opt; 
        B_n_prev = B_n_opt; 
        A_f_prev = A_f_opt; 
        B_f_prev = B_f_opt; 
        A_c_prev_n = A_c_n_opt; 
        B_c_prev_n = B_c_n_opt; 
        W_init=W_opt;
        disp(W_init);
        obj_array(m) = obj_curr;

                        
        disp(['Iteration: ', num2str(m), ' A_n_opt: ', num2str(A_n_opt')]);
        disp(['Iteration: ', num2str(m), ' B_n_opt: ', num2str(B_n_opt')]);     
        disp(['Iteration: ', num2str(m), ' A_f_opt: ', num2str(A_f_opt')]);
        disp(['Iteration: ', num2str(m), ' B_f_opt: ', num2str(B_f_opt')]);
        disp(['Iteration: ', num2str(m), ' A_c_n_opt: ', num2str(A_c_n_opt')]);
        disp(['Iteration: ', num2str(m), ' B_c_n_opt: ', num2str(B_c_n_opt')]);
        disp(['Iteration: ', num2str(m), ', Objective Value: ', sprintf('%.10f', obj_array(m))]);
        if(m>1)
            disp(['Iteration: ', num2str(m), ', Objective Change: ', sprintf('%.10f', abs(obj_array(m) - obj_array(m-1)))]);
        end 
        disp(['Iteration: ', num2str(m), ', Status: ', num2str(cvx_status)]);
        disp(['Iteration: ', num2str(m), ', Rank ', num2str(rank(W_opt(:,:,1)))]);  
        % disp(W_opt);

        if m > 1 && abs(obj_array(m) - obj_array(m-1)) < 1e-5 && strcmp(cvx_status, 'Solved')
            disp(['Convergence achieved at iteration ', num2str(m), '.']);
            break;
        end
    end 