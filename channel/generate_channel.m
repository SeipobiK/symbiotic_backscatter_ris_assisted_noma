function [H,g,f] = generate_channel(para, BS_array, RIS_array)


epsilon = para.rician; % Rician factor

BS_loc = para.BS_loc;
userloc = para.userloc;

%% BS to RIS channel
% NLOS
H_NLOS = 1/sqrt(2) .* ( randn(para.N,para.M) + 1i*randn(para.N,para.M) );
disp('norm of H');
disp(norm(H_NLOS,'fro')^2);

% LOS
a_BR = steering_vector(BS_array, -BS_loc(2), -BS_loc(3));
a_RB = steering_vector(RIS_array, BS_loc(2), BS_loc(3));
H_LOS = a_RB*a_BR.';

% pathloss
path_loss = para.pathloss(BS_loc(1))';
path_loss = sqrt(10.^(- path_loss/10));
% path_loss=para.pth;
path_loss_H=path_loss;
disp(['H pathloss ',  num2str(path_loss_H)]);
disp(['BS to RIS Path Loss: ', num2str(path_loss)]);
H = path_loss .* (sqrt(epsilon/(epsilon+1)) * H_LOS + sqrt(1/(epsilon+1)) * H_NLOS);
disp('norm of H LOSS');
disp(norm(H_LOS,'fro')^2);





% %% RIS to users channel
% % NLOS
g_NLOS = 1/sqrt(2) .* ( randn(para.N,para.K) + 1i*randn(para.N,para.K,para.K_u));



% Cluster 1 (now cluster index 1 in 2nd dimension)
g_LOS = zeros(para.N, para.K, para.K_u);
for i = 1:3
    g_LOS(:,1,i) = steering_vector(RIS_array, userloc(i, 1, 2), userloc(i, 1, 3));
    disp('G norm');
    disp(norm(g_LOS(:,1,i),'fro')^2);
end
    
% Cluster 2 (now cluster index 2 in 2nd dimension)
for k = 1:3
    g_LOS(:,2,k) = steering_vector(RIS_array, userloc(k, 2, 2), userloc(k, 2, 3)); 
    disp('G norm cluster 2');
    disp(norm(g_LOS(:,1,i),'fro')^2);

end

g=zeros(para.N, para.K, para.K_u);
for i=1:3
    path_loss = para.pathloss(userloc(i, 1, 1));
    path_loss = sqrt(10.^( (-path_loss)/10));
    disp(['Distance Cluster 1 :', num2str(userloc(i, 1, 1))]);
    % path_loss=para.pth;
    g(:,1,i)=path_loss .* (sqrt(epsilon/(epsilon+1)) * g_LOS(:,1,i) + sqrt(1/(epsilon+1)) * g_NLOS(:,1,i));
    disp(['Cluster 1 User ', num2str(i), ' Path Loss: ', num2str(path_loss)]);
end

for k=1:3
    path_loss = para.pathloss(userloc(k, 2, 1));
    path_loss = sqrt(10.^( (-path_loss)/10));
    disp(['Distance cluster 2 :', num2str(userloc(i, 2, 1))]);
    % path_loss=para.pth;
    g(:,2,k)=path_loss .* (sqrt(epsilon/(epsilon+1)) * g_LOS(:,2,k) + sqrt(1/(epsilon+1)) * g_NLOS(:,2,k));
    disp(['Cluster 2 User ', num2str(k), ' Path Loss: ', num2str(path_loss)]);
end

dist_BD_User1_cluster1 = norm(para.BD_cluster1 - para.User1_cluster1);
dist_BD_User2_cluster1 = norm(para.BD_cluster1 - para.User2_cluster1);

dist_BD_User1_cluster2 = norm(para.BD_cluster2 - para.User1_cluster2);
dist_BD_User2_cluster2 = norm(para.BD_cluster2 - para.User2_cluster2);

disp('Distances:');
disp(['Cluster 1: BD to User1: ', num2str(dist_BD_User1_cluster1)]);
disp(['Cluster 1: BD to User2: ', num2str(dist_BD_User2_cluster1)]);
disp(['Cluster 2: BD to User1: ', num2str(dist_BD_User1_cluster2)]);
disp(['Cluster 2: BD to User2: ', num2str(dist_BD_User2_cluster2)]);

f=zeros(2,2);

path_loss_1_1= para.pathloss(dist_BD_User1_cluster1);
path_loss_1_2= para.pathloss(dist_BD_User2_cluster1);
path_loss_2_1= para.pathloss(dist_BD_User1_cluster2);
path_loss_2_2= para.pathloss(dist_BD_User2_cluster2);
path_loss_1_1 = sqrt(10.^(-path_loss_1_1/10));
path_loss_1_2 = sqrt(10.^(-path_loss_1_2/10));
path_loss_2_1 = sqrt(10.^(-path_loss_2_1/10));
path_loss_2_2 = sqrt(10.^(-path_loss_2_2/10));

disp(path_loss_1_1);
% path_loss_1_2=1;
% path_loss_2_1=1;
% path_loss_2_2=1;


% BD to User1 and User2 path loss
f(1,1) = path_loss_1_1* (randn + 1i*randn)/sqrt(2);
f(1,2) = path_loss_1_2* (randn + 1i*randn)/sqrt(2);
f(2,1) = path_loss_2_1* (randn + 1i*randn)/sqrt(2);
f(2,2) = path_loss_2_2* (randn + 1i*randn)/sqrt(2);

end