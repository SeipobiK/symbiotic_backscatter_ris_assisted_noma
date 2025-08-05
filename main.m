% Clear workspace
clc;
clear;
close all;

% Define coordinates
BS = [0, 0, 15];          % Base Station (x, y, z)
RIS = [35, 20, 15];     % Ze-RIS (x, y, z)
cluster_centers = [30, 0, 0; 40, 0, 0];  % Cluster centers (x, y, z)
radii = [5, 5];           % Cluster radii (5m each)

% Define users in each cluster [x, y, z]
users_cluster1 = [27, 4, 0; 32, 0, 0; 30, 5, 0];  % User1, User2, BD
users_cluster2 = [43, 4, 0; 37, 0, 0; 40, 5, 0];  % User1, User2, BD

% Create figure
figure;
hold on;
grid on;
view(3); % 3D view
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Z-axis (m)');
title('Simulation Layout: BS, Ze-RIS, Clusters, and Users');

% Plot BS and Ze-RIS
plot3(BS(1), BS(2), BS(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(BS(1), BS(2), BS(3), 'BS', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot3(RIS(1), RIS(2), RIS(3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
text(RIS(1), RIS(2), RIS(3), 'Ze-RIS', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Plot clusters (circles in XY plane)
theta = linspace(0, 2*pi, 100);
for i = 1:size(cluster_centers, 1)
    x = cluster_centers(i, 1) + radii(i) * cos(theta);
    y = cluster_centers(i, 2) + radii(i) * sin(theta);
    z = zeros(size(theta)); % z=0 for all points
    plot3(x, y, z, 'k-', 'LineWidth', 1.5);
    text(cluster_centers(i, 1), cluster_centers(i, 2), cluster_centers(i, 3), ...
        sprintf('Cluster %d', i), 'VerticalAlignment', 'bottom');
end

% Plot users in Cluster 1 (green markers)
plot3(users_cluster1(1,1), users_cluster1(1,2), users_cluster1(1,3), 'g^', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot3(users_cluster1(2,1), users_cluster1(2,2), users_cluster1(2,3), 'g^', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot3(users_cluster1(3,1), users_cluster1(3,2), users_cluster1(3,3), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); % BD

% Plot users in Cluster 2 (magenta markers)
plot3(users_cluster2(1,1), users_cluster2(1,2), users_cluster2(1,3), 'm^', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
plot3(users_cluster2(2,1), users_cluster2(2,2), users_cluster2(2,3), 'm^', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
plot3(users_cluster2(3,1), users_cluster2(3,2), users_cluster2(3,3), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); % BD

% Label users
text(users_cluster1(1,1), users_cluster1(1,2), users_cluster1(1,3), 'User1', 'FontSize', 8);
text(users_cluster1(2,1), users_cluster1(2,2), users_cluster1(2,3), 'User2', 'FontSize', 8);
text(users_cluster1(3,1), users_cluster1(3,2), users_cluster1(3,3), 'BD1', 'FontSize', 8);

text(users_cluster2(1,1), users_cluster2(1,2), users_cluster2(1,3), 'User1', 'FontSize', 8);
text(users_cluster2(2,1), users_cluster2(2,2), users_cluster2(2,3), 'User2', 'FontSize', 8);
text(users_cluster2(3,1), users_cluster2(3,2), users_cluster2(3,3), 'BD2', 'FontSize', 8);

% Adjust axes
axis equal;
xlim([0, 50]); % Focus on cluster regions
ylim([-45, 45]);  % Users are along x-axis
zlim([0, 20]);

% Add legend
legend('BS', 'Ze-RIS', 'Clusters', 'Cluster1 Users', 'Cluster1 BD', 'Cluster2 Users', 'Cluster2 BD', 'Location', 'northeast');

hold off;