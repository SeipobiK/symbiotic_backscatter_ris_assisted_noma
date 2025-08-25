% Clear workspace
clc;
clear;
close all;


% Define coordinates
BS = [0, 0, 15];          % Base Station (x, y, z)
RIS = [45, 45, 15];   % Ze-RIS (x, y, z)
cluster_centers = [30, 0, 0; 40, 0, 0];  % Cluster centers (x, y, z)
radii = [5, 5];           % Cluster radii (5m each)

    % Cluster 1 users [User1, User2, BD]
    users_cluster1 = [27, 4, 0; 
                      32, 0, 0; 
                      30, 5, 0];
    
    % Cluster 2 users [User1, User2, BD]
    users_cluster2 = [43, 4, 0; 
                      37, 0, 0; 
                      40, 5, 0];

% ====================== Plot Configuration ======================
figure('Position', [100, 100, 800, 600], 'Color', 'w'); % Larger figure size
hold on; grid on; box on;
view(3); % 3D view
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'GridAlpha', 0.3);

% Axis labels and title
xlabel('X-axis (m)', 'FontSize', 14);
ylabel('Y-axis (m)', 'FontSize', 14);
zlabel('Z-axis (m)', 'FontSize', 14);
title('System Layout: BS, RIS, Clusters, and Users', 'FontSize', 16);

% ====================== Plot BS and RIS ======================
plot3(BS(1), BS(2), BS(3), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
text(BS(1), BS(2), BS(3), ' BS', 'FontSize', 12, 'VerticalAlignment', 'bottom');

plot3(RIS(1), RIS(2), RIS(3), 'bs', 'MarkerSize', 12, 'MarkerFaceColor', 'b', 'LineWidth', 2);
text(RIS(1), RIS(2), RIS(3), ' RIS', 'FontSize', 12, 'VerticalAlignment', 'bottom');

% ====================== Plot Clusters (Circles) ======================
theta = linspace(0, 2*pi, 100);
for i = 1:size(cluster_centers, 1)
    x = cluster_centers(i, 1) + radii(i) * cos(theta);
    y = cluster_centers(i, 2) + radii(i) * sin(theta);
    z = zeros(size(theta));
    plot3(x, y, z, 'k-', 'LineWidth', 1.5);
    text(cluster_centers(i, 1), cluster_centers(i, 2), cluster_centers(i, 3), ...
        sprintf('Cluster %d', i), 'FontSize', 12, 'HorizontalAlignment', 'center');
end

% ====================== Plot Users ======================
% Cluster 1 Users (Green markers)
plot3(users_cluster1(:,1), users_cluster1(:,2), users_cluster1(:,3), ...
    'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
text(users_cluster1(1,1), users_cluster1(1,2), users_cluster1(1,3), ...
    ' User1', 'FontSize', 10, 'Color', 'k');
text(users_cluster1(2,1), users_cluster1(2,2), users_cluster1(2,3), ...
    ' User2', 'FontSize', 10, 'Color', 'k');
text(users_cluster1(3,1), users_cluster1(3,2), users_cluster1(3,3), ...
    ' BD1', 'FontSize', 10, 'Color', 'k');

% Cluster 2 Users (Magenta markers)
plot3(users_cluster2(:,1), users_cluster2(:,2), users_cluster2(:,3), ...
    'm^', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'LineWidth', 1.5);
text(users_cluster2(1,1), users_cluster2(1,2), users_cluster2(1,3), ...
    ' User1', 'FontSize', 10, 'Color', 'k');
text(users_cluster2(2,1), users_cluster2(2,2), users_cluster2(2,3), ...
    ' User2', 'FontSize', 10, 'Color', 'k');
text(users_cluster2(3,1), users_cluster2(3,2), users_cluster2(3,3), ...
    ' BD2', 'FontSize', 10, 'Color', 'k');

% ====================== Adjust Axes and Legend ======================
axis equal;
xlim([0, 50]); ylim([-5, 50]); zlim([0, 20]);
legend({'Base Station (BS)', 'Reconfigurable RIS', 'Cluster Boundary', ...
    'Cluster 1 Users', 'Cluster 2 Users'}, ...
    'Location', 'northeast', 'FontSize', 10);

% ====================== Export for Paper ======================
exportgraphics(gcf, 'SystemLayout.png', 'Resolution', 300); % High-res PNG