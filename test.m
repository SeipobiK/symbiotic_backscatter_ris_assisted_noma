% Hard-coded convergence data from your log
iterations = 1:20;
objective_values = [1.9438, 1.2802, 0.9562, 0.7302, 0.6221, ...
                    0.5679, 0.5447, 0.5332, 0.5274, 0.5245, ...
                    0.5231, 0.5224, 0.5221, 0.5220, 0.5219, ...
                    0.5219, 0.5219, 0.5219, 0.5219, 0.5219];

% Plot
figure;
plot(iterations, objective_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Iteration');
ylabel('Objective Value');
title('Convergence of Algorithm');
