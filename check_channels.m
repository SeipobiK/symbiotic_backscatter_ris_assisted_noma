% 1. Create a completely new figure with barebones settings
f = figure('Visible','off', 'MenuBar','none', 'ToolBar','none', ...
           'Color','w', 'Renderer','painters', 'IntegerHandle','off');

% 2. Use primitive plotting (no fancy stuff)
h = plot(1:max_iter, avg_obj_history, '-o', ...  % Solid lines + circles
         'LineWidth', 3, ...                     % Thick lines
         'MarkerSize', 8, ...                    % Big markers
         'MarkerFaceColor','auto');               % Filled markers

% 3. Manual axis styling (no auto settings)
set(gca, 'XLim',[1 max_iter], ...
         'YLim',[0 max(avg_obj_history(:))*1.1], ...
         'Box','on', ...
         'FontSize',12);

% 4. Save using low-level print command
print(f, 'CONVERGENCE_PLOT_FINAL.png', '-dpng', '-r300', '-opengl');
close(f);

% 5. Verify the file exists
if exist('CONVERGENCE_PLOT_FINAL.png','file')
    disp('Plot definitely saved - check your current folder');
    winopen(pwd);  % Open folder to view
else
    disp('This should never happen - check filesystem permissions');
end