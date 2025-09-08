function demonstrate_numerical_issues()
    % This function demonstrates numerical issues due to poor scaling

    % Reset the environment
    clear all; close all; clc;

    % Create a simple problem with huge coefficient range
    m = 3; % number of constraints
    n = 2; % number of variables

    % Define data with huge spread (from 1e-9 to 1e+9)
    A = [1e-9, 1e+9;   % Constraint matrix with huge spread
         1e+0, 1e+0;
         1e+9, 1e-9];
    b = [1; 2; 3];     % Right-hand side
    c = [1e+9, 1e-9];  % Objective coefficients with huge spread

    % Define MOSEK-specific options for dumping the problem
    dump_filename = 'badly_scaled_problem.opf';

    fprintf('Solving a problem with huge coefficient range:\n');
    fprintf('Constraint matrix A:\n');
    disp(A);
    fprintf('Right-hand side b:\n');
    disp(b);
    fprintf('Objective coefficients c:\n');
    disp(c);
    fprintf('Range ratio in A: %.2e\n', max(abs(A(:)))/min(abs(A(abs(A)>0))));
    fprintf('Range ratio in c: %.2e\n', max(abs(c))/min(abs(c)));

    % Solve with CVX and MOSEK
    cvx_begin
        cvx_solver mosek
        variable x(n)
        minimize( c * x )
        subject to
            A * x == b
            x >= 0
    cvx_end

    % Display results
    fprintf('CVX status: %s\n', cvx_status);
    fprintf('Optimal value: %e\n', cvx_optval);
    fprintf('Solution x:\n');
    disp(x);

    % Dump the problem for analysis
    fprintf('Dumping problem to %s for analysis...\n', dump_filename);
    cvx_begin
        cvx_solver mosek
        cvx_solver_settings('MSK_DPAR_DATA_TOL_AIJ_LARGE', 1e20, ... % Disable large coefficient warning
                            'MSK_DPAR_DATA_TOL_AIJ_SMALL', 0)        % Disable small coefficient warning
        variable x_dump(n)
        minimize( c * x_dump )
        subject to
            A * x_dump == b
            x_dump >= 0
        % cvx_solver_settings('WRITE_DEBUG_DATA', 1, ...
        %                     'WRITE_DATA_FORMAT', 'OP', ...
        %                     'WRITE_NAME_GENERAL', dump_filename)
    cvx_end

    fprintf('Dump complete. Analyze with: mosek -anapro %s\n', dump_filename);
end