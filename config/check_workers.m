function setupParallelPool()
    % SETUPPARALLELPOOL - Optimized parallel pool setup for 12-core Intel Ultra 7
    
    try
        % Get cluster information
        cluster = parcluster('local');
        
        % System specs from your screenshot
        physicalCores = 12;    % Your CPU has 12 physical cores
        logicalCores = 14;     % And 14 logical threads
        totalRAM = 32;         % 32GB RAM
        
        % Calculate recommended workers
        % Rule of thumb: 1 worker per physical core, leave 2 cores for system
        recommendedWorkers = min([physicalCores - 2, feature('numcores')]);
        
        % Memory considerations (~2GB RAM per worker as starting point)
        ramPerWorker = 2; % GB
        maxWorkersByRAM = floor(totalRAM / ramPerWorker);
        finalWorkerCount = min([recommendedWorkers, maxWorkersByRAM]);
        
        fprintf('System Detected:\n');
        fprintf(' - Intel Ultra 7 155U (12 cores/14 threads)\n');
        fprintf(' - 32GB RAM\n\n');
        fprintf('Recommended parallel pool settings:\n');
        fprintf(' - %d workers (1 per physical core, minus 2 for system)\n', recommendedWorkers);
        fprintf(' - Up to %d workers possible by RAM\n', maxWorkersByRAM);
        
        % Check existing pool
        currentPool = gcp('nocreate');
        if ~isempty(currentPool)
            if currentPool.NumWorkers == finalWorkerCount
                fprintf('\nParallel pool already running with optimal %d workers.\n', finalWorkerCount);
                return;
            else
                delete(currentPool);
                fprintf('\nClosing existing pool with %d workers.\n', currentPool.NumWorkers);
            end
        end
        
        % Adjust cluster settings
        cluster.NumWorkers = finalWorkerCount;
        cluster.NumThreads = 1; % Best for CPU-bound tasks
        
        % Create pool
        fprintf('\nStarting new parallel pool with %d workers...\n', finalWorkerCount);
        parpool(cluster, finalWorkerCount);
        
        % Verify
        newPool = gcp();
        fprintf('Successfully created pool with %d workers.\n', newPool.NumWorkers);
        fprintf('Each worker has %.1fGB RAM available.\n', totalRAM/newPool.NumWorkers);
        
    catch err
        fprintf('\nError setting up parallel pool:\n');
        fprintf('%s\n', err.message);
        
        % Fallback suggestions
        fprintf('\nTry these alternatives:\n');
        fprintf('1) parpool(10) - Slightly fewer workers\n');
        fprintf('2) parpool(''Processes'', 8) - Force process-based workers\n');
        fprintf('3) Check with IT if Parallel Computing Toolbox is licensed\n');
    end
end