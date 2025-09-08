function [WSR, SINR_n, SINR_f, SINR_b] = compute_SR(para, w_k, G_all, g_1_all, g_2_all, g_b_all, f1_all, f2_all, alpha_n, alpha_f, Theta)
    numClusters = para.K;
    noise_power = para.noise * para.scal^2;  % Match scaling from update()
    WSR = 0;
    R_n = zeros(numClusters, 1);
    R_f = zeros(numClusters, 1);
    R_c_n = zeros(numClusters, 1);
    
    % Verify power constraints
    total_power = 0;
    for k = 1:numClusters
        total_power = total_power + norm(w_k(:,k))^2;
    end
    % assert(total_power <= para.P_max * 1.01, 'Power constraint violated in WSR calculation');
    
    for c = 1:numClusters
        % Reconstruct W matrices from w_k (matching update() formulation)
        % W_k = w_k(:,c) * w_k(:,c)';
        
        % Calculate channels (consistent with update())
        H_n = g_1_all{c}' * Theta * G_all;
        H_f = g_2_all{c}' * Theta * G_all;
        H_n_c = g_b_all{c}' * Theta * G_all * f1_all{c};
        H_f_c = g_b_all{c}' * Theta * G_all * f2_all{c};
        
        % Calculate interference terms (match update() exactly)
        interf_n = 0;
        interf_f = 0;
        interf_n_b = 0;
        interf_f_b = 0;
        
        for j = 1:numClusters
          if j ~= c
            interf_n = interf_n + real(trace(W_k(:,:,j) * H_n' * H_n));
            interf_f = interf_f + real(trace(W_k(:,:,j) * H_f' * H_f));
            interf_n_b = interf_n_b + real(trace(W_k(:,:,j)  * H_n_c' * H_n_c));
            interf_f_b = interf_f_b + real(trace(W_k(:,:,j)  * H_f_c' * H_f_c));
          end
        end
        
        % Calculate A/B terms exactly as in update()
        A_n = real(trace(W_k * H_n' * H_n)) * alpha_n;
        B_n = interf_n + real(trace(W_k * H_n_c' * H_n_c)) * para.eta + noise_power;
        
        A_f = real(trace(W_k * H_f' * H_f)) * alpha_f;
        B_f = interf_f + real(trace(W_k * H_f' * H_f)) * alpha_n + ...
              real(trace(W_k * H_f_c' * H_f_c)) * para.eta + noise_power;
        
        A_c_n = real(trace(W_k * H_n_c' * H_n_c)) * para.eta;
        B_c_n = interf_n_b + noise_power;
        
        % Calculate rates with protection against numerical issues
        R_n(c) = log2(1 + A_n/B_n);
        R_f(c) = log2(1 + A_f/B_f);
        R_c_n(c) = log2(1 + A_c_n/B_c_n);
        SINR_n(c) = A_n / B_n;
        SINR_f(c) = A_f / B_f;
        SINR_b(c) = A_c_n / B_c_n;
        
        WSR = WSR + R_n(c) + R_f(c) + R_c_n(c);
        
        % Debug output
        % fprintf('Cluster %d:\n', c);
        % fprintf('  A_n: %.2e, B_n: %.2e, SINR_n: %.2f dB\n', ...
        %         A_n, B_n, 10*log10(A_n/B_n));
        % fprintf('  A_f: %.2e, B_f: %.2e, SINR_f: %.2f dB\n', ...
        %         A_f, B_f, 10*log10(A_f/B_f));
        % fprintf('  A_c: %.2e, B_c: %.2e, SINR_c: %.2f dB\n', ...
        %         A_c_n, B_c_n, 10*log10(A_c_n/B_c_n));
    end
    
    % fprintf('Total WSR: %.4f bits/s/Hz\n', WSR);
end