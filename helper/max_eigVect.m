function [V_max]=max_eigVect(V_opt)
        [V,D] = eig(V_opt);
     
        [max_eigenvalue, index] = max(diag(D));
        max_Veigenvector = sqrt(max_eigenvalue) * V(:, index); 
        V_max= max_Veigenvector'*max_Veigenvector;
        disp(V_max);

end