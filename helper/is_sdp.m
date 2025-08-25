function isSDP = is_sdp(A)
    isSymmetric = issymmetric(A);
    eigenvalues = eig(A);
    isPositiveSemidefinite = all(eigenvalues >= -1e-10);  % Allow small numerical errors
    isSDP = isSymmetric && isPositiveSemidefinite;
end
