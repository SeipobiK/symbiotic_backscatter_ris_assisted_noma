function [P] = permut_JT(x)
    [~, idx] = sort(abs(x), 'ascend');  % Sort by magnitude only
    P = eye(length(x));
    P = P(idx, :);
end