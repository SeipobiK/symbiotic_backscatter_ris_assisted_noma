function [P] = permut(x)
    [~, idx] = sort(abs(x), 'descend');  % Sort by magnitude only
    P = eye(length(x));
    P = P(idx, :);
end