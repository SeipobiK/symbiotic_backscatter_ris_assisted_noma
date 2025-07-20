function [P] = permut(x)
    [~, idx] = sort(x, 'descend');
    P = eye(length(x));
    P = P(idx, :);
end 
