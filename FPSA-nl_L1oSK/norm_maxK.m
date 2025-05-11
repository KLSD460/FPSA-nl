% compute ||x||_(K)

function [y,idx] = norm_maxK(x,K)
    [temp,idx] = sort(abs(x),"descend");
    idx        = idx(1:K);
    y    = sum(temp(1:K));
end