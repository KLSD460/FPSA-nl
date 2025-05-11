% Compute the StatRes for L1oLk problem

function kkt = StatRes_L1oSK(x,K,lambda,A,b) 
    n = length(x);
    idx0 = find(abs(x) < 1e-6);
    idx1 = find(abs(x) >= 1e-6);
    [nrm_K,idx_K] = norm_maxK(x,K);
    nabla_K = zeros(n,1);
    nabla_K(idx_K) = sign(x(idx_K));

    temp = (nrm_K+1e-6)*A'*(A*x - b) - (lambda*norm(x,1) + 0.5*norm(A*x - b)^2)*nabla_K;
    temp(idx1) = temp(idx1) + lambda*(nrm_K+1e-6)*sign(x(idx1));
    temp(idx0) = max(abs(temp(idx0)) - lambda*(nrm_K+1e-6),0).*sign(temp(idx0));

    kkt = norm(temp);
end