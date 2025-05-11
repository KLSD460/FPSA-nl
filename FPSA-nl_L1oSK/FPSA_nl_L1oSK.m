function [x,output] = FPSA_nl_L1oSK(A,b,lb,ub,K,lambda,x0,pm)
%% FPSA-nl for solving
% min_x\in S (lambda*||x||_1 + 0.5*||Ax - b||^2)/(||x||_(K)+deps)
% S = { lb <= x <= ub}
% h(x) = 0.5*||Ax - b||^2, g(x) = \|lambda*||x||_1, f(x) = ||x||_(K)

% Input:
% maxit: the maximum iteration numbers for the algoritm
% tol: the termination criterion
% sigma: the parameter for extrapolation
% gamma, q, T: the parameter for non-monotone linesearch

% Output: x

%% parameter
[m,n] = size(A);
if isfield(pm,'maxit'); maxit = pm.maxit ; else; maxit = 20*m; end
if isfield(pm,'tol'); tol = pm.tol; else; tol = 1e-8;end
if isfield(pm,'sigma'); sigma = pm.sigma; else; sigma = 1.7;end
if isfield(pm,'gamma'); gamma = pm.gamma; else; gamma = 2;end
if isfield(pm,'q'); q = pm.q; else; q = 0.9;end
if isfield(pm,'T'); T = pm.T; else; T = 5;end
if isfield(pm,'deps'); deps = pm.deps; else; deps = 1e-6;end
%% initialize

c1 = 1e-3; c2 = 1e-3;
g = @(x) lambda*norm(x,1);
f = @(x) sum(maxk(abs(x),K))+deps;
h = @(x) 0.5*norm(A*x - b)^2;
gradh = @(x) A'*(A*x - b);
phi = @(x,u,delta) h(x) + g(x) + norm(x-u)^2/(2*delta);
Pi = @(x,u,delta) phi(x,u,delta)/f(x);
output.tag = 0;
%% iteration
xold = x0; 
uold = x0;
fold = ones(T,1)*inf;
theta = Pi(xold,uold,1);
dx = norm(x0);
ghold = zeros(n,1);

for i = 1:maxit
    
    % compute y
    [~,idxk] = maxk(abs(xold),K);
    y = zeros(n,1);
    y(idxk) = sign(xold(idxk));
    
    % update x
    idx = mod(i-1,T);
    fold(idx+1) =theta;
    gh     = gradh(xold);
    delta = gamma*dx/norm(gh - ghold);
    for ii = 1:250
        d1 = uold + delta*(theta*y - gh);
        x = sign(d1).*max(abs(d1) - lambda*delta,0);
        x = min(max(x,lb),ub);    
        u = uold - sigma*(uold - x);
        fnew = Pi(x,u,delta);
        if fnew < max(fold) - c1*norm(x-xold)^2 - c2*norm(u - uold)^2
            break;
        end
        delta = delta*q;
    end

    % compute theta
    theta = Pi(x,u,delta);

    % record
    dx      = norm(x-xold);
    relerrx = dx/max(norm(xold),eps);
    
    if (relerrx < tol) && i >=2
        output.tag = 1;
        break;
    end
    xold = x;
    uold = u;
    ghold = gh;

end
end