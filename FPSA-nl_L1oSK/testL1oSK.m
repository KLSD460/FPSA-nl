
% This program solves the following problem:
% min_x\in S (lambda*||x||_1 + 0.5*||Ax - b||^2)/(||x||_(K) + dpes)
% S = { lb <= x <= ub}
% and compares three algorithms: FPSA-nl, PGSA_BE, e-PSG

% References:
%
% Deren Han, Min Tao and Zihao Xia
% A Single-loop Proximal Subgradient Algorithm for A Class Structured Fractional Programs,
% https://arxiv.org/abs/2503.12176

clc;clear;
%% generate problem
n = 1024;
m = 64;
m_idx = 0;
trials = 1;
% L = [1,5,10,15,20];
% K = [4,8,12,16,20,24];
L = [1];
K = [4];

obj1 = zeros(length(L),length(K),trials);
cpu_time1 = zeros(length(L),length(K),trials);
err1 = zeros(length(L),length(K),trials);
nnz1 = zeros(length(L),length(K),trials);
kkt1 = zeros(length(L),length(K),trials);

maxit = 1000;
lambda = 1e-3;

for t = 1:trials
    m_idx = m_idx + 1;
    tt1 = 0; tt2 = 0; tt3= 0;
    for ii = 1:length(L)
        for jj = 1:length(K)
            fprintf("trail = %d, L = %d, K = %d\n",t,L(ii),K(jj));
            rng(46*t*ii*jj);
            ppm.sen_mat = 'Oversampled_DCT'; ppm.F = L(ii);
            ppm.sparsity = K(jj);
            ppm.MM = m; ppm.NN = n;
            [A] = data_generator_A(ppm);

            supp = randsample_separated(n,K(jj),2*L(ii));
            xg  = zeros(n,1); xg(supp) = sign(randn(K(jj),1))/2;
            b   = A*xg;
            lb = -1*ones(n,1); ub = 1*ones(n,1);
            deps = 1e-6;
            objval = @(x) (lambda*norm(x,1) + 0.5*norm(A*x - b)^2)/(sum(maxk(abs(x),K(jj)))+deps);
            x0    = xg + 0.20*(2*rand(n,1) - 1);

            %% FPSA-nl
            pm.maxit = 5000;
            pm.sigma = 1.35;
            pm.xg = xg;
            pm.tol = 1e-6;
            pm.gamma = 0.9;
            pm.deps = deps;
            pm.is_linesearch = 1;
            pm.q = 0.9;
            pm.T = 20;

            t1 = tic;
            [x1,output1] = FPSA_nl_L1oSK(A,b,lb,ub,K(jj),lambda,x0,pm);
            tt1 = toc(t1);

            cpu_time1(ii,jj,t) = cpu_time1(ii,jj,t) + tt1;
            err1(ii,jj,t) = norm(x1 - xg)/max(1,norm(xg));
            obj1(ii,jj,t) = objval(x1);
            nnz1(ii,jj,t) = nnz(x1);
            kkt1(ii,jj,t) = StatRes_L1oSK(x1,K(jj),lambda,A,b);
           
            %% Report
            fprintf("FPSA-nl: cpu: %.2e, \t err:%.2e,\t, obj: %.2e\t, nnz:%.1f,\n" ...
                ,cpu_time1(ii,jj,t),err1(ii,jj,t),obj1(ii,jj,t),nnz1(ii,jj,t));

        end

    end
end

ave_err1 = mean(err1,3);
ave_obj1 = mean(obj1,3);
ave_cpu_time1 = mean(cpu_time1,3);
ave_nnz1 = mean(nnz1,3);
ave_kkt1 = mean(kkt1,3);