function [Theta, Z, V, lsfnval, iter] = njgcg(S, lambda1, lambda2, ...
    omega1, omega2, omega3, show, nlist, p, K, rho, tau, max_iter)
%% ADMM algorithm parameters
if nargin < 7 || isempty(show)
    show = 0;
end
if nargin < 8 || isempty(nlist)
    nlist = ones(1, numel(S)) * 1;
end
if nargin < 9 || isempty(p)
    p = size(S{1}, 1);
end
if nargin < 10 || isempty(K)
    K = numel(S);
else
    assert(K == numel(S), "The input K is not equal to numel(S).")
end
if nargin < 11 || isempty(rho)
    rho = 2.5;
end
if nargin < 12 || isempty(tau)
    tau = 1e-10;
end
if nargin < 13 || isempty(max_iter)
    max_iter = 1500;
end

for k = 1: K
    assert(size(S{k}, 1) == size(S{k}, 2), ['A matrix in S is not a square matrix: k = ', num2str(k)])
    assert(size(S{k}, 1) == p, ['The number of rows in matrix in S is not all equal: k = ', ...
        num2str(k)])
    assert(size(S{k}, 2) == p, ['The number of columns in matrix in S is not all equal: k = ', ...
        num2str(k)])
end

%% Other parameters
iter = 1;
sc_val = 1e10;

%% Init
[Theta, Z, V, Theta_hat, Z_hat, V_hat, Wtheta, Wz, Wv] = init(p, K);
Theta_old = cell2mat(Theta);

%% ADMM begin
lsfnval = zeros(max_iter, 6);
while (( sc_val > tau) && (iter < max_iter))
    %% Update B{theta, Z, V}
    [Theta, Z, V] = nextB(S, Theta_hat, Z_hat, V_hat, Wtheta, Wz, Wv, ...
        nlist, p, K, lambda1, lambda2, omega1, omega2, omega3, rho);

    %% Update B_hat{Theta_hat, Z_hat, V_hat}
    [Theta_hat, Z_hat, V_hat] = nextB_hat(Theta, Z, V, Wtheta, Wz, Wv, ...
        K, rho);

    %% Update W{Wtheta, Wz, Wv}
    [Wtheta, Wz, Wv] = nextW(Theta, Z, V, Theta_hat, Z_hat, V_hat, ...
        Wtheta, Wz, Wv, K);

    %% Calculate loss
    lsfnval(iter, :) = lossFunction(S, Theta, Z, V, ...
        nlist, p, K, lambda1, lambda2, omega1, omega2, omega3);
    
    %% Check convergence conditions
    V_t = V;
    for k = 1: K
        V_t{k} = V_t{k}';
    end
%     ######################### Theta = Z + V + V' #########################
%     Theta_new = cell2mat(Z) + cell2mat(V) + cell2mat(V_t);
    Theta_new = cell2mat(Theta);
    sc_val = sum(sum((Theta_new - Theta_old).^2)) ./ sum(sum(Theta_old.^2));
%     Theta_old = Theta_new;
%     Theta_old = cell2mat(Theta);
    Theta_old = cell2mat(Z) + cell2mat(V) + cell2mat(V_t);
    
%     #####################################################################
    % 查看迭代过程
    if show == -2
%         Theta = mat2cell(Theta_new, ones(K, 1) .* p, p);
%         show_results(Theta, Z, V, Theta_hat);
        show_results([Theta, Z, V, Theta_hat, Z_hat, V_hat, Wtheta, Wz, Wv])
    end
%     #####################################################################
    
    if (mod(iter, 10) == 0 && show < 0)
        disp(['Current Iteration: ' , num2str(iter)]);
%         show_results(Theta, Theta, Z, V);
    end
    iter = iter + 1;
    
    if (iter == max_iter && show > 0)
        disp(['Maximum number of iteration reached, hglmv may not converge.', ...
            '-', num2str(lambda1), '-', num2str(lambda2)])
    end
end
end
