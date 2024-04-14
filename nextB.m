function [Theta, Z, V] = nextB(S, Theta_hat, Z_hat, V_hat, Wtheta, Wz, Wv, ...
    nlist, p, K, lambda1, lambda2, omega1, omega2, omega3, rho)
% B[t+1] <- argmin_B L{B, B_hat[t], W[t]}

%% Next Theta
Theta = cell(K, 1);
for k = 1: K
    nk = nlist(k);
    temp = S{k} - (rho / nk) .* (Theta_hat{k} - Wtheta{k});
    [U,D] = eig(temp);
    dhat = -D + sqrt(D.^2 + (4 * rho / nk) .* eye(p));
    Theta{k} = U * ((nk / 2 / rho) .* dhat) * U';
end

%% Next Z
Z = cell(K, 1);
for k = 1: K
    Z{k} = softThreshold(Z_hat{k} - Wz{k}, lambda1 / rho);
end

%% Next V
weight = [omega1, omega2, omega3] .* lambda2 ./ rho;
[vector1, index1] = altraAssistance(V_hat, p, K, weight); % altra 的辅助函数
[vector2, index2] = altraAssistance(Wv, p, K, weight); % altra 的辅助函数
assert(all(all(index1 == index2)), 'nextB: index1 != index2')
vector = vector1 - vector2;
Vvector = altra(vector, length(vector), index1, size(index1, 2)); % altra 的辅助函数
V = altraAssistance(Vvector, p, K);

end
