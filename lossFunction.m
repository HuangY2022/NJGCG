function lsfnval = lossFunction(S, Theta, Z, V, ...
    nlist, p, K, lambda1, lambda2, omega1, omega2, omega3)
% Calculating loss value of objective function

%% Loss function
loss = 0;
for k = 1: K
    loss = loss + nlist(k) * (-log(det(Theta{k})) + trace(S{k} * Theta{k}));
end 

%% Regularization of Z
R1 = 0;
for k = 1: K
    R1 = R1 + lambda1 * sum(sum(abs(Z{k})));
end

%% Regularization of V
R2 = 0;
R3 = 0;
R4 = 0;
for k = 1: K
    R2 = R2 + lambda2 * omega1 * sum(sum(abs(V{k})));
    R3 = R3 + lambda2 * omega2 * sum(sqrt(sum(V{k}.^2)));
end

for j = 1: p
    R4sub = 0;
    for k = 1: K
        R4sub = R4sub + sum(V{k}(:, j).^2);
    end
    R4 = R4 + lambda2 * omega3 * sqrt(R4sub);
end

%% objective function
lsfnval = [loss + R1 + R2 + R3 + R4, loss, R1, R2, R3, R4];

end
