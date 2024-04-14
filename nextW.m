function [Wtheta, Wz, Wv] = nextW(Theta, Z, V, Theta_hat, Z_hat, V_hat, ...
    Wtheta, Wz, Wv, K)
% W[t+1] <- argmin_W L{B[t+1], B_hat[t+1], W}

for k = 1: K
    Wtheta{k} = Wtheta{k} + Theta{k} - Theta_hat{k};
    Wz{k} = Wz{k} + Z{k} - Z_hat{k};
    Wv{k} = Wv{k} + V{k} - V_hat{k};
end
end
