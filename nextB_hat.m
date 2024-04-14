function [Theta_hat, Z_hat, V_hat] = nextB_hat(Theta, Z, V, Wtheta, Wz, Wv, ...
    K, rho)
% B_hat[t+1] <- argmin_B_hat L{B[t+1], B_hat, W[t]}

Theta_hat = cell(K, 1);
Z_hat = cell(K, 1);
V_hat = cell(K, 1);

for k = 1: K
    %% Next Gamma
    Gamma = (rho / 6) * ((Theta{k} + Wtheta{k}) - (V{k} + Wv{k}) - (V{k} + Wv{k})' - (Z{k} + Wz{k}));
    
    %% Next Theta_hat
    Theta_hat{k} =  Theta{k} + Wtheta{k} - Gamma / rho;
    
    %% Next Z_hat
    Z_hat{k} =  Gamma / rho + Z{k} + Wz{k};

    %% Next V_hat
    V_hat{k} =  ((Gamma + Gamma') / rho) + V{k} + Wv{k};
end
end
