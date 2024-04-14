function [Theta, Z, V, Theta_hat, Z_hat, V_hat, Wtheta, Wz, Wv] = init(p, K, initmode)
% Initializing B, B_hat and W
if nargin < 3 || isempty(initmode)
    initmode = 'zeros';
end

% 初始化为对角矩阵
Theta = generator(p, K, 'eye');
Z = generator(p, K, 'eye');
V = generator(p, K, 'eye');

% 初始化为零矩阵
Theta_hat = generator(p, K, initmode);
Z_hat = generator(p, K, initmode);
V_hat = generator(p, K, initmode);

% 初始化为零矩阵
Wtheta = generator(p, K, initmode);
Wz = generator(p, K, initmode);
Wv = generator(p, K, initmode);

end

function initcell = generator(p, K, initmode)
% 根据 initmode 初始化矩阵
initcell = cell(K, 1);
for k = 1: K
    if strcmp(initmode, 'eye')
        % 初始化为对角矩阵
        initcell{k} = eye(p);
    elseif strcmp(initmode, 'zeros')
        % 初始化为零矩阵
        initcell{k} = zeros(p);
    elseif strcmp(initmode, 'randn')
        % 初始化为服从正态分布的随机随机矩阵
        initcell{k} = randn(p);
    else
        error('init->generator: Wrong input initmode')
    end
end
end
