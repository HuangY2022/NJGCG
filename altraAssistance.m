function [outobj, index] = altraAssistance(anyobj, p, K, weight)
% Creat vector and index for function altra or let vector return to cell
if nargin == 4
    % 使用 altra 函数前调用
    [outobj, index] = altraBegin(anyobj, p, K, weight);
elseif nargin == 3
    % 使用 altra 函数后调用
    outobj = altraEnd(anyobj, p, K);
else
    error('altraAssistance: The number of input arguments is wrong')
end
end

function [vector, index] = altraBegin(anycell, p, K, weight)
% 将矩阵转成向量，并按照 altra 的规则设置群组稀疏的超参矩阵
%% Creat vector
temp = cell2mat(anycell);
vector = temp(:);
pK = p * K;
ppK = p^2 * K;
assert(all(size(vector) == [ppK, 1]), 'altraReady: Dimensions of input cell is wrong')

%% Creat index
% 三行的矩阵，第一行为某个群组的起始元素在向量中的索引位置，
% 第二行为某个群组的终结元素在向量中的索引位置，
% 第三行为该群组的稀疏惩罚超参数
index1 = [-1, -1, weight(1)]';
index2 = [1: p: ppK;  p: p: ppK; ones(1, pK) * weight(2)];
index3 = [1: pK: ppK; pK: pK: ppK; ones(1, p) * weight(3)];
index = [index1, index2, index3];

end

function [anycell] = altraEnd(vector, p, K)
% 将向量还原为矩阵元胞数组
pK = p * K;
ppK = p^2 * K;
assert(all(size(vector) == [ppK, 1]), 'altraEnd: Dimensions of input vector is wrong')
temp = reshape(vector, pK, p);
anycell = mat2cell(temp, ones(K, 1) .* p, p);
end
