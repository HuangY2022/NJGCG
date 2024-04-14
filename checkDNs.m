function checkDNs(anycell, p, K, funcname)
% Check dimension of input obj; c is for cell, d is for dict
if nargin < 4 || isempty(funcname)
    funcname = 'no';
end

if isequal(class(anycell), 'cell')
    if ~(size(anycell, 1) == K)
        disp([funcname, ': The k dimension is wrong'])
        disp(size(anycell))
    end
    for k = 1: numel(anycell)
        anymatrix = anycell{k};
        checkOne(anymatrix, p, k, funcname)
    end
else
    checkOne(anycell, p, 0, funcname)
end

end

function checkOne(anymatrix, p, k, funcname)
if ~(all(size(anymatrix) == [p, p]))
    disp([funcname, ': The i or j dimension is wrong --', num2str(k)])
    disp(size(anymatrix))
end
if ~(sum(sum(isnan(anymatrix))) == 0)
    disp([funcname, ': The matrix have nan --', num2str(k)])
    [nanhs, nanls] = find(isnan(anymatrix));
end
if ~(sum(sum(isinf(anymatrix))) == 0)
    disp([funcname, ': The matrix have inf --', num2str(k)])
    [infhs, infls] = find(isinf(anymatrix));
end
end
