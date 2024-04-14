function est_hubcol = getNodeColumns(Vcell, r)
% Get node columns in estimated V
K = numel(Vcell);
est_hubcol = cell(size(Vcell));
for k = 1: K
    V = Vcell{k};
    V(logical(eye(size(V)))) = 0;
    est_hubcol{k} = find(sum(~~real(V), 1) > r);
end 
end
