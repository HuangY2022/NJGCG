function st = softThreshold(x, thres) 
% 软阈值函数
    st = sign(x).*max(abs(x) - thres, 0);
end
