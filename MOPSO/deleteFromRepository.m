% 从保留集中删除粒子 
function rep = deleteFromRepository(rep, Ne, ngrid)
[n,M] = size(rep.fx);
fx = rep.fx;
CD = zeros(n, M);  % 每个个体各个目标函数下的拥挤度
for m = 1 : M
    [~, objfunsortidx] = sort(fx(:, m));
    fxs = fx(objfunsortidx,:);
    objfunsort_min = fxs(1, m);
    bojfunsort_max = fxs(n, m);
    CD(objfunsortidx(1), m) = inf;
    CD(objfunsortidx(n), m) = inf;
    
    for j = 2 : n - 1
        objfun_next  = fxs(j + 1, m);
        objfun_last  = fxs(j - 1, m);
        if (bojfunsort_max - objfunsort_min == 0)
            CD(objfunsortidx(j), m) = Inf;
        else
            CD(objfunsortidx(j), m) =  (objfun_next - objfun_last)/(bojfunsort_max - objfunsort_min);
        end
    end
end
CD = sum(CD,2);
% 删除多余的粒子
[~,del_idx] = sort(CD,'ascend');
del_idx = del_idx(1:Ne);
rep.X(del_idx,:) = [];
rep.fx(del_idx,:) = [];
rep = updateGrid(rep,ngrid);
end
