% 更新保留集 
function rep = updateRepository(rep, X, fx, ngrid)
% 确定非支配解
Idx  = getNondominated(fx);
rep.X    = [rep.X; X(Idx,:)];
rep.fx = [rep.fx; fx(Idx,:)];
% 合并后的非支配解
Idx  = getNondominated(rep.fx);
rep.fx= rep.fx(Idx,:);
rep.X    = rep.X(Idx,:);
% Updating the grid
rep        = updateGrid(rep,ngrid);
end
