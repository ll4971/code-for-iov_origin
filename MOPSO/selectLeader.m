function selected = selectLeader(rep) 
% 轮盘赌
Pcum = cumsum(rep.quality(:,2));     % 轮盘累积概率
sel_hyp = rep.quality(find(rand(1,1)*max(Pcum)<=Pcum,1,'first'),1); % 选择超立方体
% 在被选中的超立方体中随机选择一个粒子作为Leader
idx = 1:1:length(rep.grid_idx);
selected = idx(rep.grid_idx==sel_hyp);
selected = selected(randi(length(selected)));
end

