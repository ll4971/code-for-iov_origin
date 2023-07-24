 % 更新超立方体网格，每一个立方体网格的质量由其内部粒子数目决定
function rep = updateGrid(rep,ngrid)
% 每一个超立方体的界限
M = size(rep.fx,2);
rep.hypercube_limits = zeros(ngrid+1,M);
for j = 1:1:M
    rep.hypercube_limits(:,j) = linspace(min(rep.fx(:,j)),max(rep.fx(:,j)),ngrid+1)';
end

% 计算粒子所属网格
Nr = size(rep.fx,1);
rep.grid_idx = zeros(Nr,1);
rep.grid_subidx = zeros(Nr,M);
for i = 1:1:Nr
    idnames = [];
    for k = 1:1:M
        rep.grid_subidx(i,k) = find(rep.fx(i,k)<=rep.hypercube_limits(:,k)',1,'first')-1;
        if(rep.grid_subidx(i,k)==0)
            rep.grid_subidx(i,k) = 1;
        end
        idnames = [idnames ',' num2str(rep.grid_subidx(i,k))];
    end
    rep.grid_idx(i) = eval(['sub2ind(ngrid.*ones(1,Nr)' idnames ');']);
end

rep.quality = zeros(ngrid,2);
ids = unique(rep.grid_idx);
for i = 1:length(ids)
    rep.quality(i,1) = ids(i);  % First, the hypercube's identifier
    rep.quality(i,2) = 10/sum(rep.grid_idx==ids(i)); % Next, its quality
end
end
