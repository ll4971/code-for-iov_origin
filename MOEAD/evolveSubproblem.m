function offspring = evolveSubproblem(i, neighbors, weightVectors, F, CR, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p)
    NP = size(neighbors, 1); % 获取邻居子问题的个数
    dim = size(neighbors, 2) - 2; % 获取决策变量的个数
    
    offspring = zeros(NP, dim + 2); % 初始化子代种群
    
    for j = 1:NP % 对每个邻居子问题进行循环
        parent = neighbors(j, :); % 获取父代个体
        
        % 从邻居中随机选择三个不同的个体
        rand_indices = randperm(NP, 3);
        a = neighbors(rand_indices(1), :);
        b = neighbors(rand_indices(2), :);
        c = neighbors(rand_indices(3), :);

        trial = zeros(1, dim + 2); % 初始化试验个体
        
        % 执行差分进化的变异操作
        trial(1:dim) = parent(1:dim) + F * (a(1:dim) - b(1:dim)) + F * (c(1:dim) - parent(1:dim));
        
        % 执行差分进化的交叉操作
        mask = rand(1,dim) < CR;
        trial_vars = mask .* trial(1:dim) + (1 - mask) .* parent(1:dim);
        
        % 进行四舍五入或取整操作，将小数值转换为整数，并限制在[1,7]之间
        trial_vars = round(trial_vars);
        trial_vars = max(min(trial_vars, 7), 1);
        
        % 将变异和交叉后得到的决策变量赋值给试验个体
        trial(1:dim) = trial_vars;
        
        % 计算试验个体的目标函数值
        trial_objectives = fitness(trial_vars, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p);
        
        % 计算试验个体的约束函数值
        gfun = 0;
        cs = zeros(1, n);  % 卖家实际供应量
        for k = 1 : n
            cs(k) = sum(com(trial_vars(1:m) == k));
            gfun = gfun + max(cs(k)-COM(k),0);
        end

        ss = zeros(1, n);  % 卖家实际供应量
        for k = 1 : n
            ss(k) = sum(spc(trial_vars(m+1:m*2) == k));
            gfun = gfun + max(ss(k)-SPC(k),0);
        end
        
        % 将目标函数值和约束函数值赋值给试验个体
        trial(dim + 1:dim + 2) = trial_objectives + gfun * 1e8;
        
        % 如果试验个体支配父代个体，进行替换
        if dominates(trial_objectives, parent(dim + 1:end))
            offspring(j, :) = trial;
        else
            % 如果试验个体不支配父代个体，根据权重向量和适应度函数选择更优的一个保留
            weightVector = weightVectors(i, :); % 获取当前子问题的权重向量
            fitnessParent = fitnessFunction(parent(dim + 1:end), weightVector); % 计算父代个体的适应度函数值
            fitnessTrial = fitnessFunction(trial_objectives, weightVector); % 计算试验个体的适应度函数值
            if fitnessTrial < fitnessParent % 如果试验个体的适应度函数值更小，说明更优
                offspring(j, :) = trial;
            else
                offspring(j, :) = parent;
            end
        end
    end

% 定义适应度函数，根据权重向量和目标函数值计算一个标量值，用于评价个体的优劣
function fitnessValue = fitnessFunction(objectives, weightVector)
    numObjectives = length(objectives); % 获取目标函数的个数
    fitnessValue = 0; % 初始化适应度函数值为0
    for ii = 1:numObjectives % 对每个目标函数进行加权求和
        fitnessValue = fitnessValue + objectives(ii) * weightVector(ii);
    end
end

end
