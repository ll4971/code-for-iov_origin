function offspring = evolveSubproblem(gbest, i, neighbors, weightVectors, F, CR, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p)
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
        if gbest(:) ~= 0 % 检查变量 gbest 是否为0
            c = gbest(1,:);
        end
        
        trial = zeros(1, dim+2); % 初始化试验个体
%       trial_vars = gendifEvo(trial, parent(1:dim), a(1:dim), b(1:dim), c(1:dim), F, CR);% 差分进化函数
        trial_vars = gendifEvo(parent(1:dim), F, CR, c(1:dim), dim, n);% 差分进化函数
        % 进行四舍五入或取整操作，将小数值转换为整数，并限制在[1,7]之间

        % 将变异和交叉后得到的决策变量赋值给试验个体
        trial(1:dim) = trial_vars;
        
        % 计算试验个体的目标函数值
        trial_objectives = fitness(trial_vars, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p);
        trial(dim+1:end) = trial_objectives;
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

% 定义差分进化函数，主要实现归一化
%     function trial_vars = gendifEvo(trial, parent, a, b, c, F, CR)
%         normalizedParent = normalizeMatrix(parent);
%         normalizedA = normalizeMatrix(a);
%         normalizedB = normalizeMatrix(b);
%         normalizedC = normalizeMatrix(c);
%         % 执行差分进化的变异操作
%         trial = normalizedParent + F * (normalizedA - normalizedB) + F * (normalizedC - normalizedParent);
%         
%         % 执行差分进化的交叉操作
%         mask = rand < CR;
%         trial_vars = mask .* trial+ (1 - mask) .* normalizedParent;
%         trial_vars = round(7*trial_vars);
%         trial_vars = max(min(trial_vars, 7), 1);
%     end
    function trial_vars = gendifEvo(parent, Pm, CR, c, dim, n)
        % 执行差分进化的变异操作
        x = parent;
        % 执行差分进化的交叉操作
        mask = rand < CR;
        trial_vars = mask .* c+ (1 - mask) .* parent;
        if rand < Pm
            r0 = randi(dim);
            x(r0) = randi(n);
        end
        trial_vars = x;
    end
    function normalized_matrix = normalizeMatrix(matrix)
        min_values = min(matrix);
        max_values = max(matrix);
        normalized_matrix = (matrix - min_values) ./ (max_values - min_values);
    end

    function fitnessValue = fitnessFunction(objectives, weightVector)
        numObjectives = length(objectives); % 获取目标函数的个数
        fitnessValue = 0; % 初始化适应度函数值为0
        for ii = 1:numObjectives % 对每个目标函数进行加权求和
            fitnessValue = fitnessValue + objectives(ii) * weightVector(ii);
        end
    end
