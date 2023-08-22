function weightVectors = updateWeightVectors(weightVectors, gbest)
    numSubproblems = size(weightVectors, 1); % 获取子问题的个数
    numObjectives = size(weightVectors, 2); % 获取目标函数的个数
    numGbest = size(gbest, 1); % 获取 Pareto 前沿中解的个数
    
    gbest(:,end - numObjectives + 1) = -gbest(:,end - numObjectives + 1);
    % 计算每个子问题与 Pareto 前沿中解之间的欧几里得距离
    distances = zeros(numSubproblems, numGbest);
    for i = 1:numSubproblems
        for j = 1:numGbest
            distances(i, j) = norm(weightVectors(i, :) - gbest(j, end - numObjectives + 1:end));
        end
    end
    
    % 对每个子问题，找到距离最近的 Pareto 前沿中解，并根据其目标函数值更新权重向量
    for i = 1:numSubproblems
        [~, nearestIndex] = min(distances(i, :)); % 找到距离最近的索引
        nearestSolution = gbest(nearestIndex, end - numObjectives + 1:end); % 找到距离最近的解
        
        % 对自适应参数 alpha 进行限制和调节，使其在一个较小的常数范围内变化
        alpha = rand * 0.1; % 随机生成一个自适应参数
        
        % 对 Pareto 前沿中解的目标函数值进行归一化处理，使其与权重向量在同一数量级
         minObjectives = min(gbest(:, end - numObjectives + 1:end)); % 计算每个目标函数值的最小值
         maxObjectives = max(gbest(:, end - numObjectives + 1:end)); % 计算每个目标函数值的最大值
         nearestSolution = (nearestSolution - minObjectives) ./ (maxObjectives - minObjectives); % 使用最大最小值归一化方法
        
        weightVectors(i, :) = weightVectors(i, :) + alpha * (nearestSolution - weightVectors(i, :)); % 更新权重向量
    end
    
end
