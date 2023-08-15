% 定义一个名为 generateWeightVectors 的函数
% 输入参数为 numSubproblems 和 numObjectives
% 输出参数为 weightVectors
function weightVectors = generateWeightVectors(numSubproblems, numObjectives)
    % 初始化一个 numSubproblems 行 numObjectives 列的零矩阵
    weightVectors = zeros(numSubproblems, numObjectives);
    
    % 用一个 for 循环遍历每一行
    for i = 1:numSubproblems
        % 生成一个 1 行 numObjectives 列的随机向量 w
        w = rand(1, numObjectives);
        % 把 w 除以它的元素之和，使得 w 的元素之和为 1
        w = w / sum(w); % Normalize to ensure the weights sum up to 1
        % 把 w 赋值给 weightVectors 的第 i 行
        weightVectors(i, :) = w;
    end
end



