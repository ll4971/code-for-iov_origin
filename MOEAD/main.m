clc;
clear;
close all;

num_experiments = 1;
distanceGenerate(); % 获得距离

for times = 1:num_experiments
    rng(8);
    com = [30 33 36 46 49 50 54 57 59 59];
    x0 = [150 168 182 230 245 248 266 279 302 311];
    spc = [20 20 22 27 30 30 32 33 34 46];
    COM = [80 50 100 70 150 200 90];
    SPC = [40 60 80 60 100 140 100];
    N = [160 150 180 165 190 200 170];
    Ur = [0.1 0.6 0.1 0.1 0.6 0.1 0.1 1 0.6 0.1]*5;
    r = [0.7	0.2	0.3	0.8	0.1	0.8	0.7	0.3	1	0.7
        0.6	0.2	0.3	0.7	0.2	0.9	0.8	0.4	0.8	0.6
        0.6	0.3	0.5	0.6	0.1	0.8	0.8	0.3	0.9	0.6
        0.6	0.4	0.4	0.7	0.1	0.7	0.7	0.5	0.9	0.6
        0.7	0.4	0.2	0.8	0.1	0.8	0.9	0.4	0.9	0.5
        0.5	0.4	0.2	0.9	0.2	0.7	0.8	0.5	0.9	0.7
        0.6	0.2	0.3	0.6	0.3	0.9	0.7	0.3	0.8	0.6
        ]*10; 
    r = r';
    rho = 0.5;
    ka = 10;
    v = 0.5;
    epsilon = 1;
    sigma = 0.8;
    p = 1000;
    m = length(com);
    n = length(COM);

    % 算法参数
    NP = 80; %种群数量
    maxGenerations = 600;   %迭代次数
    subproblemSize = 8;
    numObjectives = 2; %目标函数个数
    numSubproblems = ceil(NP / subproblemSize);
    dim = m * 2;
    F = 0.8;
    CR = 0.2;

    
    % 初始化种群
    population = initPop(NP, numObjectives, dim, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p);    
    
    % 生成权重向量
    weightVectors = generateWeightVectors(numSubproblems,numObjectives);

    % 非支配排序（non-domination-sort）
    population = nonDominatedSort(population, numObjectives, dim);

    for gen = 1:maxGenerations
        for i = 1:numSubproblems
            % 根据权重向量和种群，选择邻居个体
            neighbors = selectNeighbors(i, weightVectors, population, subproblemSize); 
            % 交叉和变异操作
            newchrom = geneOperator(neighbors, numObjectives, dim, F, CR, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p); 
            % 执行子问题上的进化操作（交叉、变异等）
            offspring = evolveSubproblem(i, neighbors(:,1:dim+2), weightVectors, F, CR, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p);
            % 更新外部存档
%             population = updateExternalArchive(population, offspring);
            Nc = size(population,1);
            Nn = size(offspring,1);
            allchrom(1:Nc,:) = population;
            allchrom(Nc + 1 : Nc + Nn,1 : numObjectives+dim) = offspring;
            allchrom = nonDominatedSort(allchrom, numObjectives, dim);
            population = replace_chromosome(allchrom, numObjectives, dim, NP);
            FG1(gen,1) = -min(population(:,dim + 1));
            FG2(gen,1) = min(population(:,dim + 2));
            plot(-population(:,dim + 1) ,population(:,dim + 2),'*');
            str = sprintf('多目标遗传算法帕累托求解第%d次迭代',gen);
            title(str)
            xlabel('卖价与买价差值');
            ylabel('资源耗能和');
            pause(0.05)
            hold off
        end
        
        % 提取当前迭代的Pareto前沿
        gbest = population(:,1:dim+numObjectives);
        gbest = unique(gbest,'rows');
        % 更新权重向量
        weightVectors = updateWeightVectors(weightVectors, gbest);
    
        % 保存当前迭代的结果
%         saveResults(paretoFront, times, gen);
    end

end
