clc;
clear;
close all;

num_experiments = 10;
distanceGenerate(); % 获得距离
swt = 0; % 1：添加紧急程度；0：取消紧急程度
com = [30 33 36 46 49 50 54 57 59 59];
x0 = [150 168 182 230 245 248 266 279 302 311];
spc = [20 20 22 27 30 30 32 33 34 46];
COM = [80 50 100 70 150 200 90];
SPC = [40 60 80 60 100 140 100];
N = [160 150 180 165 190 200 170];
if swt == 1
    Ur = [0.1 0.6 0.1 0.1 0.6 0.1 0.1 1 0.6 0.1]*5; %紧急程度
elseif swt ==0
    Ur = ones(1,10);
end
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
NP = 200; %种群数量
maxGenerations = 100;   %迭代次数
subproblemSize = 10;
numObjectives = 2; %目标函数个数
numSubproblems = ceil(NP / subproblemSize);
dim = m * 2;
F0 = 1;
CR0 = 0.5;
gbest = zeros(1,22);


for times = 1:num_experiments
    rng(times);
    % 初始化种群
    population = initPop(NP, numObjectives, dim, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p);    
    
    % 生成权重向量
    weightVectors = generateWeightVectors(numSubproblems,numObjectives);

    % 非支配排序（non-domination-sort）
    population = nonDominatedSort(population, numObjectives, dim);

    for gen = 1:maxGenerations
        times
        gen
        for i = 1:numSubproblems
            lamda = exp(1-maxGenerations/(maxGenerations+1-gen));
            F = F0*2^lamda;
            CR = CR0*(1-lamda);
            % 根据权重向量和种群，选择邻居个体
            neighbors = selectNeighbors(i, weightVectors, population, subproblemSize); 
            % 交叉和变异操作
%             offspring = geneOperator(neighbors, numObjectives, dim, F, CR, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p); 
            % 执行子问题上的进化操作（交叉、变异等）
            offspring = evolveSubproblem(gbest, i, neighbors(:,1:dim+2), weightVectors, F, CR, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p);
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
            figure(1);
            plot(-population(:,dim + 1) ,population(:,dim + 2),'*');
            str = sprintf('多目标遗传算法帕累托求解第%d次迭代',gen);
            title(str)
            xlabel('卖价与买价差值');
            ylabel('资源耗能和');
            pause(0.05)

            % 提取当前迭代的Pareto前沿

        end
            gbest = population(:,1:dim+numObjectives);
            gbest = unique(gbest,'rows');
            % 更新权重向量
            weightVectors = updateWeightVectors(weightVectors, gbest);
        % 保存当前迭代的结果
%         saveResults(paretoFront, times, gen);

    end
        figure(1)
        plot(-population(:,dim + 1),population(:,dim + 2),'ko')
        xlabel('市场总价值')
        ylabel('资源耗能和')
        grid on
        title('帕累托解集')
        figure(2)
        plot(FG1,'k-')
        xlabel('迭代次数')
        ylabel('市场总价值')
        grid on
        
        figure(3)
        plot(FG2,'k-')
        xlabel('迭代次数')
        ylabel('资源耗能和')
        grid on
        fprintf('帕累托解集数量为 %d\n',size(gbest,1))

        if swt == 0
            file_path = '../MOEAD_results'; % 修改为你希望保存的文件夹路径
        elseif swt == 1
            file_path = '../MOEAD_results_Ur'; % 修改为你希望保存的文件夹路径
        end
        %保存帕累托解集
        file_name_01 = 'Pareto_results.xlsx'; % 修改为你希望保存的文件名
        gbest(:, 21) = -gbest(:, 21);   % 21列取反
        file_restore_01 = gbest(:,[21:22]);   %保存21和22列
        
        %保存总价值
        file_name_02 = 'Revenue_results.xlsx'; % 修改为你希望保存的文件名
        file_restore_02 = FG1;   %保存21和22列
        
        %保存总能耗
        file_name_03 = 'Consumption_results.xlsx'; % 修改为你希望保存的文件名
        file_restore_03 = FG2;   %保存21和22列

        % 使用 xlswrite 函数保存数据到 Excel 文件中
        xlswrite(fullfile(file_path, file_name_01), file_restore_01, times, 'A1'); % 将数据从 A1 单元格开始保存
        xlswrite(fullfile(file_path, file_name_02), file_restore_02, times, 'A1'); % 将数据从 A1 单元格开始保存
        xlswrite(fullfile(file_path, file_name_03), file_restore_03, times, 'A1'); % 将数据从 A1 单元格开始保存
end
