%% 设置多次实验
clc
clear
close all
num_experiments = 1;
distanceGenerate(); % 获得距离
for times = 1:num_experiments
    rng(8);
    %% 输入数据
    com = [30	33	36	46	49	50	54	57	59	59];   % 买家需求计算资源
    x0 = [150 168 182 230 245 248 266 279 302 311]; % 买家需求任务量
    spc = [20	20	22	27	30	30	32	33	34	46];   % 买家需求频谱资源
    COM = [80	50	100	70	150	200	90];   % 卖家供给计算资源
    SPC = [40	60	80	60	100	140	100];  % 卖家供给传输资源
    N = [160	150	180	165	190	200	170];  % 卖价
    
    swt = 1; % 1：添加紧急程度；0：取消紧急程度
    if swt == 1
        Ur = [0.1 0.6 0.1 0.1 0.6 0.1 0.1 1 0.6 0.1]*5; %紧急程度
    elseif swt ==0
        Ur = ones(1,10);
    end
    
    %信誉值
    %行代表V，列代表RSU，[1,2]代表V1对RSU2
    % r = (0.5 + 0.5 * randn(10, 7));
    
    rep = 1; % 1:存在信誉变化；0：不存在信誉变化
    if rep == 0
        r = [0.7	0.2	0.3	0.8	0.1	0.8	0.7	0.3	1	0.7
        0.6	0.2	0.3	0.7	0.2	0.9	0.8	0.4	0.8	0.6
        0.6	0.3	0.5	0.6	0.1	0.8	0.8	0.3	0.9	0.6
        0.6	0.4	0.4	0.7	0.1	0.7	0.7	0.5	0.9	0.6
        0.7	0.4	0.2	0.8	0.1	0.8	0.9	0.4	0.9	0.5
        0.5	0.4	0.2	0.9	0.2	0.7	0.8	0.5	0.9	0.7
        0.6	0.2	0.3	0.6	0.3	0.9	0.7	0.3	0.8	0.6
        ]*10; 
        r = r';
    elseif rep == 1
        r = [0.7	0.2	0.3	0.8	0.1	0.8	0.7	0.3	0.1	0.7
        0.6	0.2	0.3	0.7	0.2	0.9	0.8	0.4	0.1	0.6
        0.6	0.3	0.5	0.6	0.1	0.8	0.8	0.3	0.1	0.6
        0.6	0.4	0.4	0.7	0.1	0.7	0.7	0.5	0.1	0.6
        0.7	0.4	0.2	0.8	0.1	0.8	0.9	0.4	0.1	0.5
        0.5	0.4	0.2	0.9	0.2	0.7	0.8	0.5	0.1	0.7
        0.6	0.2	0.3	0.6	0.3	0.9	0.7	0.3	0.1	0.6
        ]*10; 
        r = r';
    end
    % D = [15	18	25	27	30	17	10	18	7	5
    %     20	13	17	19	5	7	17	30	13	12
    %     12	30	12	24	17	19	23	2	7	17
    %     8	20	20	5	14	15	29	10	18	7
    %     23	8	15	21	27	29	16	14	14	8
    %     11	25	18	16	28	9	23	23	11	17
    %     12	12	16	23	28	14	17	14	77	10
    %     ]; 

    rho = 0.5; % 最大收益目标函数买家花费占比
    ka = 10;   %
    v = 0.5;    % 最大收益目标函数卖家收入占比
    epsilon = 1;
    sigma = 0.8;
    p = 1000; % 传输功率
    m = length(com);  % 买家数量
    n = length(COM);  % 卖家数量
    %% 算法数据
    NP = 80;          % 种群数量
    maxgen = 600;     % 迭代次数
    Pc = 0.8;
    Pm = 0.2;
    M = 2;            % 目标函数个数
    dim = m * 2 ;      % 决策变量维数
    %% 初始化种群
    chrom = initpop(NP, M, dim, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 ,rho, v,ka, epsilon, sigma, p);
    
    %% 非支配排序（non-domination-sort）
    chrom = nonDominatedSort(chrom, M, dim );
    
    %% 进化过程
    figure(1);
    for gen = 1 : maxgen
        times
        gen
        % 选择父代用于繁殖后代。原始NSGA-II采用基于拥挤度比较算子的二进制锦标赛进行选择的
        pool = round(NP/2);
        tour = 2;
        % 选择操作    
        parentchrom = tournamentSelect(chrom, pool, tour);   
        % 交叉和变异操作
        newchrom = geneOperator(parentchrom, M, dim, Pm, Pc, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p);   
        % 子代个体和父代个体融合
        Nc = size(chrom,1);
        Nn = size(newchrom,1);
        allchrom(1:Nc,:) = chrom;
        allchrom(Nc + 1 : Nc + Nn,1 : M+dim) = newchrom;
      
        % 融合种群非支配排序
        allchrom = nonDominatedSort(allchrom, M, dim);
        chrom = replace_chromosome(allchrom, M, dim, NP);
        FG1(gen,1) = -min(chrom(:,dim + 1));
        FG2(gen,1) = min(chrom(:,dim + 2));
        plot(-chrom(:,dim + 1) ,chrom(:,dim + 2),'*');
        str = sprintf('多目标遗传算法帕累托求解第%d次迭代',gen);
        title(str)
        xlabel('卖价与买价差值');
        ylabel('资源耗能和');
        pause(0.05)
        hold off
    end
    %% 结果输出
    clc
    gbest = chrom(:,1:dim+M);
    gbest = unique(gbest,'rows');
    
    figure(1)
    plot(-chrom(:,dim + 1),chrom(:,dim + 2),'ko')
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
    
    
    %%保存至excel
    if swt == 0
        file_path = '../NAGA-II_results'; % 修改为你希望保存的文件夹路径
    elseif swt == 1
        file_path = '../NAGA-II_results_Ur'; % 修改为你希望保存的文件夹路径
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
%     xlswrite(fullfile(file_path, file_name_01), file_restore_01, times, 'A1'); % 将数据从 A1 单元格开始保存
%     xlswrite(fullfile(file_path, file_name_02), file_restore_02, times, 'A1'); % 将数据从 A1 单元格开始保存
%     xlswrite(fullfile(file_path, file_name_03), file_restore_03, times, 'A1'); % 将数据从 A1 单元格开始保存
end

