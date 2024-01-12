%% 设置多次实验
clc
clear
close all
num_experiments = 50;
%% 算法数据
NP = 100;          % 种群数量
maxgen_base = 100;     % 基础迭代次数
Pc = 0.8;
Pm = 0.2;
M = 2;            % 目标函数个数
rho = 0.5; % 最大收益目标函数买家花费占比
ka = 10;  
v = 0.5;    % 最大收益目标函数卖家收入占比
epsilon = 1;
sigma = 6;
p = 1000; % 传输功率
T_delay = zeros(num_experiments,1); % 总延迟时间
delay_average_results = zeros(num_experiments,1); % 平均延迟时间
x = zeros(num_experiments,1); % 车辆数横坐标
%% 主循环
for times = 1:num_experiments
    rng(times);
    %% 输入数据
    % 输入买家和卖家的数量
    m = 2*times; % 买家数量
    n = 3; % 卖家数量
    swt = 1; % 1：添加紧急程度；0：取消紧急程度
    rep = 0; % 1:存在信誉变化；0：不存在信誉变化
    % 调用 generate_data 函数生成需求和供给数据
    [com, spc, COM, SPC, Ur, r, N, D, x0] = generate_data(m, n, swt, rep);
    dim = m * 2 ;      % 决策变量维数
    %% 初始化种群
    [chrom, com, spc] = initpop(NP, M, dim, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 ,rho, v,ka, epsilon, sigma, p);
    
    %% 非支配排序（non-domination-sort）
    chrom = nonDominatedSort(chrom, M, dim );
    %% 记录延迟的数组
    maxgen = maxgen_base * m; %实际迭代次数
    delay = zeros(maxgen, 1);
    start_time = tic;
    elapsed_time = zeros(maxgen, 1);
    %% 进化过程
    figure(1);
    for gen = 1 : maxgen
        tic; % 记录迭代开始时间
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
        allchrom = chrom;
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
        toc
        delay(gen) = toc; % 记录迭代结束时间，并计算出延迟
        %% 计算并打印累积用时
        elapsed_time(gen) = toc(start_time);
        fprintf('累积用时：%.2f秒\n', toc(start_time));
    end

    %% 传输延迟 + 计算延迟
    y = chrom(1, m + 1: dim);
    for i = 1:m
        h = 1e-8 + exp(2-5*log10(D(i,y(i))));
        T_delay(times) = T_delay(times) + sum(x0 ./ com) + x0(i) / (spc(i) * log2(1 + p * h /sigma^2));
    end
    %% 平均延迟
    delay_average_results(times) = elapsed_time(maxgen)/maxgen;

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
    figure(4) 
    plot(elapsed_time) 
    xlabel('迭代次数') 
    ylabel('累积用时') 
    title('累积用时随迭代次数的变化')
    x(times) = m; 
    figure(5) 
    plot(x, T_delay)
    xlabel('车辆数') 
    ylabel('总延迟') 
    title('总延迟随车辆数的变化')
    %% 保存至excel
    if swt == 0
        file_path = '../NSGA-II_results'; % 修改为你希望保存的文件夹路径
    elseif swt == 1
        file_path = '../NSGA-II_results_Ur'; % 修改为你希望保存的文件夹路径
    end
    %保存帕累托解集
    file_name_01 = 'Pareto_results.xlsx'; % 修改为你希望保存的文件名
    gbest(:, dim+1) = -gbest(:, dim+1);   % 21列取反
    file_restore_01 = gbest(:,[dim+1:dim+2]);   %保存21和22列
    
    %保存总价值
    file_name_02 = 'Revenue_results.xlsx'; % 修改为你希望保存的文件名
    file_restore_02 = FG1;   %保存21和22列
    
    %保存总能耗
    file_name_03 = 'Consumption_results.xlsx'; % 修改为你希望保存的文件名
    file_restore_03 = FG2;   %保存21和22列

    %保存总延迟
    file_name_04 = 'T_delay_results.xlsx'; % 修改为你希望保存的文件名
    file_restore_04 = [x, T_delay];

    %保存算法平均延迟
    file_name_05 = 'delay_average_results.xlsx'; % 修改为你希望保存的文件名
    file_restore_05 = delay_average_results;
    
    %保存算法总延迟
    file_name_06 = 'elapsed_time_results.xlsx'; % 修改为你希望保存的文件名
    file_restore_06 = elapsed_time;



    % 使用 xlswrite 函数保存数据到 Excel 文件中
    xlswrite(fullfile(file_path, file_name_01), file_restore_01, times, 'A1'); % 将数据从 A1 单元格开始保存
    xlswrite(fullfile(file_path, file_name_02), file_restore_02, times, 'A1'); % 将数据从 A1 单元格开始保存
    xlswrite(fullfile(file_path, file_name_03), file_restore_03, times, 'A1'); % 将数据从 A1 单元格开始保存
    xlswrite(fullfile(file_path, file_name_04), file_restore_04, times, 'A1'); % 将数据从 A1 单元格开始保存
    xlswrite(fullfile(file_path, file_name_05), file_restore_05, times, 'A1'); % 将数据从 A1 单元格开始保存
    xlswrite(fullfile(file_path, file_name_06), file_restore_06, times, 'A1'); % 将数据从 A1 单元格开始保存
end

