%% 设置多次实验
clc
clear
close all
num_experiments = 1;
%% 输入数据
% 输入买家和卖家的数量
m = 50; % 买家数量
n = 7; % 卖家数量
swt = 1; % 1：添加紧急程度；0：取消紧急程度
rep = 0; % 1:存在信誉变化；0：不存在信誉变化
% 调用 generate_data 函数生成需求和供给数据
[com, spc, COM, SPC, Ur, r, N, D, x0] = generate_data(m, n, swt, rep);


    rho = 0.5; % 最大收益目标函数买家花费占比
    ka = 10;   %
    v = 0.5;    % 最大收益目标函数卖家收入占比
    epsilon = 1;
    sigma = 6;
    p = 1000; % 传输功率
    %% 算法数据
    NP = 80;          % 种群数量
    maxgen = 60;     % 迭代次数
    Pc = 0.8;
    Pm = 0.2;
    M = 2;            % 目标函数个数
    dim = m * 2 ;      % 决策变量维数
%% 主循环
for times = 1:num_experiments
    rng(8);

    %% 记录延迟的数组
    delay = zeros(maxgen, 1);
    start_time = tic;
    elapsed_time = zeros(maxgen, 1);
    figure(1);
    for gen = 1 : maxgen
        tic; % 记录迭代开始时间
        times
        gen
        %% 初始化种群
        chrom = initpop(NP, M, dim, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 ,rho, v,ka, epsilon, sigma, p);
    
        %% 非支配排序（non-domination-sort）
        chrom = nonDominatedSort(chrom, M, dim );
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

