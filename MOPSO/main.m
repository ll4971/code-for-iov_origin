clc
clear
close all
%% 输入数据
com = [30	33	36	46	49	50	54	57	59	59];   % 买家需求计算资源
spc = [20	20	22	27	30	30	32	33	34	46];   % 买家需求频谱资源
COM = [80	50	100	70	150	200	90];   % 卖家供给计算资源
SPC = [40	60	80	60	100	140	100];  % 卖家供给传输资源
N = [160	150	180	165	190	200	170];  % 卖价
r = [0.7	0.2	0.3	0.8	0.1	0.8	0.7	0.3	1	0.7
    0.6	0.2	0.3	0.7	0.2	0.9	0.8	0.4	0.8	0.6
    0.6	0.3	0.5	0.6	0.1	0.8	0.8	0.3	0.9	0.6
    0.6	0.4	0.4	0.7	0.1	0.7	0.7	0.5	0.9	0.6
    0.7	0.4	0.2	0.8	0.1	0.8	0.9	0.4	0.9	0.5
    0.5	0.4	0.2	0.9	0.2	0.7	0.8	0.5	0.9	0.7
    0.6	0.2	0.3	0.6	0.3	0.9	0.7	0.3	0.8	0.6
    ] * 10;
r = r';
D = [15	18	25	27	30	17	10	18	7	5
    20	13	17	19	5	7	17	30	13	12
    12	30	12	24	17	19	23	2	7	17
    8	20	20	5	14	15	29	10	18	7
    23	8	15	21	27	29	16	14	14	8
    11	25	18	16	28	9	23	23	11	17
    12	12	16	23	28	14	17	14	77	10
    ];
D = D';
x0 = [150  160 185 220	222	230	250	275	290	300];
rho = 0.5;
ka = 10;
v = 0.5;
epsilon = 1;
sigma = 1;
p = 0.5;
m = length(com);  % 买家数量
n = length(COM);  % 卖家数量
dim = m * 2;
% MOPSO
%% 粒子群参数
NP = 80;                   % 种群数量
M = 2;
maxgen = 600;               % 最大迭代次数
w = 0.4;                    % 惯性系数
c1 = 2;
c2 = 2;
ngrid = 20;
Xmin = ones(1, dim);
Xmax = n * ones(1, dim);
Vmax = (Xmax - Xmin) / 2;
Vmin = - Vmax;
%% 初始化
X = zeros(NP, dim);
V = zeros(NP,dim);
fx = zeros(NP, M);
for i = 1:NP
    X(i,:) = initpop(m, n, com, spc, COM, SPC);
    V(i,:) = Vmin + (Vmax - Vmin) .* rand(1,dim);
    fx(i,:) = fitness(X(i,:), m, n, com, spc, COM, SPC, N, r ,D, x0 , rho, v, ka, epsilon, sigma, p);
end

pbest = X;
fpbest = fx;
Idx = getNondominated(fx);
rep.X  = X(Idx,:);         % 非劣解集
rep.fx = fx(Idx,:);
rep = updateGrid(rep,ngrid);
figure(1)
for gen = 1:maxgen
    % 选择leader
    h = selectLeader(rep);
    gbest = rep.X(h,:);
    for i = 1 : NP
        % 粒子速度与位置更新
        V(i,:) = w .* V(i,:) + c1 * rand(1,dim) .* (pbest(i,:) - X(i,:)) + c2 * rand(1,dim) .* (gbest - X(i,:));
        % 保证粒子速度位于界内
        index = (V(i,:) > Vmax);
        V(i,index) = Vmax(index);
        index = (V(i,:) < Vmin);
        V(i,index) = Vmin(index);
        
        % 种群更新
        X(i,:) = X(i,:)+V(i,:);
        % 保证粒子位置在界内
        index = (X(i,:) > Xmax);
        X(i,index) = Xmax(index);
        index = (X(i,:) < Xmin);
        X(i,index) = Xmin(index);
        
        % 评价适应度
        fx(i,:) = fitness(X(i,:), m, n, com, spc, COM, SPC, N, r ,D, x0 , rho, v,ka, epsilon, sigma, p);
    end
    % 更新repository
    rep = updateRepository(rep,X,fx,ngrid);
    if(size(rep.X,1) > NP)
        rep = deleteFromRepository(rep,size(rep.X,1)-NP,ngrid);
    end
    
    % 更新个体最优pbest
    if dominates(fx(i,:), fpbest(i,:))
        fpbest(i,:) = fx(i,:);
        pbest(i,:) = X(i,:);
    elseif dominates(fpbest(i,:), fx(i,:))
        
    else
        if rand < 0.5
            fpbest(i,:) = fx(i,:);
            pbest(i,:) = X(i,:);
        end
    end
    FG1(gen,1) = -min(rep.fx(:,1));
    FG2(gen,1) = min(rep.fx(:,2));
    % 绘图显示
    plot(-rep.fx(:,1),rep.fx(:,2),'ok'); hold on;
    hold off
    grid on; xlabel('f1'); ylabel('f2');
    pause(0.01)
    drawnow;
    axis square;
end
gbest = rep.X;
fgbest = rep.fx;
[fgbest,IA] = unique(fgbest,'rows');
gbest = round(gbest(IA,:));
fgbest(:,1) = -fgbest(:,1);
figure(1)
plot(fgbest(:,1),fgbest(:,2),'ko')
xlabel('Overall Revenue')
ylabel('Resource Consumption')
grid on
title('Pareto solution set')
figure(2)
plot(FG1,'k-')
xlabel('Iterations Times')
ylabel('Overall Revenue')
grid on

figure(3)
plot(FG2,'k-')
xlabel('Iterations Times')
ylabel('Resource consumption')
grid on
fprintf('帕累托解集数量为 %d\n',size(gbest,1))