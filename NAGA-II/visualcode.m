clc
clear
close all

distanceGenerate(); % 获得距离
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
    
    if rep == 0 && swt == 1
        ch = [1     6     3     6     7     5     5     6     6     3     6     5     2     1     6     3     2     6     3     7];
    elseif rep == 0 && swt == 0
        ch = [5     6     3     6     7     5     5     4     6     6     1     4     2     6     6     3     5     6     3     6];
    elseif rep == 1 && swt == 1
        ch = [1	6	3	6	7	5	5	6	3	6	2	7	2	1	6	3	5	6	3	7];
    elseif rep == 1 && swt == 0
        ch = [5	6	3	6	2	5	5	6	7	6	1	5	2	6	6	3	5	6	3	6];
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




x = ch(1:m);
y = ch(m+1:m*2);
lambda = zeros(m, n); %计算资源购买
for i = 1 : m
   lambda(i,x(i)) = 1; 
end
u = zeros(m, n); %频谱资源购买
for i = 1 : m
   u(i,y(i)) = 1; 
end
C = zeros(1, m);
for i = 1 : m
   C(i) = 0;
   for j = 1 : n
      C(i) = C(i) +  N(j) / r(i,j) /Ur(i)* (lambda(i,j) * com(i) + u(i,j) * spc(i)); % 第i辆V购买花费
   end
end

R = zeros(1, n);
for j = 1 : n
   R(j) = 0;
   rRSU(j)=0;
   for i = 1 : m
      R(j) = R(j) +  (N(j) * r(i,j)* Ur(i) - 0.5) * (lambda(i,j) * com(i) + u(i,j) * spc(i)); % 第j个RSU出售收益
      rRSU(j) = r(i,j)+ rRSU(j);
   end
end

u0 = zeros(1, n);
cs = zeros(1, n);
ss = zeros(1, n);
% 卖家实际供应量
for k = 1 : n
    cs(k) = sum(com(x == k));
    ss(k) = sum(spc(y == k));
    u0(k)=rRSU(k)*(COM(k)-cs(k)+ SPC(k)-ss(k))+u0(k);   
end
U = ka *(-rho * sum(C) + v * sum(R))+ 1/ka *sum(u0) ;

%能耗计算
E = epsilon * sum(com .* x0);
for i = 1 : m
    h = 1e-8 + exp(2-5*log10(D(i,y(i))));
    E = E + p * x0(i) / (spc(i) * log2(1 + p * h /sigma^2));
end

%% 约束处理
gfun = 0;
cs = zeros(1, n);  % 卖家实际供应量
for k = 1 : n
    cs(k) = sum(com(x == k));
    gfun = gfun + max(cs(k)-COM(k),0);
end


ss = zeros(1, n);  % 卖家实际供应量
for k = 1 : n
    ss(k) = sum(spc(y == k));
    gfun = gfun + max(ss(k)-SPC(k),0);
end


disp('gbest=1时，计算量分配：')
for j = 1 : m
    fprintf('买家 %d —— 卖家 %d \n',j, x(j))
end
disp('频谱量分配：')
for j = 1 : m
    fprintf('买家 %d —— 卖家 %d \n',j, y(j))
end
fprintf('C = %.3f\n', sum(C))
fprintf('R = %.3f\n', sum(R))
fprintf('U = %.3f\n', U)
fprintf('E = %.3f\n', E)
% C = C';
% R = R';
% disp('紧急程度最高买家(8)对应的卖家资源分配情况：')
buyer = 8; %买家8
seller = 6; %卖家6

% disp('计算资源：')
% x_all = sum(com(x == seller)); % 卖家t实际出售计算资源量
% x_percentage = x_all/(COM(seller));
% xt_percentage = com(buyer)/(COM(seller));
% fprintf('卖家的出售资源比例为：%.3f\n',x_percentage)
% fprintf('紧急程度最高买家(8)占卖家的出售资源比例为：%.3f\n',xt_percentage)
% 
% disp('频谱资源：')
% y_all = sum(spc(y == seller)); % 卖家t实际出售计算资源量
% y_percentage = y_all/(SPC(seller));
% yt_percentage = spc(buyer)/(SPC(seller));
% fprintf('卖家的出售资源比例为：%.3f\n',y_percentage)
% fprintf('紧急程度最高买家(8)占卖家的出售资源比例为：%.3f\n',yt_percentage)

disp('卖家s收入情况')
fprintf('总收入：%.3f\n',R(seller));
for j = 1 : m
    r_buyer(j) = (N(seller) * r(j,seller)* Ur(j) - 0.5) * (lambda(j,seller) * com(j) + u(j,seller) * spc(j));
    p_buyer(j) = r_buyer(j)/R(seller);
    fprintf('来自%d买家收入：%.3f\n',j,r_buyer(j));
    fprintf('%d买家占比：%.3f\n', j, p_buyer(j));
end




