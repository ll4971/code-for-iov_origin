function [com, spc, COM, SPC, Ur, r, N, D, x0] = generate_data(m, n, swt, rep)
% 生成服从高斯分布的需求和供给数据
% 输入参数：
% m - 买家的数量
% n - 卖家的数量
% 输出参数：
% com - 买家的计算资源需求，1 x m 的向量
% spc - 买家的频谱资源需求，1 x m 的向量
% COM - 卖家的计算资源供给，1 x n 的向量
% SPC - 卖家的频谱资源供给，1 x n 的向量
% swt - 买家的紧急程度，1 x m 的向量
% r - 买家对卖家的信誉值，10 x 7 的矩阵
% N - 卖家的卖价，1 x n 的向量
% D - 买家和卖家之间的距离矩阵，10 x 7 的矩阵
% x0 - 买家的任务量，1 x m 的向量

% 设置随机数种子
rng(8);

% 设置需求和供给的均值和标准差
mu_com = 50; % 买家计算资源需求的均值
sigma_com = 10; % 买家计算资源需求的标准差
mu_spc = 30; % 买家频谱资源需求的均值
sigma_spc = 10; % 买家频谱资源需求的标准差
mu_COM = 1000; % 卖家计算资源供给的均值
sigma_COM = 20; % 卖家计算资源供给的标准差
mu_SPC = 800; % 卖家频谱资源供给的均值
sigma_SPC = 20; % 卖家频谱资源供给的标准差
mu_N = 150; % 卖家卖价的均值
sigma_N = 20; % 卖家卖价的标准差
mu_x0 = 200; % 买家任务量的均值
sigma_x0 = 50; % 买家任务量的标准差

% 生成服从高斯分布的需求和供给数据
com = round(normrnd(mu_com, sigma_com, 1, m)); % 买家计算资源需求
spc = round(normrnd(mu_spc, sigma_spc, 1, m)); % 买家频谱资源需求
COM = round(normrnd(mu_COM, sigma_COM, 1, n)); % 卖家计算资源供给
SPC = round(normrnd(mu_SPC, sigma_SPC, 1, n)); % 卖家频谱资源供给
N = round(normrnd(mu_N, sigma_N, 1, n)); % 卖家卖价
x0 = round(normrnd(mu_x0, sigma_x0, 1, m)); % 买家任务量

% 确保需求和供给都是正数
com = max(com, 0);
spc = max(spc, 0);
COM = max(COM, 0);
SPC = max(SPC, 0);
N = max(N, 0);
x0 = max(x0, 0);

%% 生成 Ur 变量
if swt == 1
    Ur = [0.1 0.6 1]*5; % 可选的紧急程度
    Ur = Ur (randi (3,1,m)); % 随机为每个买家分配一个紧急程度
elseif swt ==0
    Ur = ones(1,m);
end

if rep == 0
%% 生成 r 变量
r = randi([1, 10], m, n); % 生成一个 10 x 7 的矩阵，服从 [0, 10] 区间的均匀分布
elseif rep == 1
end
%% 生成服从对数正太1-80的D
    % 生成一个大小为 10x7 的矩阵，元素为随机数
    % 使用 lognrnd 函数生成对数正态分布随机数
    % lognrnd(mu, sigma, m, n)生成大小为 m x n 的对数正态分布随机数矩阵
mu_dis = log(20); % 设置对数正态分布的均值，即 log(40) ≈ 3.68888
sigma_dis = log(2); % 设置对数正态分布的标准差，即 log(2) ≈ 0.693147
D = round(lognrnd(mu_dis, sigma_dis, m, n),2); % 距离保留两位小数
end
