function chrom = initPop(NP, M, dim, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p)
    % K 是数组元素总个数将决策变量和目标函数值组合一个数组
    K = M + dim;
    
    % 初始化每条染色体
    chrom = zeros(NP, K);
    for i = 1 : NP
        % 初始化种群
        while 1
            x = randi(n, 1, m); % 生成 m 个 1~n 的随机数，索引买家对应的卖家
            cs = zeros(1, n); % 卖家计算资源实际供应量
            for k = 1 : n
                cs(k) = sum(com(x == k)); % 卖家 k 的供应量
            end
            if all(cs <= COM)
                break
            end
        end
        while 1
            y = randi(n, 1, m);
            ss = zeros(1, n); % 卖家频谱资源实际供应量
            for k = 1 : n
                ss(k) = sum(spc(y == k));
            end
            if all(ss <= SPC)
                break
            end
        end
        chrom(i, 1:dim) = [x y];
        
        % 初始化染色体上的目标函数值
        objectives = fitness(chrom(i, :), m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p);
        chrom(i, dim + 1 : K) = objectives;
    end
end
