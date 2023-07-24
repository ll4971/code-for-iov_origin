function ch = initpop(m, n, com, spc, COM, SPC)

% 初始化种群
while 1
    x = 1 + (n - 1) * rand(1, m);
    x0 = round(x);
    cs = zeros(1, n);  % 卖家实际供应量
    for k = 1 : n
        cs(k) = sum(com(x0 == k));
    end
    if all(cs <= COM)
        break
    end
end
while 1
    y = 1 + (n - 1) * rand(1, m);
    y0 = round(y);
    ss = zeros(1, n);  % 卖家实际供应量
    for k = 1 : n
        ss(k) = sum(spc(y == k));
    end
    if all(ss <= SPC)
        break
    end
end
ch = [x y];

