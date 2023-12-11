function newchrom  = geneOperator(parentchrom, M, dim, Pm, Pc, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v, ka, epsilon, sigma, p)
% 函数根据父代个体产生子代个体
NP = size(parentchrom,1);
newchrom = parentchrom(:,1:M+dim);
% 交叉操作
for i = 1 :2: (NP-rem(NP,2))
    % 父代染色体
    ch1 = parentchrom(i,1:dim);
    ch2 = parentchrom(i+1,1:dim);
    if rand < Pc
        r0 = randi(dim);
        newx1 = [ch1(1:r0) ch2(r0+1:end)];
        newx2 = [ch2(1:r0) ch1(r0+1:end)];
        ch1 = newx1;
        ch2 = newx2;
       
    end
    newchrom(i,1:dim) = ch1;
    newchrom(i+1,1:dim) = ch2;
    % 计算子代个体的目标函数值
    newchrom(i,dim+1:M+dim) = fitness(newchrom(i,1:dim), m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p);
    newchrom(i+1,dim+1:M+dim) = fitness(newchrom(i+1,1:dim), m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p);
end
for i = 1 : NP
    x = parentchrom(i,1:dim);
    if rand < Pm
        r0 = randi(dim);
        x(r0) = randi(n);
    end
    ch(1:dim) = x;
    ch(1,dim+1:M+dim) = fitness(ch, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p);
    newchrom(i,:) = ch;
end
end
