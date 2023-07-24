function fit = fitness(ch, m, n, com, spc, COM, SPC, N, r ,D, x0 , rho, v, ka, epsilon, sigma, p)
x = ch(1:m);
y = ch(m+1:m*2);
lambda = zeros(m, n); %计算资源购买比例
for i = 1 : m
   lambda(i,x(i)) = 1; 
end
u = zeros(m, n); %频谱资源购买比例
for i = 1 : m
   u(i,y(i)) = 1; 
end
C = zeros(1, m);
for i = 1 : m
   C(i) = 0;
   for j = 1 : n
      C(i) = C(i) +  N(j) / r(i,j) /ceil(i*3/7)* (lambda(i,j) * com(i) + u(i,j) * spc(i));
   end
end

R = zeros(1, n);
for j = 1 : n
   R(j) = 0;
   rnew(j)=0;
   for i = 1 : m
      R(j) = R(j) +  (N(j) * r(i,j)/ceil(i*3/7) - 0.5) * (lambda(i,j) * com(i) + u(i,j) * spc(i));
      rnew(j)=r(i,j)*10/ceil(i*3/7)+ rnew(j);
   end
end
u = zeros(1, n);
cs = zeros(1, n);
ss = zeros(1, n);
% 卖家实际供应量
for k = 1 : n
    cs(k) = sum(com(x == k));
    ss(k) = sum(spc(y == k));
    u(k)=rnew(k)*(max(cs(k)-COM(k),0)+ max(ss(k)-SPC(k),0))+u(k);   
end
U = ka *(-rho * sum(C) + v * sum(R))+ 1/ka *sum(u) ;
%
E = epsilon * sum(com .* x0);
for i = 1 : m
    h = 35.2 + 35 * log(D(i,y(i))) / log(10);
    E = E + p * x0(i) / (spc(i) * log(1 + p * h /sigma^2) / log(2));
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


fit = [-U E] + gfun * 1e8;