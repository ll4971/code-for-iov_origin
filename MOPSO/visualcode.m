K = 1; 
ch = gbest(K,1:dim);
x = ch(1:m);        % 计算量分配
y = ch(m+1:m*2);    % 传输量分配
lambda = zeros(m, n);
for i = 1 : m
   lambda(i,x(i)) = 1; 
end
u = zeros(m, n);
for i = 1 : m
   u(i,y(i)) = 1; 
end
C = zeros(1, m);
for i = 1 : m
   C(i) = 0;
   for j = 1 : n
      C(i) = C(i) +  N(j) / r(i,j) * (lambda(i,j) * com(i) + u(i,j) * spc(i));
   end
end
% 
R = zeros(1, n);
for j = 1 : n
   R(j) = 0;
   for i = 1 : m
      R(j) = R(j) +  (N(j) * r(i,j) - 0.5) * (lambda(i,j) * com(i) + u(i,j) * spc(i));
   end
end

U = -rho * sum(C) + v * sum(R);
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
disp('计算量分配：')
for j = 1 : m
    fprintf('买家 %d —— 卖家 %d \n',j, x(j))
end
disp('传输量分配：')
for j = 1 : m
    fprintf('买家 %d —— 卖家 %d \n',j, x(j))
end
fprintf('C = %.3f\n', sum(C))
fprintf('R = %.3f\n', sum(R))
fprintf('U = %.3f\n', U)
fprintf('E = %.3f\n', E)
fit = [-U E] + gfun * 1e8;