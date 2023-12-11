function f = nonDominatedSort(x, M, dim)
% 该函数是对当前种群进行非支配排序。第一层(Front)的个体记为rank 1，第二层个
% 体对应rank 2，依次类推。指定rank后，计算每一层的拥挤度. 
N = size(x,1);
front = 1;    % 当前计算层为第一层.
F(front).f = [];
individual = [];
%% 非支配排序
for i = 1 : N
    individual(i).np = 0;   % 支配该个体 i 的个体数量
    individual(i).Sp = [];  % 被该个体支配的个体集合
    for j = 1 : N
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M      % 对每一个目标函数值进行比较
            if (x(i,dim + k) < x(j,dim + k))
                dom_less = dom_less + 1;
            elseif (x(i,dim + k) == x(j,dim + k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M    % 后面这个条件应该是去除本身的意思，就是i~=j，应该是这样的
            individual(i).np = individual(i).np + 1;   % q 支配 p
        elseif dom_more == 0 && dom_equal ~= M
            individual(i).Sp = [individual(i).Sp j];   % p 支配 q
        end
    end   
    if individual(i).np == 0
        x(i,M + dim + 1) = 1;          % 将层号写入染色体的最后一位
        F(front).f = [F(front).f i]; % 第一层个体集合
    end
end

% Find the subsequent fronts
while ~isempty(F(front).f)
   Q = [];
   for i = 1 : length(F(front).f)    % 对front层中所有个体进行如下循环
       if ~isempty(individual(F(front).f(i)).Sp)   % front层中第i个个体的Sp不为空
        	for j = 1 : length(individual(F(front).f(i)).Sp) % 对front层中i个个体的Sp中所有个体进行如下循环
            	individual(individual(F(front).f(i)).Sp(j)).np = individual(individual(F(front).f(i)).Sp(j)).np - 1;
        	   	if individual(individual(F(front).f(i)).Sp(j)).np == 0
               		x(individual(F(front).f(i)).Sp(j),M + dim + 1) = front + 1;
                    Q = [Q individual(F(front).f(i)).Sp(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end

[temp,index_of_fronts] = sort(x(:,M + dim + 1));
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
end
current_index = 0;

%% 计算拥挤度crowding distance
for front = 1 : (length(F) - 1)
%    objective = [];
    distance = 0;
    y = [];
    previous_index = current_index + 1;
    for i = 1 : length(F(front).f)
        y(i,:) = sorted_based_on_front(current_index + i,:);
    end
    current_index = current_index + i;
    sorted_based_on_objective = [];
    for i = 1 : M
        [sorted_based_on_objective, index_of_objectives] = sort(y(:,dim + i));
        sorted_based_on_objective = [];
        for j = 1 : length(index_of_objectives)
            sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
        end
        f_max =  sorted_based_on_objective(length(index_of_objectives), dim + i);
        f_min = sorted_based_on_objective(1, dim + i);
        y(index_of_objectives(length(index_of_objectives)),M + dim + 1 + i) = Inf;
        y(index_of_objectives(1),M + dim + 1 + i) = Inf;
         for j = 2 : length(index_of_objectives) - 1
            next_obj  = sorted_based_on_objective(j + 1,dim + i);
            previous_obj  = sorted_based_on_objective(j - 1,dim + i);
            if (f_max - f_min == 0)
                y(index_of_objectives(j),M + dim + 1 + i) = Inf;
            else
                y(index_of_objectives(j),M + dim + 1 + i) =  (next_obj - previous_obj)/(f_max - f_min);
            end
         end
    end
    distance = [];
    distance(:,1) = zeros(length(F(front).f),1);
    for i = 1 : M
        distance(:,1) = distance(:,1) + y(:,M + dim + 1 + i);
    end
    y(:,M + dim + 2) = distance;
    y = y(:,1 : M + dim + 2);
    z(previous_index:current_index,:) = y;
end
f = z();
