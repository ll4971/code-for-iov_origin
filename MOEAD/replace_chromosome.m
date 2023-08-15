function f  = replace_chromosome(allchrom, M, D,NP)
N = size(allchrom); 
% 根据rank值进行排序
[~,index] = sort(allchrom(:,M + D + 1));
sortedchrom = allchrom(index,:);
max_rank = max(allchrom(:,M + D + 1));
previous_index = 0;
for i = 1 : max_rank
    current_index = max(find(sortedchrom(:,M + D + 1) == i));
    if current_index > NP
        remaining = NP - previous_index;
        temp_pop = ...
            sortedchrom(previous_index + 1 : current_index, :);
        [temp_sort,temp_sort_index] = ...
            sort(temp_pop(:, M + D + 2),'descend');
        for j = 1 : remaining
            f(previous_index + j,:) = temp_pop(temp_sort_index(j),:);
        end
        return;
    elseif current_index < NP
        f(previous_index + 1 : current_index, :) = ...
            sortedchrom(previous_index + 1 : current_index, :);
    else
        f(previous_index + 1 : current_index, :) = ...
            sortedchrom(previous_index + 1 : current_index, :);
        return;
    end
    previous_index = current_index;
end
