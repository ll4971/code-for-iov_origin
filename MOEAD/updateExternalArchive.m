function archive = updateExternalArchive(archive, offspring)
    numArchive = size(archive, 1);
    numOffspring = size(offspring, 1);
    dim = size(archive, 2) - 2;

    combinedPopulation = [archive; offspring]; % 将存档和子代个体合并
    
    % 计算支配关系和非支配排序
    [dominated, dominatedMatrix, rank] = calculateDominance(combinedPopulation, dim);

    % 选择非支配解形成新的存档
    newArchive = combinedPopulation(dominated, :);
    
    % 如果新存档超过预设的存档容量，则根据非支配排序选择一部分
    if size(newArchive, 1) > numArchive
        selectedFronts = selectFronts(newArchive, rank, numArchive);
        newArchive = combinedPopulation(selectedFronts, :);
    end
    
    archive = newArchive;
end


function [dominated, dominatedMatrix, rank] = calculateDominance(combinedPopulation, dim)
    numCombined = size(combinedPopulation, 1);
    dominated = false(numCombined, 1);
    dominatedMatrix = false(numCombined, numCombined);
    rank = zeros(numCombined, 1);

    for i = 1:numCombined
        for j = i+1:numCombined % 只需要比较后面的个体，避免重复计算
            if dominates(combinedPopulation(i, end - 1:end), combinedPopulation(j, end - 1:end)) % 只比较目标函数值
                dominatedMatrix(i, j) = true;
                dominatedMatrix(i, dominatedMatrix(j, :)) = true; % 如果 i 支配 j ，那么 i 也支配 j 支配的所有个体
            elseif dominates(combinedPopulation(j, end - 1:end), combinedPopulation(i, end - 1:end))
                dominatedMatrix(j, i) = true;
                dominatedMatrix(j, dominatedMatrix(i, :)) = true; % 如果 j 支配 i ，那么 j 也支配 i 支配的所有个体
            else
                dominatedMatrix(i, dominatedMatrix(j, :)) = false; % 如果 i 和 j 互不支配，那么 i 和 j 支配的所有个体也互不支配
                dominatedMatrix(j, dominatedMatrix(i, :)) = false;
            end
        end
        if sum(dominatedMatrix(i, :)) == 0 && sum(dominatedMatrix(:, i)) == 0 % 需要同时判断行和列是否全为0
            rank(i) = 1;
            dominated(i) = true;
        end
    end

    % 根据每个个体被支配的次数来确定它们的非支配层次
    currentRank = 2;
    while any(rank == 0)
        currentFront = find(rank == currentRank - 1); % 获取上一层的个体
        for i = find(rank == 0) % 获取还未分层的个体
            if sum(dominatedMatrix(currentFront, i)) == 0 % 如果 i 没有被上一层中的任何个体支配，那么它属于当前层
                rank(i) = currentRank;
            end
        end
        currentRank = currentRank + 1;
    end
end



function selectedFronts = selectFronts(newArchive, rank, numArchive)
    currentRank = 1;
    selectedFronts = [];

    while length(selectedFronts) < numArchive
        currentFront = find(rank == currentRank);
        remainingSpace = numArchive - length(selectedFronts);

        if length(currentFront) <= remainingSpace
            selectedFronts = [selectedFronts currentFront];
        else
            [~, frontOrder] = sortrows(newArchive(currentFront, end), 'descend'); % 根据最后一个目标值降序排序
            selectedFronts = [selectedFronts currentFront(frontOrder(1:remainingSpace))];
            break;
        end

        currentRank = currentRank + 1;
    end
end
