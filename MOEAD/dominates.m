function isDominating = dominates(individualA, individualB)
    % 检查个体 A 是否支配个体 B
    objectivesA = individualA;
    objectivesB = individualB;

    % 检查在每个目标函数上的优劣关系
    betterInSomeObjectives = any(objectivesA < objectivesB); % 至少有一个目标函数值更小
    worseInNoObjectives = all(objectivesA <= objectivesB); % 没有一个目标函数值更大
    isDominating = betterInSomeObjectives && worseInNoObjectives; % 同时满足两个条件
end
