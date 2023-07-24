 % 确定种群中的非支配解 
function Idx = getNondominated(fx)
NP = size(fx,1);
Idx = ones(NP,1);
for i = 1 : NP
    for j = 1 : NP
        if i ~= j
            if dominates(fx(j,:), fx(i,:))
                Idx(i) = 0;
            end
        end
    end
end
Idx = (Idx == 1);

