function neighbors = selectNeighbors(currentIndex, weightVectors, population, subproblemSize)
    numSubproblems = size(weightVectors, 1);
    numObjectives = size(weightVectors, 2);
    
    % Calculate distances to other weight vectors
    distances = zeros(numSubproblems, 1);
    for i = 1:numSubproblems
        distances(i) = norm(weightVectors(currentIndex, :) - weightVectors(i, :));
    end
    
    % Sort subproblems by distance in ascending order
    [~, sortedIndices] = sort(distances);
    
    % Select the closest neighbors based on subproblemSize
    neighbors = population(sortedIndices(1:subproblemSize), :);
end
