function weightVectors = updateWeightVectors(weightVectors, paretoFront)
    numSubproblems = size(weightVectors, 1);
    numObjectives = size(weightVectors, 2);
    
    % Calculate the centroid of Pareto front
    centroid = mean(paretoFront(:, 1:numObjectives), 1);
    
    % Update weight vectors based on the centroid
    newWeightVectors = repmat(centroid, numSubproblems, 1);
    
    % Perturb the new weight vectors
    perturbation = randn(numSubproblems, numObjectives) * 0.1; % Adjust the perturbation factor as needed
    newWeightVectors = newWeightVectors + perturbation;
    
    % Normalize the new weight vectors
    newWeightVectors = newWeightVectors ./ sum(newWeightVectors, 2);
    
    % Replace the old weight vectors with the new ones
    weightVectors = newWeightVectors;
end
