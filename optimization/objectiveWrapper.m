function [fval] = objectiveWrapper(params, target_population, alpha, beta, current_params)
    % Global variable to store parameters and objective values
    global paramHistory;
   

    % Call the actual objective function
    fval = objectiveFun(params, target_population, alpha, beta, current_params);

    % Record the parameters and objective value
    paramHistory = [paramHistory; params, fval];
    writeErrors('/Users/4477116/Documents/projects/Ploidy_abm/output/errors.txt',paramHistory)

end 

