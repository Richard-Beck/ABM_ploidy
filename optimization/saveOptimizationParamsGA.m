function [updatedState, options, optchanged] = saveOptimizationParamsGA(options, state, flag)
    % This output function saves the current parameters to a file at each iteration
    

    % Initialize optchanged to false
    optchanged = false;

    % Define the file to write parameters
    fileName = '/Users/4477116/Documents/projects/Ploidy_abm/input/parameters.txt';

    % Check the state of the optimization
    switch flag
        case 'iter'
            % Find the index of the best solution (individual with the lowest fitness value)
            [~, bestIdx] = min(state.Score);
            
            % Get the parameters of the current best solution
            currentParams = state.Population(bestIdx,:);

            % Convert parameters to a structure for saving
            paramCellArray = struct('natural_growth_rate', currentParams(1),...
                                    'natural_death_rate', currentParams(2),...
                                    'initialTumorSize',currentParams(3),...
                                    'vssl_bdary_value',currentParams(4));
            
            % Write parameters to Parameters.tx file to be used by ABM
            writeParameters(fileName, paramCellArray); 
           
            
            %
            disp('Parameters saved at iteration.');
            
        case 'init'
            % Some message 
            disp('GA optimization has started.');
            
        case 'done'
            % Some message
            disp('GA optimization Done.');
    end

end
