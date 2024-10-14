function [stop, optimValues, changed]= saveOptimizationParams(optimValues, vargin,state)



 disp('optim')
vargin

allParams=struct();
disp('Output fun')
    % Define the file to write parameters
    fileName = '/Users/4477116/Documents/projects/Ploidy_abm/input/parameters.txt';

    % Initialize stop flag
     stop = false;
     changed=false;


    % Check the state of the optimization
    switch state
        case 'iter'
            % Write the parameters to the file at each iteration
             %currentParams =optimValues.x;  % When using patternsearch
             %currentParams =optimValues;   % When using fmincon
              %currentParams=updatedOptimValues; %
               currentParams=vargin.x; %When running imulated Annealing solver

                paramCellArray = struct('natural_growth_rate', currentParams(1),...
                                    'natural_death_rate', currentParams(2),...
                                    'initialTumorSize',currentParams(3),...
                                    'vssl_bdary_value',currentParams(4));
              
            
            % paramCellArray = struct('natural_growth_rate', currentParams(1), 'natural_death_rate', currentParams(2));
            writeParameters(fileName, paramCellArray);
    end


end
