function [stop ,optimValues , changed ]  = saveOptimizationParamsSA(optold, optimValues, state)
    stop = false;  % Continue the optimization by default
    changed=false;
    disp('called')
    optimValues
    optold
    state

    global paramHistory;

    % Check the state and log data if necessary
    if strcmp(state, 'iter') || strcmp(state, 'done')
        % Open the file in append mode

        paramHistory

      paramCellArray =optimValues.x
         paramCellArray = struct('natural_growth_rate', currentParams(1),...
                                    'natural_death_rate', currentParams(2),...
                                    'initialTumorSize',currentParams(3),...
                                    'vssl_bdary_value',currentParams(4));



        % Write all recorded parameters and objective values to the file
        for i = 1:size(paramHistory, 1)
            fprintf(fileID, 'Iteration: %d, Objective: %f, Parameters: %s\n', ...
                    i, paramHistory(i, end), mat2str(paramHistory(i, 1:end-1)));
        end

        % Clear the history after saving
        paramHistory = [];

        % Close the file
        fclose(fileID);
    end
end
