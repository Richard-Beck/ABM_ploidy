function writeParameters(filename, params)


    % Read the existing file
    fileID = fopen(filename, 'r');
    
    % Check if file was opened successfully
    if fileID == -1
        error('Error: Could not open the file. Check the file path and name.');
    end
    
    % Read the entire file as a set of strings and numerical values (as strings)
    data = textscan(fileID, '%s %s');
    fclose(fileID);
    
    % Convert numerical strings to their respective data types
    values = cell(size(data{2}));
    originalData = data{2};  % Keep original data for formatting
    
    for i = 1:length(data{2})
        values{i} = str2double(data{2}{i});  % Convert to double
    end
    
    % % Update the values for growth_rate and death_rate
    % for i = 1:length(data{1})
    %     if strcmp(data{1}{i}, 'natural_growth_rate')
    %         values{i} = params.natural_growth_rate; % Update with the new parameter value
    %         originalData{i} = sprintf('%.6f', params.natural_growth_rate); % Preserve formatting
    %     elseif strcmp(data{1}{i}, 'natural_death_rate')
    %         values{i} = params.natural_death_rate; % Update with the new parameter value
    %         originalData{i} = sprintf('%.6f', params.natural_death_rate); % Preserve formatting
    %     elseif strcmp(data{1}{i}, 'initialTumorSize')
    %         values{i} = params.initialTumorSize; % Update with the new parameter value
    %         originalData{i} = sprintf('%d', params.initialTumorSize); % Preserve formatting
    %     elseif strcmp(data{1}{i}, 'vssl_bdary_value')
    %         values{i} = params.vssl_bdary_value; % Update with the new parameter value
    %         originalData{i} = sprintf('%d', params.vssl_bdary_value); % Preserve formatting
    %     end
    % end

    
for i = 1:length(data{1})
    if strcmp(data{1}{i}, 'natural_growth_rate')
        values{i} = params.natural_growth_rate; % Update with the new parameter value
        originalData{i} = sprintf('%.6f', params.natural_growth_rate); % Preserve formatting
    elseif strcmp(data{1}{i}, 'natural_death_rate')
        values{i} = params.natural_death_rate; % Update with the new parameter value
        originalData{i} = sprintf('%.6f', params.natural_death_rate); % Preserve formatting
    elseif strcmp(data{1}{i}, 'initialTumorSize')
        % Round and format as an integer
        roundedValue = round(params.initialTumorSize);
        values{i} =round(params.initialTumorSize); % Update with the new parameter value
        originalData{i} = sprintf('%d', roundedValue); % Preserve formatting as integer
    elseif strcmp(data{1}{i}, 'vssl_bdary_value')
        % Round and format as an integer
        roundedValue = round(params.vssl_bdary_value);
        values{i} = roundedValue; % Update with the new parameter value
        originalData{i} = sprintf('%d', roundedValue); % Preserve formatting as integer
    end
end




    
    % Write the updated contents back to the file
    fileID = fopen(filename, 'w');
    
    % Check if file was opened successfully
    if fileID == -1
        error('Error: Could not open the file for writing.');
    end
    
    % Write each parameter and its value to the file, preserving original formatting
    for i = 1:length(data{1})
        % Use originalData to maintain formatting
        fprintf(fileID, '%s %s\n', data{1}{i}, originalData{i});
    end
    
    fclose(fileID);
end
