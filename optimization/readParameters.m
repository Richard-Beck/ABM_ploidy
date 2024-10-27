
function params = readParameters(filename)
%Code only reads birth,death and tEnd values for nwo
    % Open the file for reading
    fileID = fopen(filename, 'r');
    
    % Read the entire file as a set of strings and numerical values
    data = textscan(fileID, '%s %f');
    
    % Close the file after reading
    fclose(fileID);
    
    % Initialize the parameters struct
    params = struct();
    
    % Loop through all the read data and assign the values to params
    for i = 1:length(data{1})
        if strcmp(data{1}{i}, 'natural_growth_rate')
            params.natural_growth_rate = data{2}(i);
        elseif strcmp(data{1}{i}, 'natural_death_rate')
            params.natural_death_rate = data{2}(i);
         elseif strcmp(data{1}{i}, 'tEnd')
            params.tEnd = data{2}(i);
        elseif strcmp(data{1}{i}, 'vssl_bdary_value')
            params.vssl_bdary_value = data{2}(i);
        elseif strcmp(data{1}{i}, 'initialTumorSize')
            params. initialTumorSize = data{2}(i);
        end
    end
    
  
    % Provide default values if they were not found in the file
    if ~isfield(params, 'natural_growth_rate')
        params.growth_rate = NaN; % or some default value
    end
    if ~isfield(params, 'natural_death_rate')
        params.death_rate = NaN; % or some default value
    end
     if ~isfield(params, 'tEnd')
        params.tEnd = NaN; % or some default value
    end
end
