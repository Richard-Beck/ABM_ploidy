function writeParamsHistory(fileName, params)


% Check if the file exists
if isfile(fileName)
    % Read the existing data from the file
    existingData = dlmread(fileName);
    
    % Append the new vector to the existing data
    updatedData = [existingData; params];
else
    % If the file does not exist, start with the new vector
    updatedData = params;
end

% Write the updated data back to the file
dlmwrite(fileName, updatedData, 'delimiter', '\t', 'precision', '%.6f');

% Display the updated content
disp('Updated content of the file:');
disp(updatedData);

end