function mesh = getComsolmesh(filename)

% Open the file for reading
fileID = fopen(filename, 'r');

% Read the preamble until the key is found
key1 = '# Mesh vertex coordinates';
key2 = '# number of mesh vertices';
key3 = '# Type #1';
key4 = '# number of elements';
key5 = '# Elements';
key6 = '# Type #2';
key7 = '# number of elements';
key8 = '# Elements';

%% Get node coordinates

while true
    line = fgetl(fileID);
    index = contains(line, key2);
    if index ==1
        % Extract the numeric value using regular expressions
        pattern = '(\d)';  % Regular expression pattern for decimal numbers
        matches = regexp(line, pattern, 'match');
        numericValue = str2double(strcat(matches{1:end}));
    end
    if strcmp(line, key1)
        break; % Exit the loop when the key is found
    end
    % Process the preamble data (e.g., extract numbers, strings, etc.)
    % Do something with the preamble data
end

% Read the number of rows in the array from the preamble
% numRows = fscanf(fileID, '%d', 1);
numRows_p = numericValue;

% Read the number array
dataArray = fscanf(fileID, '%f %f\n', [3, numRows_p]);
p = dataArray; % Transpose the array if needed

mesh.p = p;

%% Get boundary nodes indices

while true
    line = fgetl(fileID);
    if strcmp(line, key3)
        break; % Exit the loop when the key is found
    end
end

while true
    line = fgetl(fileID);
    index = contains(line, key4);
    if index ==1
        % Extract the numeric value using regular expressions
        pattern = '(\d)';  % Regular expression pattern for decimal numbers
        matches = regexp(line, pattern, 'match');
        numericValue = str2double(strcat(matches{1:end}));
    end
    if strcmp(line, key5)
        break; % Exit the loop when the key is found
    end
    % Process the preamble data (e.g., extract numbers, strings, etc.)
    % Do something with the preamble data
end

% Read the number array
dataArray = fscanf(fileID, '%f', [2, numericValue]);
boundary_indices = dataArray + 1; % Transpose the array if needed

boundary_indices = unique(boundary_indices);
mesh.boundary_indices = boundary_indices;


 %% Get element matrix


while true
    line = fgetl(fileID);
    if strcmp(line, key6)
        break; % Exit the loop when the key is found
    end
end

while true
    line = fgetl(fileID);
    index = contains(line, key7);
    if index ==1
        % Extract the numeric value using regular expressions
        pattern = '(\d)';  % Regular expression pattern for decimal numbers
        matches = regexp(line, pattern, 'match');
        numericValue = str2double(strcat(matches{1:end}));
    end
    if strcmp(line, key8)
        break; % Exit the loop when the key is found
    end
    % Process the preamble data (e.g., extract numbers, strings, etc.)
    % Do something with the preamble data
end

numRows_t = numericValue;

% Read the number array
dataArray = fscanf(fileID, '%f', [3, numRows_t]);
t = dataArray + 1; % Transpose the array if needed
t(4,:) = 1;

mesh.t = t;

% Close the file
fclose(fileID);
 

end