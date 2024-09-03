% Main Script for Processing Cell Images
TYPE = 'cito'; 
disp('Reading input files...');

% Load input file paths
TRAIN_INPUT = ['input_' TYPE '.txt'];
datapath = loadInputPaths(TRAIN_INPUT);
celldisp(datapath);

% Open results file for writing
resultsFile = '...\results.txt';
res = fopen(resultsFile, 'w+');

% Process each image in the input list
for i = 1:length(datapath)
    processImage(datapath{i}, res);
end

% Close the results file
fclose(res);

% Helper Functions

function datapath = loadInputPaths(inputFile)
    % Load the list of image file paths from the input file
    datapath = textread(inputFile, '%s');
end

function processImage(imagePath, res)
    % Display the current image being processed
    disp('Reading image...');
    disp(['Processing image: ', imagePath]);

    % Read the image
    dataIn = imread(string(imagePath));

    % Perform cell migration analysis
    [Res_stats, Res_colour, Res_gray] = cellMigration(dataIn);

    % Display and log results
    displayResults(imagePath, Res_stats, res);

    % Display the processed image with color overlay
    showImageWithOverlay(imagePath, Res_colour);
end

function displayResults(imagePath, Res_stats, res)
    % Print and log statistics for the current image
    fprintf('On Cell: %s\n', imagePath);
    disp(Res_stats);
    fprintf(res, '%f\n', Res_stats.avDist);
end

function showImageWithOverlay(imagePath, Res_colour)
    % Display the image with color overlay in a new figure
    figure('Name', string(imagePath));
    imagesc(Res_colour);
    axis off;
end
