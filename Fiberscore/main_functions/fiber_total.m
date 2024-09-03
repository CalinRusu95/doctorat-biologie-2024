clc; clear all;

% Set the type and file paths
TYPE = 'cito';

% Define paths for models, input files, and output files
modelPaths = defineModelPaths(TYPE);
TRAIN_INPUT = ['include/input_' TYPE '.txt'];
datapath = textread(TRAIN_INPUT, '%s');

% Load the required model data
[lengths, angles, index_clean, orient_clean, mean_clean] = loadModelData(modelPaths);

% Prepare for error handling and open result files
dbstop if error
[fileID, fileIDQV] = openResultFiles();

disp('Processing images...');

% Write the header for the CSV file
fprintf(fileIDQV, 'image_name|number_fibers|total_length|mean_length|polarity \n');

% Initialize result arrays
total_fiber_length = zeros(1, length(datapath));
mean_fiber_length = zeros(1, length(datapath));
polarity = zeros(1, length(datapath));

tic;

% Process each image
for i = 1:length(datapath)
    processImage(i, datapath{i}, mean_clean, lengths, angles, orient_clean, fileIDQV, total_fiber_length, mean_fiber_length, polarity);
end

% Close the result files
fclose(fileIDQV);

fprintf('Took %g minutes to obtain the totals for the images.\n\n', toc / 60);

% Helper Functions

function modelPaths = defineModelPaths(TYPE)
    % Define the paths for the required models
    basePath = 'temps/';
    modelPaths.INDEX_MODEL = fullfile(basePath, TYPE, 'index_clean.mat');
    modelPaths.ORIEN_MODEL = fullfile(basePath, TYPE, 'orient_clean.mat');
    modelPaths.MEAN_MODEL = fullfile(basePath, TYPE, 'mean_clean.mat');
    modelPaths.FIBER = fullfile(basePath, TYPE, 'lengths.mat');
    modelPaths.ORIENT = fullfile(basePath, TYPE, 'angles.mat');
end

function [lengths, angles, index_clean, orient_clean, mean_clean] = loadModelData(modelPaths)
    % Load the necessary model data
    load(modelPaths.FIBER, 'lengths');
    load(modelPaths.ORIENT, 'angles');
    load(modelPaths.INDEX_MODEL, 'index_clean');
    load(modelPaths.ORIEN_MODEL, 'orient_clean');
    load(modelPaths.MEAN_MODEL, 'mean_clean');
end

function [fileID, fileIDQV] = openResultFiles()
    % Open result files for writing
    fileID = fopen('Total_results_bend.txt', 'w');
    fprintf('Results are displayed in Total_results_bend.txt\n');
    fileIDQV = fopen('fiberscore_data.csv', 'w');
end

function processImage(i, image_name, mean_clean, lengths, angles, orient_clean, fileIDQV, total_fiber_length, mean_fiber_length, polarity)
    fprintf('Processing cell number: %d\n', i);

    % Read and preprocess the image
    original = imread(image_name);    
    orig_eq = adapthisteq(original);
    
    % Process the binary mask
    BW = mean_clean{i} ~= 0;   

    % Write image name to the CSV file
    fprintf(fileIDQV, '%10s|', image_name);

    % Calculate fiber statistics
    fibers_number = length(lengths{i});
    total_fiber_length(i) = sum(lengths{i});
    mean_fiber_length(i) = mean(lengths{i}); 

    % Calculate fiber orientation and polarity
    total_mean = sum(mean_clean{i}(BW));
    Dx = sum(mean_clean{i}(BW) .* cos(2 * orient_clean{i}(BW))) / total_mean;
    Dy = sum(mean_clean{i}(BW) .* sin(2 * orient_clean{i}(BW))) / total_mean;
    polarity(i) = sqrt(Dx^2 + Dy^2); 
    
    % Write results to the CSV file
    fprintf(fileIDQV, '%6.2f|', fibers_number);
    fprintf(fileIDQV, '%6.2f|', total_fiber_length(i));
    fprintf(fileIDQV, '%6.2f|', mean_fiber_length(i));
    fprintf(fileIDQV, '%.15f|\n', polarity(i));
end
