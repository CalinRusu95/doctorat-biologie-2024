clc; clear;

% Set parameters and configurations
TYPE = 'cito';
tic;

% Configuration Parameters
params = defineParameters();

disp('Reading input files...');
TRAIN_INPUT = ['...\include\input_' TYPE '.txt'];
datapath = textread(TRAIN_INPUT, '%s');

% Initialize cell arrays to store results
[index_cell, orientation_cell, mean_cell] = initializeCellArrays(datapath);

% Process each image
for imgid = 1:length(datapath)
    processImage(datapath{imgid}, imgid, params, index_cell, orientation_cell, mean_cell, TYPE);
end

fprintf('All images processed. Total time: %g minutes.\n', toc / 60);

% Helper Functions

function params = defineParameters()
    % Define the parameters used for processing
    params.L = 6;      % Kernel size
    params.K = 10;     % Angular resolution
    params.TC = 0.55;  % Threshold for correlation coefficient with a Gaussian profile
    params.M = 0.3;    % Threshold for NSD
    params.N = 2.1;    % Threshold for ratio of NSD between perpendicular rods
    params.T = 0.25;   % Threshold for background-subtracted intensity value for a fiber
    params.w = 1;      % One pixel value
end

function [index_cell, orientation_cell, mean_cell] = initializeCellArrays(datapath)
    % Initialize cell arrays to store results
    index_cell = cell(size(datapath));
    orientation_cell = cell(size(datapath));
    mean_cell = cell(size(datapath));
end

function processImage(imageName, imgid, params, index_cell, orientation_cell, mean_cell, TYPE)
    % Display the current image being processed
    fprintf('Processing Image %d: %s\n', imgid, imageName);
    
    % Read and preprocess the image
    original = readAndPreprocessImage(imageName);
    
    % Calculate fiber index, orientation, and mean fiber strength
    [index, orientation, mean_fs] = fiberscore_conv(imageName, original, params.K, params.L, params.TC, params.M, params.N, params.T);
    
    % Store results in cell arrays
    index_cell{imgid} = index;
    orientation_cell{imgid} = orientation;
    mean_cell{imgid} = mean_fs;
    
    % Save the processed image and parameters
    saveProcessedImageAndData(imageName, imgid, index, datapath, index_cell, orientation_cell, mean_cell, params, TYPE);
    
    fprintf('Took %g minutes to extract the fibers and find the orientation.\n\n', toc / 60);
end

function original = readAndPreprocessImage(imagePath)
    % Read the image
    IMAGE = imread(['...\cito\', imagePath]);
    [x, y, D] = size(IMAGE);

    % Convert to grayscale if needed
    if D == 1
        original = IMAGE;
    elseif D >= 3
        original = rgb2gray(IMAGE(:, :, 1:3));
    end

    % Normalize the image
    original = double(original) / max(max(double(original)));
end

function saveProcessedImageAndData(imageName, imgid, index, datapath, index_cell, orientation_cell, mean_cell, params, TYPE)
    % Save the processed image
    imwrite(index, fullfile('...\temp_img\conv\', [datapath{imgid}, '.tif']), 'tif');
    
    % Save the results and parameters
    save(fullfile('...\temps\', TYPE, 'index_cell.mat'), 'index_cell');
    save(fullfile('...\temps\', TYPE, 'orientation_cell.mat'), 'orientation_cell');
    save(fullfile('...\temps\', TYPE, 'mean_cell.mat'), 'mean_cell');
    save(fullfile('...\temps\', TYPE, 'param.mat'), 'params');
end
