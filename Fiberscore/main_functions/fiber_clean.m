clc; clear all;

% Set the type and file paths
TYPE = 'cito';

% Define paths for models and input files
modelPaths = defineModelPaths(TYPE);
TRAIN_INPUT = ['include/input_' TYPE '.txt'];
datapath = textread(TRAIN_INPUT, '%s');

% Load the model data
[index_cell, orientation_cell, mean_cell] = loadModelData(modelPaths);

% Initialize parameters
T = 0;
min_length = 5;
K = 10;

disp('Cleaning images...');
tic;

% Process each cell image
for cell_id = 1:numel(index_cell)
    processCellImage(cell_id, datapath, index_cell, mean_cell, orientation_cell, T, min_length, K, TYPE);
end

fprintf('All images processed. Total time: %g minutes.\n', toc / 60);

% Helper Functions

function modelPaths = defineModelPaths(TYPE)
    basePath = '...\temps\';
    modelPaths.INDEX_MODEL = fullfile(basePath, TYPE, 'index_cell.mat');
    modelPaths.ORIEN_MODEL = fullfile(basePath, TYPE, 'orientation_cell.mat');
    modelPaths.MEAN_MODEL = fullfile(basePath, TYPE, 'mean_cell.mat');
end

function [index_cell, orientation_cell, mean_cell] = loadModelData(modelPaths)
    % Load the models
    load(modelPaths.INDEX_MODEL, 'index_cell');
    load(modelPaths.ORIEN_MODEL, 'orientation_cell');
    load(modelPaths.MEAN_MODEL, 'mean_cell');
end

function processCellImage(cell_id, datapath, index_cell, mean_cell, orientation_cell, T, min_length, K, TYPE)
    fprintf('Processing Cell: %d\n', cell_id);

    % Extract data for the current cell
    indexfs = index_cell{cell_id};
    meanfs = mean_cell{cell_id};
    orientfs = orientation_cell{cell_id};
    image_name = datapath{cell_id};

    % Clean the cell image
    [indexfs, meanfs, orientfs] = cleanCellImage(indexfs, meanfs, orientfs, T, min_length, K);

    % Convert and save the cleaned data
    saveCleanedData(cell_id, indexfs, meanfs, orientfs, image_name, TYPE);
end

function [indexfs, meanfs, orientfs] = cleanCellImage(indexfs, meanfs, orientfs, T, min_length, K)
    BW = indexfs > 0;
    indexfs(~BW) = 0;

    % Perform directional morphological cleaning
    BW_final = directional_morph_clean(indexfs, min_length, T, K);

    % Update index, mean, and orientation fields
    orientfs(~BW_final) = -100;
    meanfs(~BW_final) = 0;
    indexfs(~BW_final) = 0;

    % Remove small objects
    BW = bwareaopen(indexfs > 0, min_length);
    orientfs(~BW) = -100;
    meanfs(~BW) = 0;
    indexfs(~BW) = 0;
end

function saveCleanedData(cell_id, indexfs, meanfs, orientfs, image_name, TYPE)
    % Convert index to binary
    indexBW = im2bw(indexfs, 0.1);

    % Save the cleaned data
    index_clean{cell_id} = indexBW;
    mean_clean{cell_id} = meanfs;
    orient_clean{cell_id} = orientfs;

    % Save the binary image
    imwrite(indexBW, fullfile('...\temp_img\thin\', [image_name, '.jpg']));

    % Save the cleaned data to MAT files
    save(fullfile('...\temps\', TYPE, 'index_clean.mat'), 'index_clean');
    save(fullfile('...\temps\', TYPE, 'mean_clean.mat'), 'mean_clean');
    save(fullfile('...\temps\', TYPE, 'orient_clean.mat'), 'orient_clean');

    fprintf('Finished processing Cell %d. Time elapsed: %g minutes.\n\n', cell_id, toc / 60);
end
