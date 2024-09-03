clc;

% Set the type and file paths
TYPE = 'cito';

% Define paths for models and input files
modelPaths = defineModelPaths(TYPE);
TRAIN_INPUT = ['input_' TYPE '.txt'];

% Load the model data
[index_clean, orient_clean, mean_clean] = loadCleanData(modelPaths);

% Load the input data path
datapath = textread(TRAIN_INPUT, '%s');

% Initialize cell arrays to store results
[lengths, angles, number, maxnumlen, sd_density] = initializeResultsArrays(datapath);

disp('Running quantitative analysis...');

tic;

% Process each cell image
for cell_id = 1:numel(index_clean)
    processCellData(cell_id, orient_clean, mean_clean, lengths, angles, number, maxnumlen, sd_density);
end

% Save the results
saveResults(modelPaths, lengths, angles, number, maxnumlen, sd_density);

fprintf('Took %g minutes to quantify the fibers and find the orientation.\n\n', toc / 60);

% Helper Functions

function modelPaths = defineModelPaths(TYPE)
    % Define the paths for the cleaned data models
    basePath = 'temps/';
    modelPaths.INDEX_MODEL = fullfile(basePath, TYPE, 'index_clean.mat');
    modelPaths.ORIEN_MODEL = fullfile(basePath, TYPE, 'orient_clean.mat');
    modelPaths.MEAN_MODEL = fullfile(basePath, TYPE, 'mean_clean.mat');
end

function [index_clean, orient_clean, mean_clean] = loadCleanData(modelPaths)
    % Load the cleaned data models
    load(modelPaths.INDEX_MODEL, 'index_clean');
    load(modelPaths.ORIEN_MODEL, 'orient_clean');
    load(modelPaths.MEAN_MODEL, 'mean_clean');
end

function [lengths, angles, number, maxnumlen, sd_density] = initializeResultsArrays(datapath)
    % Initialize cell arrays to store the results
    numCells = length(datapath);
    lengths = cell(1, numCells);
    angles = cell(1, numCells);
    number = cell(1, numCells);
    maxnumlen = cell(1, numCells);
    sd_density = cell(1, numCells);
end

function processCellData(cell_id, orient_clean, mean_clean, lengths, angles, number, maxnumlen, sd_density)
    % Display the current cell being processed
    fprintf('Processing cell number: %d\n', cell_id);
    
    % Extract fiber characteristics
    [fiber_length, orient, num, max_num_fiber] = angle_length_extract(orient_clean{cell_id}, mean_clean{cell_id});
    
    % Compute the standard deviation of the density
    SD_density = std(mean_clean{cell_id});
    
    % Store the results in the respective arrays
    lengths{cell_id} = fiber_length;
    angles{cell_id} = orient;
    number{cell_id} = num;
    maxnumlen{cell_id} = max_num_fiber;
    sd_density{cell_id} = SD_density;
end

function saveResults(modelPaths, lengths, angles, number, maxnumlen, sd_density)
    % Save the results to MAT files
    save(fullfile(modelPaths.INDEX_MODEL, '../lengths.mat'), 'lengths');
    save(fullfile(modelPaths.INDEX_MODEL, '../angles.mat'), 'angles');
    save(fullfile(modelPaths.INDEX_MODEL, '../number.mat'), 'number');
    save(fullfile(modelPaths.INDEX_MODEL, '../maxnumlen.mat'), 'maxnumlen');
    save(fullfile(modelPaths.INDEX_MODEL, '../sd_density.mat'), 'sd_density');
end
