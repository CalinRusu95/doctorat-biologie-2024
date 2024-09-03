clc; clear all;
TYPE = 'cito';

INDEX_MODEL = ['...\temps\' TYPE '\index_cell.mat'];
ORIEN_MODEL = ['...\temps\' TYPE '\orientation_cell.mat'];
MEAN_MODEL = ['...\temps\' TYPE '\mean_cell.mat'];
MIN_FIBER = ['...\temps\' TYPE '\fiber_length.mat'];
INTENSITY = ['...\temps\' TYPE '\intensity.mat'];
TRAIN_INPUT = ['include/input_' TYPE '.txt'];
datapath = textread(TRAIN_INPUT,'%s');

load(INDEX_MODEL)
load(ORIEN_MODEL)
load(MEAN_MODEL)

T = 0;
min_length = 5;
K = 10;

disp('Cleaning images...');

tic;
for cell_id = 1:size(index_cell)
    fprintf('On Cell: %d\n',cell_id);
        
    indexfs = index_cell{cell_id};
    meanfs = mean_cell{cell_id};
    orientfs = orientation_cell{cell_id};
    image_name = datapath{cell_id};


    [MM NN] = size(indexfs);
    BW = indexfs >0;
    indexfs(BW==0) = 0;
        
    [BW_final] = directional_morph_clean(indexfs, min_length, T, K);

    orientfs(BW_final==0) = -100; 
    meanfs(BW_final==0) = 0;  
    indexfs(BW_final==0) = 0; 

   
    BW = bwareaopen(indexfs > 0, min_length);
    orientfs(BW==0) = -100;
    meanfs(BW==0) = 0;  
    indexfs(BW==0) = 0;

    indexBW = im2bw(indexfs,0.1);
    index_clean{cell_id} = indexBW;
    mean_clean{cell_id} = meanfs;
    orient_clean{cell_id} = orientfs;

    imwrite(indexBW, strcat('...\temp_img\thin\',image_name),'jpg');

    save (['...\temps\' TYPE '\index_clean.mat'], 'index_clean')
    save (['...\temps\' TYPE '\mean_clean.mat'], 'mean_clean')
    save (['...\temps\' TYPE '\orient_clean.mat'], 'orient_clean')
    fprintf('Took %g minutes to obtain the clean images.\n\n',toc/60);
end