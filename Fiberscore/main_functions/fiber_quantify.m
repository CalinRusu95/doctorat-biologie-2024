clc;
TYPE = 'cito';

INDEX_MODEL = ['temps/' TYPE '/index_clean.mat'];
ORIEN_MODEL = ['temps/' TYPE '/orient_clean.mat'];
MEAN_MODEL = ['temps/' TYPE '/mean_clean.mat'];
TRAIN_INPUT = ['input_' TYPE '.txt'];

load(INDEX_MODEL)
load(ORIEN_MODEL)
load(MEAN_MODEL)
datapath = textread(TRAIN_INPUT,'%s');

lengths = cell(1,length(datapath));
angles = cell(1,length(datapath));
number = cell(1,length(datapath));
maxnumlen = cell(1,length(datapath));

disp('Running quantitative analysis...');

tic;
for cell_id=1:size(index_clean,2)
    fprintf('On cell number: %d\n',cell_id)
    [fiber_length orient num max_num_fiber] = angle_length_extract(orient_clean{cell_id}, mean_clean{cell_id});
    SD_density = std(mean_clean{cell_id});
    lengths{cell_id} = fiber_length;
    angles{cell_id} = orient;
    number{cell_id} = num;
    maxnumlen{cell_id} = max_num_fiber;
    sd_density{cell_id} = SD_density;
end

save (['temps/' TYPE '/lengths.mat'], 'lengths')
save (['temps/' TYPE '/angles.mat'], 'angles')
save (['temps/' TYPE '/number.mat'], 'number')
save (['temps/' TYPE '/maxnumlen.mat'], 'maxnumlen')
save (['temps/' TYPE '/sd_density.mat'], 'sd_density')
fprintf('Took %g minutes to quantify the fibers and find the orientation.\n\n',toc/60);