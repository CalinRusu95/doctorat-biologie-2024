clc; clear all;
TYPE = 'cito';

INDEX_MODEL = ['temps/' TYPE '/index_clean.mat'];
ORIEN_MODEL = ['temps/' TYPE '/orient_clean.mat'];
MEAN_MODEL = ['temps/' TYPE '/mean_clean.mat'];

FIBER = ['temps/' TYPE '/lengths.mat'];
ORIENT = ['temps/' TYPE '/angles.mat'];
TRAIN_INPUT = ['include/input_' TYPE '.txt'];
datapath = textread(TRAIN_INPUT,'%s');

load(FIBER);
load(ORIENT);
load(INDEX_MODEL);
load(ORIEN_MODEL);
load(MEAN_MODEL);

dbstop if error
fileID = fopen('Total_results_bend.txt','w');
fprintf('Results are displayed in Total_results_bend.txt\n');
fileIDQV = fopen('fiberscore_data.csv','w');

disp('Processing images...')

j=1;
polm = 0;
fprintf(fileIDQV,'image_name|number_fibers|total_length|mean_length|polarity \n');

tic;
for i=1:length(datapath)
    fprintf('On cell number: %d\n',i)
    image_name = datapath{i};
    original= imread(image_name);    

    orig_eq = adapthisteq(original);
    BW = mean_clean{i} ~=0;   

    fprintf(fileIDQV,'%10s|', image_name);
    fibers_number = length(lengths{i});
    total_fiber_length(i) = sum(lengths{i});
    mean_fiber_length(i) = mean(lengths{i}); 
    total_fiberscore = sum(sum(orig_eq(BW==1)));
    total_mean = sum(mean_clean{i}(BW==1));
    Dx=sum(mean_clean{i}(BW==1).*cos(2*orient_clean{i}(BW==1)))/total_mean;
    Dy=sum(mean_clean{i}(BW==1).*sin(2*orient_clean{i}(BW==1)))/total_mean;
    
    polarity(i) = sqrt(Dx^2+Dy^2); 
    mean_angle = 1/2*atan(Dy/Dx)*180/pi;
    
    fprintf(fileIDQV,'%6.2f|',fibers_number);
    fprintf(fileIDQV,'%6.2f|',total_fiber_length(i));
    fprintf(fileIDQV,'%6.2f|',mean_fiber_length(i));
    fprintf(fileIDQV,'%.15f|\n',polarity(i));

    S = std(angles{i});
    fmt = '%3d          %4d         %6.2f          \n\';
end

fclose(fileIDQV);
fprintf('Took %g minutes to obtain the totals for the images.\n\n',toc/60);
