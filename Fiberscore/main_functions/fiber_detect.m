clc; clear;
TYPE = 'cito'; 
tic;

L = 6;      % Kernal size. Needs to be an even number so that there will be a center when building
K = 10;     % Angular resolution
TC = .55;   % threshold for correlation coefficient with a Gaussian profile
M = 0.3;    % threshold for NSD
N = 2.1;    % threshold for ratio of NSD between the perpendicular rods
T = .25;    % threshold for background-subracted intensity value for a fiber in the acquired image
w = 1;      % one pixel value


disp('Reading input files...')
TRAIN_INPUT = ['...\include\input_' TYPE '.txt'];
datapath = textread(TRAIN_INPUT,'%s');
  
index_cell = cell(size(datapath));
orientation_cell =cell(size(datapath));
mean_cell =cell(size(datapath)) ; 

for imgid=1:length(datapath)

    disp('Reading image...')
    image_path = ['...\cito\' datapath{imgid}];
    image_name = image_path;

    IMAGE = imread(image_name);
    [x,y,D] = size(IMAGE);

    if D==1
            original = IMAGE;
        elseif D >=3
            RGB = IMAGE;
        original = rgb2gray(RGB(:,:,1:3));
    end

    original = double(original)/max(max(double(original)));
    orig_eq = double(original);   
    modi_act = .5;
    [index,orientation,mean_fs] = fiberscore_conv(image_name,orig_eq,K,L,TC,M,N,T);


    index_cell{imgid} = index;
    orientation_cell{imgid} = orientation;
    mean_cell{imgid} = mean_fs;

    imwrite(index, strcat('...\temp_img\conv\',datapath{imgid}),'tif');

    param.K=K;
    param.L=L;
    param.w=w;
    param.TC=TC;
    param.M=M;
    param.N=N;
    param.T=T;

    save (['...\temps\' TYPE '\index_cell.mat'], 'index_cell')
    save (['...\temps\' TYPE '\orientation_cell.mat'], 'orientation_cell')
    save (['...\temps\' TYPE '\mean_cell.mat'], 'mean_cell')
    save (['...\temps\' TYPE '\param.mat'], 'param')
    fprintf('Took %g minutes to extract the fibers and find the orientation\n\n',toc/60);
end