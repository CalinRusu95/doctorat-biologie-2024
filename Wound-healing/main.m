TYPE = 'cito'; 
disp('Reading input files...')
TRAIN_INPUT = ['input_' TYPE '.txt'];
datapath = textread(TRAIN_INPUT,'%s');
celldisp(datapath);
res = fopen('...\results.txt', 'w+');

for i=1:length(datapath)
    disp('Reading image...');
    image_name = datapath(i);

    dataIn=imread(string(image_name));
    [Res_stats,Res_colour,Res_gray]=cellMigration(dataIn);

    fprintf('On Cell: %s\n',datapath{i});
    disp(Res_stats);
    fprintf(res, '%f\n', Res_stats.avDist);

    figure('Name',string(image_name))
    imagesc(Res_colour)
    axis off;
end

fclose(res);