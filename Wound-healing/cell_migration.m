function [Res_stats, Res_colour, Res_gray, Res_Cells] = cellMigration(dataIn)
    % Validate input arguments
    if nargin < 1
        help cellMigration;
        [Res_stats, Res_colour, Res_gray, Res_Cells] = deal([]);
        return;
    end

    % Load or assign input image data
    b = loadImage(dataIn);

    % Normalize image if needed
    if min(b(:)) < 0
        b = b - min(b(:));
    end

    % Preprocess image with high-pass filtering
    bHP1 = preprocessImage(b);

    % Resize image
    reduced_bHP1 = resizeImage(bHP1);

    % Determine threshold for cell segmentation
    thres = determineThreshold(reduced_bHP1);

    % Segment cells in the image
    bCELLS = segmentCells(bHP1, thres);

    % Analyze cell orientation and morphology
    [angOrientation, SE1, SE2] = analyzeOrientation(bCELLS, b);

    % Process segmented cells
    [bCELLS3, bAreas4, numAreas] = processSegmentedCells(bCELLS, SE1, SE2, size(b));

    % Calculate cell statistics and common areas
    Res_stats = calculateStatistics(bCELLS3, bAreas4, numAreas, b);

    % Generate the output images
    [Res_colour, Res_gray, Res_Cells] = generateOutputImages(b, bCELLS3, Res_stats, numAreas, size(b));

    % Cleanup variables
    cleanupWorkspace();
end

% Helper Functions

function b = loadImage(dataIn)
    if isa(dataIn, 'char')
        b = imread(dataIn);
    else
        b = dataIn;
    end
end

function bHP1 = preprocessImage(b)
    [rows, cols, levs] = size(b);
    f1 = gaussF(rows, cols, 1, floor(rows / 50), floor(cols / 50), 1, floor(rows / 2), floor(cols / 2), 1);
    f1 = f1 / max(f1(:));
    f2 = gaussF(15, 15, 1);
    
    bf = fftshift(fft2(double(b)));
    if levs > 1
        bHP1 = processColorChannels(bf, f1, f2);
    else
        bHP1 = processSingleChannel(bf, f1, f2);
    end
    
    bHP1 = bHP1(20:end-19, 20:end-19);
end

function bHP1 = processColorChannels(bf, f1, f2)
    bHP1a = abs(ifft2(((1-f1) .* bf(:, :, 1))));
    bHP2a = abs(ifft2(((1-f1) .* bf(:, :, 2))));
    bHP3a = abs(ifft2(((1-f1) .* bf(:, :, 3))));
    bHP1 = conv2(bHP1a, f2, 'same') + conv2(bHP2a, f2, 'same') + conv2(bHP3a, f2, 'same');
    bHP1 = bHP1 - min(bHP1(:));
end

function bHP1 = processSingleChannel(bf, f1, f2)
    bHP1a = abs(ifft2(((1-f1) .* bf(:, :, 1))));
    bHP1 = conv2(bHP1a, f2, 'same');
    bHP1 = bHP1 - min(bHP1(:));
end

function reduced_bHP1 = resizeImage(bHP1)
    [rowsbHP1, colsbHP1] = size(bHP1);
    rowVector = linspace(1, rowsbHP1, 2^floor(log2(rowsbHP1)));
    colVector = linspace(1, colsbHP1, 2^floor(log2(colsbHP1)));
    reduced_bHP1 = interp2(bHP1, colVector, rowVector');
    reduced_bHP1 = normalizeImage(reduced_bHP1);
end

function normImg = normalizeImage(img)
    normImg = img - min(img(:));
    normImg = normImg / max(normImg(:));
end

function thres = determineThreshold(reduced_bHP1)
    for k = 3:6
        thres2(k-2) = graythresh(reduceu(reduced_bHP1, k));
    end
    thres = min(0.5 * (thres2 * max(reduced_bHP1(:)) + min(reduced_bHP1(:))));
end

function bCELLS = segmentCells(bHP1, thres)
    bCELLS(:, :, 1) = (bHP1 > thres);
end

function [angOrientation, SE1, SE2] = analyzeOrientation(bCELLS, b)
    angleRange = 0:10:179;
    TraceCells = traceTransform_old(bCELLS(1:4:end, 1:4:end), angleRange, [-1 1]);
    
    tracesStd = calculateTraceStd(TraceCells);
    [~, indVar] = max(tracesStd);
    angOrientation = angleRange(indVar);
    
    sizeOperator = max(10, round(min(size(b)) / 40));
    SE1 = strel('line', sizeOperator, angOrientation);
    SE2 = strel('line', sizeOperator, 90 + angOrientation);
end

function tracesStd = calculateTraceStd(TraceCells)
    tracesStd = zeros(1, 18);
    for k = 1:18
        ttt = squeeze(TraceCells(:, k, 9));
        ttt(squeeze(TraceCells(:, k, 11)) == 0) = [];
        tracesStd(k) = std(ttt);
    end
end

function [bCELLS3, bAreas4, numAreas] = processSegmentedCells(bCELLS, SE1, SE2, imgSize)
    bCELLS2 = imclose(bCELLS, SE2);
    bCELLS2 = imclose(bCELLS2, SE1);
    bCELLS2a = imfill(bCELLS2, 'holes');
    
    bCELLS2a = imdilate(bCELLS2a, ones(3));
    if sum(bCELLS2a(:)) > (0.95 * prod(imgSize))
        bCELLS3 = imopen(bCELLS2, SE2);
    else
        bCELLS3 = imopen(bCELLS2a, SE2);
    end
    
    bAreas1 = bwlabel(bCELLS3);
    bAreas2 = regionprops(bAreas1, 'area');
    bAreas3 = ismember(bAreas1, find(([bAreas2.Area] / prod(imgSize)) > 0.01));
    [bAreas4, numAreas] = bwlabel(bAreas3);
end

function Res_stats = calculateStatistics(bCELLS3, bAreas4, numAreas, imgSize)
    if numAreas == 1
        Res_stats = calculateSingleAreaStats(bCELLS3, bAreas4, imgSize);
    else
        Res_stats = calculateMultipleAreaStats(bCELLS3, bAreas4, imgSize);
    end
end

function Res_stats = calculateSingleAreaStats(bCELLS3, bAreas4, imgSize)
    bCELLS7a = bwmorph(bCELLS3 == max(bAreas4(:)), 'skel', Inf);
    complementArea = calculateComplementArea(bAreas4);
    
    totArea = sum(complementArea(:));
    relArea = totArea / prod(imgSize);
    
    Res_stats.area = [totArea; relArea];
end

function complementArea = calculateComplementArea(bAreas4)
    complementArea = bwlabel(1 - bAreas4);
    complementSizes = regionprops(complementArea, 'Area');
    [~, temp2] = max([complementSizes.Area]);
    complementArea = (complementArea == temp2);
end

function Res_stats = calculateMultipleAreaStats(bCELLS3, bAreas4, imgSize)
    % Skeletonize the largest two areas
    bCELLS7a = bwmorph(bCELLS3 == max(bAreas4(:)), 'skel', Inf);
    bCELLS7b = bwmorph(bCELLS3 == secondMax(bAreas4(:)), 'skel', Inf);
    
    % Calculate distance between skeleton points
    [rows, cols] = size(bCELLS7a);
    Res_stats = calculateDistanceStats(bCELLS7a, bCELLS7b, rows, cols);

    % Calculate common area between the two largest regions
    commonArea = calculateCommonArea(bCELLS7a, bCELLS7b);
    
    totArea = sum(commonArea(:));
    relArea = totArea / prod(imgSize);
    
    Res_stats.area = [totArea; relArea];
end

function Res_stats = calculateDistanceStats(bCELLS7a, bCELLS7b, rows, cols)
    [rowsA, colsA, rowsB, colsB] = extractSkeletonPoints(bCELLS7a, bCELLS7b, rows);
    
    distBetPoints = sqrt((rowsA - rowsB').^2 + (colsA - colsB').^2);
    minimumDist1 = min(distBetPoints);
    minimumDist2 = min(distBetPoints, [], 2);
    
    Res_stats.minimumDist = min(minimumDist1);
    Res_stats.maxDist = max([max(min(distBetPoints)) max(min(minimumDist2))]);
    Res_stats.avDist = mean([minimumDist1, minimumDist2']);
end

function [rowsA, colsA, rowsB, colsB] = extractSkeletonPoints(bCELLS7a, bCELLS7b, rows)
    xyCELLSa = find(bCELLS7a);
    rowsA = 1 + floor(xyCELLSa / rows);
    colsA = rem(xyCELLSa - 1, rows) + 1;
    
    xyCELLSb = find(bCELLS7b);
    rowsB = 1 + floor(xyCELLSb / rows);
    colsB = rem(xyCELLSb - 1, rows) + 1;
end

function commonArea = calculateCommonArea(bCELLS7a, bCELLS7b)
    Area2A = bwlabel(1 - (conv2(double(bCELLS7a), gaussF(3, 3, 1), 'same') > 0));
    Area2B = bwlabel(1 - (conv2(double(bCELLS7b), gaussF(3, 3, 1), 'same') > 0));
    
    commonArea = calculateCommonIntersection(Area2A, Area2B);
end

function commonArea = calculateCommonIntersection(Area2A, Area2B)
    A1_B1 = sum(sum((Area2A == 1) & (Area2B == 1)));
    A1_B2 = sum(sum((Area2A == 1) & (Area2B == 2)));
    A2_B1 = sum(sum((Area2A == 2) & (Area2B == 1)));
    A2_B2 = sum(sum((Area2A == 2) & (Area2B == 2)));
    
    if A1_B1 == 0
        commonArea = (Area2A == 2) & (Area2B == 2);
    elseif A1_B2 == 0
        commonArea = (Area2A == 2) & (Area2B == 1);
    elseif A2_B1 == 0
        commonArea = (Area2A == 1) & (Area2B == 2);
    elseif A2_B2 == 0
        commonArea = (Area2A == 1) & (Area2B == 1);
    end
end

function [Res_colour, Res_gray, Res_Cells] = generateOutputImages(b, bCELLS3, Res_stats, numAreas, imgSize)
    kernelDilation = ones(max(5, ceil(imgSize(1) / 200)));
    
    if size(b, 3) > 1
        [Res_colour, Res_gray] = generateColorImages(b, bCELLS3, Res_stats, numAreas, kernelDilation);
    else
        [Res_colour, Res_gray] = generateGrayImages(b, bCELLS3, Res_stats, numAreas, kernelDilation);
    end
    
    Res_Cells = storeCellStages(bCELLS3, numAreas, imgSize);
end

function [Res_colour, Res_gray] = generateColorImages(b, bCELLS3, Res_stats, numAreas, kernelDilation)
    if numAreas > 1
        finalBoundaries = uint8(imdilate(bCELLS3, kernelDilation));
    else
        finalBoundaries = uint8(imdilate(bCELLS3, kernelDilation));
    end
    
    finalBoundaries = repmat(finalBoundaries, [1 1 3]);
    Res_colour = 255 * finalBoundaries + (1 - finalBoundaries) .* b;
    Res_colour(:, :, 3) = Res_colour(:, :, 3) .* (1 + 0.5 * uint8(Res_stats.area > 0));
    Res_gray = double(sum(Res_colour(:, :, 2), 3)) .* (1 - 0.5 * (Res_stats.area > 0));
end

function [Res_colour, Res_gray] = generateGrayImages(b, bCELLS3, Res_stats, numAreas, kernelDilation)
    if numAreas > 1
        finalBoundaries = uint8(imdilate(bCELLS3, kernelDilation));
    else
        finalBoundaries = uint8(imdilate(bCELLS3, kernelDilation));
    end
    
    Res_colour = 255 * finalBoundaries + (1 - finalBoundaries) .* b;
    Res_colour(:, :, 2) = Res_colour(:, :, 1);
    Res_colour(:, :, 3) = 0.9 * Res_colour(:, :, 1) .* (1 + 0.5 * uint8(Res_stats.area > 0));
    Res_gray = double(sum(Res_colour(:, :, 2), 3)) .* (1 - 0.5 * (Res_stats.area > 0));
end

function Res_Cells = storeCellStages(bCELLS3, numAreas, imgSize)
    Res_Cells = cell(9, 1);
    Res_Cells{1} = bCELLS3;
    
    if numAreas > 1
        Res_Cells{10} = bCELLS3;
    end
end

function cleanupWorkspace()
    clearvars -except Res_stats Res_colour Res_gray Res_Cells;
end
