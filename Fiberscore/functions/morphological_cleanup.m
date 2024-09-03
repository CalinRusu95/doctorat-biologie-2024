function [BW_final] = directional_morph_clean(index, LEN, gdir_thresh, ang_res)
    % Initialize the image and parameters
    img = index;
    [x, y] = size(img);
    
    % Define the range of angles for the structuring element
    ang = linspace(0, 180, ang_res); 
    D = length(ang);
    
    % Initialize a 3D array to store the results of the morphological operations
    union_set = zeros(x, y, D);
    
    % Perform directional morphological cleaning for each angle
    for i = 1:D
        union_set(:, :, i) = performMorphologicalOperation(img, LEN, ang(i));
    end
    
    % Calculate the gradient direction and threshold the result
    gdir = max(union_set, [], 3) - min(union_set, [], 3);
    BW_final = gdir > gdir_thresh;
    
    % Perform skeletonization and clean-up operations
    BW_final = cleanSkeleton(BW_final);
end

% Helper Functions

function IMG = performMorphologicalOperation(img, LEN, angle)
    % Create a linear structuring element at the given angle
    se = strel('line', LEN, angle);
    lambda = length(se.getneighbors);
    
    % Perform morphological operations based on the structuring element
    if lambda == length(se.getneighbors)
        IMG = imopen(img, se);
    else
        IMG = ordfilt2(ordfilt2(img, lambda - length(se.getneighbors) + 1, se.getnhood), lambda, se.getnhood);
    end
    
    % Minimize the result with the original image
    IMG = min(img, IMG);
end

function BW_final = cleanSkeleton(BW_final)
    % Perform skeletonization and other morphological cleaning operations
    BW_final = bwmorph(BW_final, 'skel', inf);
    BW_final = bwmorph(BW_final, 'hbreak', inf);
    BW_final = bwmorph(BW_final, 'spur', 5);
end
