function [fiber_length, orient, num, max_num_fiber] = angle_length_extract(orient_clean, mean_clean)
    % Initialize output variables
    fiber_length = [];
    max_num_fiber = 0;

    % Create a binary mask where mean_clean is non-zero
    BW = mean_clean ~= 0;

    % Label connected components in the binary mask
    [L, num] = bwlabel(BW, 8);
    
    % Initialize length_label matrix
    [N, M] = size(mean_clean);
    length_label = zeros(N, M);

    % Extract fiber lengths and find the maximum fiber length
    fiber_length = extractFiberLengths(L, num, length_label, max_num_fiber);

    % Extract the orientations of the fibers
    orient = extractOrientations(orient_clean, BW);
end

% Helper Functions

function fiber_length = extractFiberLengths(L, num, length_label, max_num_fiber)
    fiber_length = [];
    ii = 1;

    for i = 1:num
        len = sum(L(:) == i);
        id = find(L == i);
        length_label(id) = len;

        if len >= 0
            fiber_length(ii) = len;
            ii = ii + 1;
        end

        if len > max_num_fiber
            max_num_fiber = len;
        end
    end
end

function orient = extractOrientations(orient_clean, BW)
    orient = orient_clean(BW);
    orient = orient(orient ~= -100) * 180 / pi;
end
