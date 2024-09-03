function dataPadded = padData(qtdata, numPadPixels, dimsToPad, padWith)
    % Set default padding value if not provided
    if nargin < 4
        padWith = 1;
    end
    
    % Determine the dimensions to pad
    if nargin < 3 || isempty(dimsToPad)
        [rows, cols, levs, numFeats] = size(qtdata);
    else
        [rows, cols, levs, numFeats] = parseDimsToPad(dimsToPad);
    end
    
    % Pad the rows and columns as needed
    qtdata = padRowsAndCols(qtdata, numPadPixels, padWith, rows, cols);
    
    % Pad the levels (depth) if needed
    if levs > 1
        qtdata = padLevels(qtdata, numPadPixels, padWith, levs);
    end
    
    % Return the padded data
    dataPadded = qtdata;
end

% Helper Functions

function [rows, cols, levs, numFeats] = parseDimsToPad(dimsToPad)
    % Extend dimsToPad to handle up to four dimensions
    dimsToPad = [dimsToPad ones(1, 4 - numel(dimsToPad))];
    rows = dimsToPad(1);
    cols = dimsToPad(2);
    levs = dimsToPad(3);
    numFeats = dimsToPad(4);
end

function paddedData = padRowsAndCols(data, numPadPixels, padWith, rows, cols)
    % Pad rows if more than one row is present
    if rows > 1
        padTop = padWith * repmat(data(1, :, :, :), [numPadPixels, 1, 1, 1]);
        padBottom = padWith * repmat(data(end, :, :, :), [numPadPixels, 1, 1, 1]);
        data = [padTop; data; padBottom];
    end
    
    % Pad columns if more than one column is present
    if cols > 1
        padLeft = padWith * repmat(data(:, 1, :, :), [1, numPadPixels, 1, 1]);
        padRight = padWith * repmat(data(:, end, :, :), [1, numPadPixels, 1, 1]);
        data = [padLeft, data, padRight];
    end
    
    paddedData = data;
end

function paddedData = padLevels(data, numPadPixels, padWith, levs)
    % Initialize the padded data array
    [rows, cols, ~, numFeats] = size(data);
    paddedData = zeros(rows, cols, levs + 2 * numPadPixels, numFeats);
    
    % Copy the original data into the middle of the padded array
    paddedData(:, :, numPadPixels+1:numPadPixels+levs, :) = data;
    
    % Pad the front and back of the levels dimension
    padFront = padWith * repmat(data(:, :, 1, :), [1, 1, numPadPixels, 1]);
    padBack = padWith * repmat(data(:, :, end, :), [1, 1, numPadPixels, 1]);
    paddedData(:, :, 1:numPadPixels, :) = padFront;
    paddedData(:, :, numPadPixels+levs+1:end, :) = padBack;
end
