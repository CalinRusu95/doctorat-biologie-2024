function [redData] = reduceu(data, numReductions, reducein3D)
    % Check input arguments and set defaults
    if nargin < 1
        help reduceu;
        redData = [];
        return;
    end
    
    if nargin < 2 || isempty(numReductions)
        numReductions = 1;
    end
    
    if nargin < 3 || isempty(reducein3D)
        reducein3D = 0;
    end
    
    % Ensure data is in double format
    data = ensureDouble(data);
    
    % Recursively reduce data if multiple reductions are requested
    if numReductions > 1
        data = reduceu(data, numReductions - 1, reducein3D);
    end
    
    % Perform the reduction
    redData = performReduction(data, reducein3D);
end

% Helper Functions

function data = ensureDouble(data)
    if ~isa(data, 'double')
        data = double(data);
    end
end

function redData = performReduction(data, reducein3D)
    [rows, cols, levels] = size(data);

    if levels == 1 || reducein3D == 0
        % 2D reduction
        redData = convn(data, [1 1; 1 1]);
        redData = redData(2:2:end, 2:2:end, :) / 4;
    else
        % 3D reduction
        redData = convn(data, ones(2, 2, 2));
        redData = redData(2:2:end, 2:2:end, 2:2:end) / 8;
    end
end
