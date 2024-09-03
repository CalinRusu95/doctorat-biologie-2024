function bor = zerocross(imag)
    % Constants and initial setup
    delta = 1e-5;
    [lins, cols, levels] = size(imag);
    
    % Convert image to logical sign matrix if necessary
    imag = convertToSign(imag);
    
    % Determine the type of zero-crossing detection based on dimensions
    if isVector(lins, cols, levels)
        bor = detectZerocrossVector(imag, cols, lins, delta);
    elseif isSingleColumnOrRowVector(lins, cols, levels)
        bor = detectZerocrossSingleVector(imag, delta);
    elseif is2DImage(lins, cols, levels)
        bor = detectZerocross2D(imag, lins, cols);
    else
        bor = detectZerocross3D(imag, lins, cols, levels, delta);
    end
end

% Helper Functions

function signImg = convertToSign(imag)
    if ~isa(imag, 'logical')
        signImg = sign(imag);
    else
        signImg = imag;
    end
end

function isVec = isVector(lins, cols, levels)
    isVec = (cols == 1 || lins == 1) && levels == 1;
end

function isSingleVec = isSingleColumnOrRowVector(lins, cols, levels)
    isSingleVec = (cols == 1 && lins == 1 && levels ~= 1);
end

function is2D = is2DImage(lins, cols, levels)
    is2D = lins ~= 1 && cols ~= 1 && levels == 1;
end

function bor = detectZerocrossVector(imag, cols, lins, delta)
    if cols >= lins
        yy = [0 imag(1:cols-1)];
    else
        yy = [0 imag(1:lins-1)']';
    end
    bor = abs(sign(imag + delta) - sign(yy + delta));
end

function bor = detectZerocrossSingleVector(imag, delta)
    imag = permute(imag, [2 3 1]);
    yy = [0 imag(1:end-1)];
    bor = sign(imag + delta) - sign(yy + delta);
end

function bor = detectZerocross2D(imag, lins, cols)
    diffVer = diff(imag, 1, 1);
    diffHor = diff(imag, 1, 2);
    
    qq1 = [zeros(1, cols); diffVer > 0];
    qq2 = [diffVer < 0; zeros(1, cols)];
    qq3 = [(diffHor < 0) zeros(lins, 1)];
    qq4 = [zeros(lins, 1) (diffHor > 0)];
    
    bor = qq1 | qq2 | qq3 | qq4;
end

function bor = detectZerocross3D(imag, lins, cols, levels, delta)
    yy5 = [zeros(1, cols, levels); imag(1:lins-1, :, :)];
    yy6 = [imag(2:lins, :, :); zeros(1, cols, levels)];
    yy7 = [zeros(lins, 1, levels) imag(:, 1:cols-1, :)];
    yy8 = [imag(:, 2:cols, :) zeros(lins, 1, levels)];
    
    bor5 = calculateBorder(imag, yy5, delta);
    bor6 = calculateBorder(imag, yy6, delta);
    bor7 = calculateBorder(imag, yy7, delta);
    bor8 = calculateBorder(imag, yy8, delta);
    
    bor = sign(bor5 + bor6 + bor7 + bor8);
end

function bor = calculateBorder(imag, yy, delta)
    bor = fix(delta + (sign(imag) - sign(yy)) / 2);
end
