function [gauss] = gaussF(rowDim, colDim, levDim, rowSigma, colSigma, levSigma, rowMiu, colMiu, levMiu, rho)
    % Validate and set default input arguments
    if nargin < 1
        help gaussF;
        gauss = [];
        return;
    end
    
    if nargin < 10
        rho = 0;
    end
    
    [rowDim, colDim, levDim] = parseDimensions(rowDim, colDim, levDim);
    
    % Create meshgrid for Gaussian calculation
    [rr, cc, dd] = createMeshGrid(rowDim, colDim, levDim);
    
    % Set default values for means and sigmas if not provided
    [rowMiu, colMiu, levMiu] = setDefaultMeans(nargin, rr, cc, dd, rowMiu, colMiu, levMiu);
    [rowSigma, colSigma, levSigma] = setDefaultSigmas(nargin, rowMiu, colMiu, levMiu, rowSigma, colSigma, levSigma);
    
    % Ensure sigma values are positive
    [rowSigma, colSigma, levSigma] = ensurePositiveSigmas(rowSigma, colSigma, levSigma);
    
    % Calculate Gaussian filter based on input arguments
    if ismatrix(rho)
        gauss = calculateCorrelatedGaussian(rr, cc, dd, rowSigma, colSigma, levSigma, rowMiu, colMiu, levMiu, rho);
    else
        gauss = calculateStandardGaussian(rr, cc, dd, rowSigma, colSigma, levSigma, rowMiu, colMiu, levMiu, rho, rowDim, colDim, levDim);
    end
    
    % Replace NaN values with zeros
    gauss(isnan(gauss)) = 0;
end

% Helper Functions

function [rowDim, colDim, levDim] = parseDimensions(rowDim, colDim, levDim)
    if nargin == 1
        [wRow, wCol, ~] = size(rowDim);
        if wCol == 3
            levDim = rowDim(3);
            colDim = rowDim(2);
            rowDim = rowDim(1);
        elseif wCol == 2
            colDim = rowDim(2);
            rowDim = rowDim(1);
            levDim = 1;
        elseif wCol == 1
            colDim = 1;
            levDim = 1;
        end
    elseif nargin == 2
        levDim = 1;
    end
end

function [rr, cc, dd] = createMeshGrid(rowDim, colDim, levDim)
    [rr, cc, dd] = meshgrid(1:ceil(rowDim), 1:ceil(colDim), 1:ceil(levDim));
end

function [rowMiu, colMiu, levMiu] = setDefaultMeans(nargin, rr, cc, dd, rowMiu, colMiu, levMiu)
    if nargin <= 6
        rowMiu = sum(rr(:)) / numel(rr);
        colMiu = sum(cc(:)) / numel(cc);
        levMiu = sum(dd(:)) / numel(dd);
    end
end

function [rowSigma, colSigma, levSigma] = setDefaultSigmas(nargin, rowMiu, colMiu, levMiu, rowSigma, colSigma, levSigma)
    sigmaFactor = 1.1774;
    if nargin <= 3
        rowSigma = (rowMiu - 1) / sigmaFactor;
        colSigma = (colMiu - 1) / sigmaFactor;
        levSigma = (levMiu - 1) / sigmaFactor;
    end
end

function [rowSigma, colSigma, levSigma] = ensurePositiveSigmas(rowSigma, colSigma, levSigma)
    rowSigma = max(rowSigma, 1e-6);
    colSigma = max(colSigma, 1e-6);
    levSigma = max(levSigma, 1e-6);
end

function gauss = calculateCorrelatedGaussian(rr, cc, dd, rowSigma, colSigma, levSigma, rowMiu, colMiu, levMiu, rho)
    invSigma = inv(rho);
    [Srr, Scc, Sdd, Src, Srd, Scd] = extractSigmaComponents(invSigma, size(rho, 1));
    
    exp_r = (1 / rowSigma^2) * (rr - rowMiu).^2;
    exp_c = (1 / colSigma^2) * (cc - colMiu).^2;
    exp_d = (1 / levSigma^2) * (dd - levMiu).^2;
    exp_rc = (1 / (rowSigma * colSigma)) * (rr - rowMiu) .* (cc - colMiu);
    exp_rd = (1 / (rowSigma * levSigma)) * (rr - rowMiu) .* (dd - levMiu);
    exp_cd = (1 / (levSigma * colSigma)) * (dd - levMiu) .* (cc - colMiu);
    
    gauss = exp(-(Srr * exp_r + Scc * exp_c + Sdd * exp_d + Src * exp_rc + Srd * exp_rd + Scd * exp_cd));
end

function [Srr, Scc, Sdd, Src, Srd, Scd] = extractSigmaComponents(invSigma, sizeRho)
    Srr = invSigma(1, 1);
    Scc = invSigma(2, 2);
    if sizeRho == 2
        Src = 2 * invSigma(2, 1);
        Srd = 0;
        Scd = 0;
        Sdd = 0;
    else
        Src = 2 * invSigma(2, 1);
        Srd = 2 * invSigma(1, 3);
        Scd = 2 * invSigma(2, 3);
        Sdd = invSigma(3, 3);
    end
end

function gauss = calculateStandardGaussian(rr, cc, dd, rowSigma, colSigma, levSigma, rowMiu, colMiu, levMiu, rho, rowDim, colDim, levDim)
    rho = clamp(rho, -0.999999, 0.999999);
    
    filter.x2 = create1DGaussianFilter(rr, rowSigma, rowMiu);
    filter.y2 = create1DGaussianFilter(cc, colSigma, colMiu);
    filter.z2 = create1DGaussianFilter(dd, levSigma, levMiu);
    
    if colDim > 1 && rowDim > 1
        gauss = create2DOr3DGaussian(filter, rho, colDim, rowDim, levDim);
    else
        gauss = select1DGaussian(filter);
    end
end

function rho = clamp(value, minVal, maxVal)
    rho = max(min(value, maxVal), minVal);
end

function filter1D = create1DGaussianFilter(dim, sigma, miu)
    filter1D = (1 / (sqrt(2 * pi) * sigma)) * exp(-((dim - miu).^2) / (2 * sigma^2));
    filter1D = filter1D / sum(filter1D);
end

function gauss = create2DOr3DGaussian(filter, rho, colDim, rowDim, levDim)
    rhoExponent = (-(rho * (filter.x2 - mean(filter.x2))' * (filter.y2 - mean(filter.y2))) / (std(filter.x2) * std(filter.y2)));
    filter.rho = (1 / sqrt(1 - rho^2)) * exp(rhoExponent);
    
    twoDFilter = (filter.x2' * filter.y2) .* filter.rho;
    if ceil(levDim) > 1
        gauss = create3DGaussian(twoDFilter, filter.z2, levDim);
    else
        gauss = twoDFilter;
    end
end

function gauss3D = create3DGaussian(twoDFilter, filterZ, levDim)
    gauss3D = zeros([size(twoDFilter), ceil(levDim)]);
    for ii = 1:ceil(levDim)
        gauss3D(:, :, ii) = twoDFilter .* filterZ(ii);
    end
end

function gauss = select1DGaussian(filter)
    if length(filter.x2) > length(filter.y2)
        gauss = filter.x2;
    else
        gauss = filter.y2;
    end
end
