function [index, orientation, mean_intensity] = fiberscore_conv(image_name, f, K, L, TC, M, N, T)
    % Initialize variables and parameters
    disp('Initializing image matrix and parameters...');
    [x, y] = size(f);
    [index, orientation, mean_intensity] = initializeOutputMatrices(x, y);
    
    CMP = zeros(x, y);
    CM = zeros(x, y);
    
    f = normalizeImage(f);
    P = padImage(f, L);
    
    [theta_k, theta_k_perp] = calculateAngles(K);
    y_I = createGaussianProfile(L);

    % Iterate through each angle to evaluate the image
    for l = 1:length(theta_k)
        evaluateAngle(image_name, theta_k, theta_k_perp, l, f, L, y_I, TC, M, N, T, ...
            CMP, CM, index, orientation, mean_intensity);
    end
end

% Helper Functions

function [index, orientation, mean_intensity] = initializeOutputMatrices(x, y)
    % Initialize the output matrices
    index = zeros(x, y);
    orientation = ones(x, y) * -100;
    mean_intensity = ones(x, y) * -1;
end

function f_normalized = normalizeImage(f)
    % Normalize the image to a range of [0, 1]
    f_normalized = double(f) / max(f(:));
end

function P = padImage(f, L)
    % Pad the image to accommodate kernel operations
    pad = L + L / 2 + 2;
    P = zeros(size(f, 1) + pad, size(f, 2) + pad);
    P(round(pad/2) + (1:size(f, 1)), round(pad/2) + (1:size(f, 2))) = f;
end

function [theta_k, theta_k_perp] = calculateAngles(K)
    % Calculate the angles for the main directions and their perpendiculars
    min_angle = 0; 
    max_angle = pi / (2 * K) * K; 
    theta_k = linspace(min_angle, max_angle, K+1); 
    theta_k_perp = theta_k + pi / 2;
end

function y_I = createGaussianProfile(L)
    % Create a Gaussian profile for the kernels
    std_dev = 2; 
    mean_profile = 0;
    xx = linspace(-6 + mean_profile, 6 + mean_profile, L + 1);
    y_I = exp(-0.5 * ((xx - mean_profile) / std_dev).^2);
end

function evaluateAngle(image_name, theta_k, theta_k_perp, l, f, L, y_I, TC, M, N, T, ...
    CMP, CM, index, orientation, mean_intensity)
    % Evaluate the image at a specific angle
    
    fprintf('Reading image: %s\n', image_name);
    fprintf('Evaluating the image with theta_k = %.2f degrees\n', theta_k(l) * 180 / pi);
    fprintf('Angles left to evaluate: %d\n', length(theta_k) - l);

    [col, row] = pix_displace(theta_k(l), (-L/2):(L/2));
    [col_perp, row_perp] = pix_displace(theta_k_perp(l), (-L/2):(L/2));

    % Create the convolution kernels
    [kernel, kernel_perp, kernel_y, kernel_y_perp] = createKernels(L, col, row, col_perp, row_perp, y_I);

    % Evaluate the image with the kernels
    [mean_fs, corr_fs] = evaluateKernels(f, kernel, kernel_y, L);
    [mean_fs_perp, corr_fs_perp] = evaluateKernels(f, kernel_perp, kernel_y_perp, L);
    
    % Identify good pixels based on the thresholds
    [idxgood, idxgoodperp] = identifyGoodPixels(corr_fs, nstd_fs, corr_fs_perp, nstd_fs_perp, mean_fs_perp, TC, M, N, T);
    
    % Update CMP and CM based on the identified pixels
    CMP(idxgood) = corr_fs{l}(idxgood); 
    CM(idxgoodperp) = corr_fs_perp{l}(idxgoodperp); 

    % Update the final index, orientation, and mean intensity
    updateFinalValues(index, orientation, mean_intensity, CM, CMP, idxgood, idxgoodperp, ...
        mean_fs, mean_fs_perp, theta_k, theta_k_perp, l);
end

function [kernel, kernel_perp, kernel_y, kernel_y_perp] = createKernels(L, col, row, col_perp, row_perp, y_I)
    % Create the convolution kernels for the given angles
    shift_center = round(length(col) / 2);
    kernel = createKernel(L, col, row, shift_center);
    kernel_perp = createKernel(L, col_perp, row_perp, shift_center);
    kernel_y = createKernelWithProfile(L, col, row, shift_center, y_I);
    kernel_y_perp = createKernelWithProfile(L, col_perp, row_perp, shift_center, y_I);
end

function kernel = createKernel(L, col, row, shift_center)
    kernel = zeros(L + 1);
    for i = 1:length(col)
        kernel(row(i) + shift_center, -col(i) + shift_center) = 1;
    end
end

function kernel_y = createKernelWithProfile(L, col, row, shift_center, y_I)
    kernel_y = zeros(L + 1);
    for i = 1:length(col)
        kernel_y(row(i) + shift_center, -col(i) + shift_center) = y_I(i);
    end
end

function [mean_fs, corr_fs] = evaluateKernels(f, kernel, kernel_y, L)
    term1 = imfilter(f, kernel_y, 'conv'); 
    term2 = imfilter(f, kernel, 'conv') * sum(kernel_y(:)) / (L + 1);
    term3 = (L + 1) / L * (imfilter(f .^ 2, kernel, 'conv') / (L + 1) - (imfilter(f, kernel, 'conv') / (L + 1)).^2);
    term4 = std(kernel_y);
    
    std_fs = stdfilt(f, kernel);
    mean_fs = imfilter(f, kernel, 'conv') / (L + 1);    
    corr_fs = (term1 - term2) ./ (sqrt(term3) * term4 * L);
end

function [idxgood, idxgoodperp] = identifyGoodPixels(corr_fs, nstd_fs, corr_fs_perp, nstd_fs_perp, mean_fs_perp, TC, M, N, T)
    idxgood = find(corr_fs > TC & nstd_fs > M & nstd_fs ./ nstd_fs_perp > N & mean_fs_perp > T);
    idxgoodperp = find(corr_fs_perp > TC & nstd_fs_perp > M & nstd_fs_perp ./ nstd_fs > N & mean_fs > T);
end

function updateFinalValues(index, orientation, mean_intensity, CM, CMP, idxgood, idxgoodperp, ...
    mean_fs, mean_fs_perp, theta_k, theta_k_perp, l)
    idxfinal = find(CM > CMP & CM > index); 
    idxfinal_perp = find(CMP > CM & CMP > index);
    
    index(idxfinal) = CM(idxfinal); 
    orientation(idxfinal) = theta_k(l);
    mean_intensity(idxfinal) = mean_fs(idxfinal);
    
    index(idxfinal_perp) = CMP(idxfinal_perp);
    orientation(idxfinal_perp) = theta_k_perp(l);
    mean_intensity(idxfinal_perp) = mean_fs_perp(idxfinal_perp);
end
