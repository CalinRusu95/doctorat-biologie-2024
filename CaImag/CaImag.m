clear all;
close all;

% Open result files for writing
res = fopen('...\results.txt', 'w+');
specs = fopen('...\spectre.txt', 'w+');

% Define the input path and read the list of files
InputPath = '...\Input\...\';
files = getFilesList(InputPath);
FilesNum = numel(files);
cplot = ceil(sqrt(FilesNum));

% Initialize storage for signals and results
corrRes = [];

% Process each file
for i = 1:FilesNum
    % Process the file and extract parameters
    [param, corrRes] = processFile(i, files(i), InputPath, corrRes);

    % Save results to the result file
    writeResultsToFile(res, files(i).name, param);
    
    % Plot results
    plotResults(i, FilesNum, cplot, files(i).name, corrRes);
end

% Close the result files
fclose(res);
fclose(specs);

% Helper Functions

function files = getFilesList(InputPath)
    cd(InputPath);
    files = dir('**');
    files(1:2) = [];  % Remove '.' and '..' entries
end

function [param, corrRes] = processFile(i, file, InputPath, corrRes)
    % Load the file data
    FileAdrs = fullfile(InputPath, file.name);
    opts = detectImportOptions(FileAdrs);
    A = readtable(FileAdrs, opts);
    
    % Extract signal data
    X = table2array(A(:, 1)); 
    spectra = table2array(A(:, 2:end)); 
    spectra = smoothSpectra(spectra);

    % Analyze the spectrum and extract peaks
    [tmax, peakmax, vmax, vmin, tini, tfin, aria, pks, locs] = analyzeSpectrum(X, spectra);
    
    % Calculate additional parameters
    [param, corrRes] = calculateParameters(i, X, tmax, tini, tfin, peakmax, vmax, vmin, pks, locs, aria, corrRes);
end

function spectra = smoothSpectra(spectra)
    % Smooth the spectra using Savitzky-Golay filter
    for j = 1:size(spectra, 2)
        spectra(:, j) = smooth(spectra(:, j), 11, 'sgolay');
    end
end

function [tmax, peakmax, vmax, vmin, tini, tfin, aria, pks, locs] = analyzeSpectrum(X, spectra)
    % Initialize parameters
    fracini = 0.2;
    fracfin = 0.2;
    pkfrac = 0.05;
    promfrac = 0.1;
    MPD = 2000;

    % Analyze each column of the spectrum
    for j = 1:size(spectra, 2)
        Yout = spectra(:, j);
        intd = gradient(Yout);  
        [c, l] = max(Yout); 
        
        % Extract peak-related data
        [tini(j), tfin(j), aria(j), pks{j}, locs{j}] = extractPeakData(X, Yout, c, intd, l, fracini, fracfin, pkfrac, promfrac, MPD);
        
        % Calculate additional parameters
        tmax(j) = X(l);
        peakmax(j) = c;
        vmax(j) = max(intd);
        vmin(j) = mean(intd(l:end));
    end
end

function [tini, tfin, aria, pks, locs] = extractPeakData(X, Yout, c, intd, l, fracini, fracfin, pkfrac, promfrac, MPD)
    % Initialize variables
    tini = NaN;
    tfin = NaN;
    
    % Find the initial and final times of the peak
    for k = 1:length(Yout) - 1
        if Yout(k + 1) > fracini * c && isnan(tini)
            tini = (X(k) + X(k + 1)) / 2;
        end
        
        if X(k) > X(l) && Yout(k + 1) < fracfin * c && isnan(tfin)
            tfin = (X(k) + X(k + 1)) / 2;
        end
    end
    
    % Extract peak area and find peaks
    xpeak = X(find(X == tini):find(X == tfin));
    peak = Yout(find(X == tini):find(X == tfin));
    aria = trapz(xpeak, peak);
    [pks, locs] = findpeaks(peak, xpeak, 'MinPeakHeight', pkfrac * c, 'MinPeakDistance', MPD, 'MinPeakProminence', promfrac * c);
end

function [param, corrRes] = calculateParameters(i, X, tmax, tini, tfin, peakmax, vmax, vmin, pks, locs, aria, corrRes)
    % Calculate temporal and amplitude-related parameters
    Tincrease = tmax - tini;
    Tdekay = tfin - tmax;
    Ttotal = tfin - tini; 
    Latency = tini - 360966;
    Asymmetry = Tdekay ./ Tincrease; 
    Increase1Speed = vmax;
    Decrease1Speed = abs(vmin);
    
    % Determine peak characteristics
    PeakNumber = length(pks);
    Peak1Amplitude = peakmax; 
    Peak2Amplitude = 0;
    PeaksDistance = 0;
    Type = 1;
    
    if PeakNumber == 2
        [Type, Peak1Amplitude, Peak2Amplitude, PeaksDistance, Increase2Speed, Decrease2Speed] = analyzePeakTypes(PeakNumber, pks, locs, intd, X);
    else
        Increase2Speed = 0;
        Decrease2Speed = 0;
    end
    
    % Store calculated parameters
    param = [Type, Ttotal, Latency, Asymmetry, Increase1Speed, Increase2Speed, Decrease1Speed, Decrease2Speed, aria, Peak1Amplitude, Peak2Amplitude, PeaksDistance];
    corrRes(i, :) = [Ttotal, Latency, Asymmetry, Increase1Speed, Decrease1Speed, aria, Peak1Amplitude];
end

function [Type, Peak1Amplitude, Peak2Amplitude, PeaksDistance, Increase2Speed, Decrease2Speed] = analyzePeakTypes(PeakNumber, pks, locs, intd, X)
    % Determine peak type and calculate related parameters
    if pks(1) > pks(2)
        Type = 2;
        PeaksDistance = abs(locs(1) - locs(2));
        Peak1Amplitude = pks(1);
        Peak2Amplitude = pks(2);
    else
        Type = 3;
        PeaksDistance = abs(locs(1) - locs(2));
        Peak1Amplitude = pks(1);
        Peak2Amplitude = pks(2); 
    end
    
    % Calculate secondary increase and decrease speeds
    gradis = sort(intd);
    Increase2Speed = gradis(end - 1);
    dec_index = find(X == locs(2));
    Decrease2Speed = mean(intd(dec_index:end));
end

function writeResultsToFile(res, name, param)
    % Write the calculated parameters to the result file
    fprintf(res, '%9s %3.0f %10.1f %10.1f %15.6f %15.6f %15.6f %15.6f %15.6f %12.4f %11.4f %11.4f %13.0f \n', name, param);
end

function plotResults(i, FilesNum, cplot, name, corrRes)
    % Plot the spectrum and peak results
    figure(FilesNum + 2)
    subplot(cplot, cplot, i)
    plot(X, spectra)
    title(name)
    
    figure(FilesNum + 3)
    subplot(cplot, cplot, i)
    plot(X, spectra)
    hold on
    plot(xpeak, peak)
    title(name)
    hold off
    
    % Plot the correlation matrix
    if i == FilesNum
        figure(FilesNum + 4)
        plotCorrelationMatrix(corrRes);
    end
end

function plotCorrelationMatrix(corrRes)
    % Plot the correlation matrix of the results
    [H, AX, BigAx, P, PAx] = plotmatrix(corrRes);
    set(gca, 'xlim', [0 5])
    labels = {'Ttotal', 'Latency', 'Asymmetry', 'Inc. Speed', 'Dec. Speed', 'Area', 'Amplitude'};
    for idx = 1:7
        ylabel(AX(idx, 1), labels{idx}, 'FontSize', 10, 'Rotation', 90);
        xlabel(AX(7, idx), labels{idx}, 'FontSize', 10, 'Rotation', 0);
    end
end
