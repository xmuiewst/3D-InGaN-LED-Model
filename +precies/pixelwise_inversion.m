function varargout = pixelwise_inversion(action, varargin)

    switch lower(string(action))
        case "run"
            [varargout{1:nargout}] = runInversion(varargin{:});
        case "metrics"
            [varargout{1:nargout}] = calculateMetrics(varargin{:});
        case "defaultoptions"
            varargout{1} = defaultOptions();
        case "selftest"
            varargout{1} = selfTest();
        otherwise
            error('precies:pixelwise_inversion:InvalidAction', ...
                'Unsupported action: %s', action);
    end
end

function result = runInversion(experimentalData, baselineSimulationData, options)
    if nargin < 3
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);

    [wavelengthNm, expCube, simCube] = alignInputData(experimentalData, baselineSimulationData, opts);
    [expProcMat, expNormMat, rows, cols] = preprocessCubeToMatrix(expCube, opts);
    [simProcMat, simNormMat] = preprocessCubeToMatrix(simCube, opts);
    validPixelMask = buildValidPixelMask(expProcMat, simProcMat, rows, cols, opts);

    initialMetrics = calculateMetricsFromMatrices( ...
        wavelengthNm, expProcMat, expNormMat, simProcMat, simNormMat, rows, cols, opts, validPixelMask);
    initialMetrics.voltage = opts.voltage;
    initialMetrics.summaryTable = buildSummaryTable(initialMetrics);
    [inversionPixelMask, selectionInfo] = selectPixelsForInversion(initialMetrics, validPixelMask, opts);

    nPix = rows * cols;
    fittedProcMat = simProcMat;
    parameterValues = nan(nPix, 7);
    objectiveMap = nan(nPix, 1);
    exitFlagMap = zeros(nPix, 1);
    iterationMap = zeros(nPix, 1);

    bounds = parameterBounds(opts);
    optimOptions = optimset( ...
        'Display', opts.optimizerDisplay, ...
        'MaxIter', opts.maxIterations, ...
        'MaxFunEvals', opts.maxFunctionEvaluations, ...
        'TolX', opts.toleranceX, ...
        'TolFun', opts.toleranceFun);

    inversionIndices = find(inversionPixelMask(:));
    for progressIdx = 1:numel(inversionIndices)
        pixelIdx = inversionIndices(progressIdx);

        targetProc = expProcMat(pixelIdx, :)';
        targetNorm = expNormMat(pixelIdx, :)';
        baselineProc = simProcMat(pixelIdx, :)';
        p0 = initialPixelParameters(wavelengthNm, targetProc, baselineProc, opts, bounds);

        if opts.enableOptimization
            u0 = paramsToUnconstrained(p0, bounds);
            objective = @(u) pixelObjective( ...
                unconstrainedToParams(u, bounds), p0, wavelengthNm, targetNorm, targetProc, baselineProc, opts);
            [uBest, objectiveValue, exitFlag, optimOutput] = fminsearch(objective, u0, optimOptions);
            pBest = unconstrainedToParams(uBest, bounds);
            iterationCount = getFieldOrDefault(optimOutput, 'iterations', 0);
        else
            pBest = p0;
            objectiveValue = pixelObjective(pBest, p0, wavelengthNm, targetNorm, targetProc, baselineProc, opts);
            exitFlag = 1;
            iterationCount = 0;
        end

        fittedProcMat(pixelIdx, :) = synthesizeCorrectedSpectrum(wavelengthNm, baselineProc, pBest, opts)';
        parameterValues(pixelIdx, :) = pBest(:)';
        objectiveMap(pixelIdx) = objectiveValue;
        exitFlagMap(pixelIdx) = exitFlag;
        iterationMap(pixelIdx) = iterationCount;

        if opts.progressEveryPixels > 0 && ...
                (mod(progressIdx, opts.progressEveryPixels) == 0 || progressIdx == numel(inversionIndices))
            fprintf('pixelwise inversion: %d/%d selected pixels\n', progressIdx, numel(inversionIndices));
        end
    end

    fittedNormMat = normalizeRowsByMax(fittedProcMat);
    result = calculateMetricsFromMatrices( ...
        wavelengthNm, expProcMat, expNormMat, fittedProcMat, fittedNormMat, rows, cols, opts, validPixelMask);

    result.voltage = opts.voltage;
    result.initialMetrics = initialMetrics;
    result.parameterMaps = buildParameterMaps(parameterValues, objectiveMap, exitFlagMap, iterationMap, rows, cols);
    result.validPixelMask = reshape(validPixelMask, rows, cols);
    result.inversionPixelMask = reshape(inversionPixelMask, rows, cols);
    result.selectionInfo = selectionInfo;
    result.fittedData = buildSpectraData(wavelengthNm, fittedProcMat, rows, cols);
    result.baselineData = buildSpectraData(wavelengthNm, simProcMat, rows, cols);
    result.options = opts;
    result.summaryTable = buildSummaryTable(result);

    if opts.writeOutputs && ~isempty(opts.outputDirectory)
        writeInversionOutputs(result, opts);
    end
end

function result = calculateMetrics(experimentalData, simulationData, options)
    if nargin < 3
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);
    [wavelengthNm, expCube, simCube] = alignInputData(experimentalData, simulationData, opts);
    [expProcMat, expNormMat, rows, cols] = preprocessCubeToMatrix(expCube, opts);
    [simProcMat, simNormMat] = preprocessCubeToMatrix(simCube, opts);
    validPixelMask = buildValidPixelMask(expProcMat, simProcMat, rows, cols, opts);
    result = calculateMetricsFromMatrices( ...
        wavelengthNm, expProcMat, expNormMat, simProcMat, simNormMat, rows, cols, opts, validPixelMask);
    result.voltage = opts.voltage;
    result.summaryTable = buildSummaryTable(result);
end

function opts = defaultOptions()
    opts = struct();
    opts.ganWindow = [350, 400];
    opts.mqwWindow = [400, 550];
    opts.defaultWavelengthRange = [300, 800];
    opts.wavelengthAxis = [];
    opts.voltage = NaN;
    opts.roi = [];
    opts.smoothingWindow = 5;
    opts.minPeakRelativeHeight = 0.03;
    opts.minPixelSignal = 1e-9;
    opts.softWindowEdgeNm = 8;
    opts.enableOptimization = true;
    opts.maxIterations = 90;
    opts.maxFunctionEvaluations = 360;
    opts.toleranceX = 1e-4;
    opts.toleranceFun = 1e-5;
    opts.optimizerDisplay = 'off';
    opts.progressEveryPixels = 0;
    opts.writeOutputs = false;
    opts.outputDirectory = '';
    opts.storeResidualCube = false;
    opts.autoSelectBadPixels = true;
    opts.badPixelNccThreshold = 0.985;
    opts.badPixelNRMSEThreshold = 0.060;
    opts.badPixelPeakShiftThresholdNm = 6;
    opts.badPixelFwhmShiftThresholdNm = 14;
    opts.badPixelRatioErrorThresholdPercent = 25;
    opts.badPixelMaxFraction = 0.45;
    opts.badPixelMinCount = 1;
    opts.objectiveNccWeight = 0.18;
    opts.objectivePeakWeight = 0.08;
    opts.objectiveFwhmWeight = 0.04;
    opts.objectiveRatioWeight = 0.03;
    opts.regularizationWeight = 0.010;
    opts.peakToleranceNm = 8;
    opts.fwhmToleranceNm = 18;
    opts.ratioTolerancePercent = 35;
    opts.parameterBounds = struct( ...
        'ganGain', [0.25, 4.0], ...
        'mqwGain', [0.25, 4.0], ...
        'backgroundGain', [0.35, 2.8], ...
        'ganShiftNm', [-18, 18], ...
        'mqwShiftNm', [-18, 18], ...
        'ganWidthScale', [0.65, 1.75], ...
        'mqwWidthScale', [0.65, 1.75]);
end

function [wavelengthNm, cube1, cube2] = alignInputData(data1, data2, opts)
    [cube1Raw, wl1] = dataToCube(data1, opts, []);
    [cube2Raw, wl2] = dataToCube(data2, opts, wl1);

    targetRows = min(size(cube1Raw, 1), size(cube2Raw, 1));
    targetCols = min(size(cube1Raw, 2), size(cube2Raw, 2));
    cube1Raw = cube1Raw(1:targetRows, 1:targetCols, :);
    cube2Raw = cube2Raw(1:targetRows, 1:targetCols, :);

    if ~isempty(opts.wavelengthAxis)
        wavelengthNm = double(opts.wavelengthAxis(:));
    else
        wavelengthNm = clipWavelengthAxis(wl1(:), opts);
    end
    cube1 = interpSpectralCube(cube1Raw, wl1, wavelengthNm);
    cube2 = interpSpectralCube(cube2Raw, wl2, wavelengthNm);
end

function [cube, wavelengthNm] = dataToCube(data, opts, fallbackAxis)
    if isstruct(data) && isfield(data, 'wavelengthIntensityMaps') && ~isempty(data.wavelengthIntensityMaps)
        cube = double(data.wavelengthIntensityMaps);
        if isfield(data, 'wavelengthAxis') && ~isempty(data.wavelengthAxis)
            wavelengthNm = double(data.wavelengthAxis(:));
        else
            wavelengthNm = defaultAxisForCube(size(cube, 3), opts);
        end
        return;
    end

    if isstruct(data) && isfield(data, 'totalSpectra') && ~isempty(data.totalSpectra)
        [cube, wavelengthNm] = spectraCellToCube(data.totalSpectra, opts, fallbackAxis);
        return;
    end

    if iscell(data)
        [cube, wavelengthNm] = spectraCellToCube(data, opts, fallbackAxis);
        return;
    end

    if isnumeric(data) && ndims(data) == 3
        cube = double(data);
        if ~isempty(opts.wavelengthAxis)
            wavelengthNm = double(opts.wavelengthAxis(:));
        elseif ~isempty(fallbackAxis)
            wavelengthNm = double(fallbackAxis(:));
        else
            wavelengthNm = defaultAxisForCube(size(cube, 3), opts);
        end
        return;
    end

    error('precies:pixelwise_inversion:InvalidData', ...
        'Input data must contain wavelengthIntensityMaps, totalSpectra, a cell array of spectra, or a 3-D numeric cube.');
end

function [cube, wavelengthNm] = spectraCellToCube(spectraCell, opts, fallbackAxis)
    [rows, cols] = size(spectraCell);
    if ~isempty(opts.wavelengthAxis)
        wavelengthNm = double(opts.wavelengthAxis(:));
    elseif ~isempty(fallbackAxis)
        wavelengthNm = double(fallbackAxis(:));
    else
        wavelengthNm = [];
        for idx = 1:numel(spectraCell)
            spectrum = spectraCell{idx};
            if ~isempty(spectrum) && size(spectrum, 2) >= 2
                wavelengthNm = clipWavelengthAxis(double(spectrum(:, 1)), opts);
                break;
            end
        end
        if isempty(wavelengthNm)
            wavelengthNm = defaultAxisForCube(512, opts);
        end
    end

    cube = zeros(rows, cols, numel(wavelengthNm));
    for idx = 1:numel(spectraCell)
        spectrum = spectraCell{idx};
        if isempty(spectrum) || size(spectrum, 2) < 2
            continue;
        end
        [r, c] = ind2sub([rows, cols], idx);
        cube(r, c, :) = interpolateSpectrumToAxis(spectrum, wavelengthNm);
    end
end

function wavelengthNm = defaultAxisForCube(nWavelengths, opts)
    wavelengthNm = linspace(opts.defaultWavelengthRange(1), opts.defaultWavelengthRange(2), nWavelengths)';
end

function axisOut = clipWavelengthAxis(axisIn, opts)
    axisOut = double(axisIn(:));
    axisOut = axisOut(isfinite(axisOut));
    if isempty(axisOut)
        axisOut = defaultAxisForCube(512, opts);
        return;
    end
    mask = axisOut >= opts.defaultWavelengthRange(1) & axisOut <= opts.defaultWavelengthRange(2);
    if nnz(mask) >= 16
        axisOut = axisOut(mask);
    end
    axisOut = unique(axisOut, 'stable');
    if numel(axisOut) < 16
        axisOut = linspace(max(opts.defaultWavelengthRange(1), min(axisOut)), ...
            min(opts.defaultWavelengthRange(2), max(axisOut)), 512)';
    end
end

function spectrumOut = interpolateSpectrumToAxis(spectrum, wavelengthNm)
    wl = double(spectrum(:, 1));
    y = double(spectrum(:, 2));
    valid = isfinite(wl) & isfinite(y);
    wl = wl(valid);
    y = y(valid);
    if isempty(wl)
        spectrumOut = zeros(size(wavelengthNm));
        return;
    end
    [wl, order] = sort(wl);
    y = y(order);
    [wlUnique, ~, groupIdx] = unique(wl);
    if numel(wlUnique) ~= numel(wl)
        y = accumarray(groupIdx, y, [], @mean);
        wl = wlUnique;
    end
    if isscalar(wl)
        spectrumOut = zeros(size(wavelengthNm));
        [~, nearestIdx] = min(abs(wavelengthNm - wl));
        spectrumOut(nearestIdx) = max(y, 0);
    else
        spectrumOut = interp1(wl, y, wavelengthNm, 'linear', 0);
        spectrumOut = max(spectrumOut, 0);
    end
end

function cubeOut = interpSpectralCube(cubeIn, wlIn, wlOut)
    wlIn = double(wlIn(:));
    wlOut = double(wlOut(:));
    [wlIn, uniqueIdx] = unique(wlIn, 'stable');
    cubeIn = cubeIn(:, :, uniqueIdx);
    [rows, cols, nW] = size(cubeIn);
    matIn = reshape(double(cubeIn), [], nW);
    matOutT = interp1(wlIn, matIn.', wlOut, 'linear', 0);
    cubeOut = reshape(matOutT.', rows, cols, numel(wlOut));
end

function [procMat, normMat, rows, cols, nW] = preprocessCubeToMatrix(cube, opts)
    [rows, cols, nW] = size(cube);
    procMat = double(reshape(cube, [], nW));
    procMat(~isfinite(procMat)) = 0;
    procMat = procMat - min(procMat, [], 2);
    procMat(procMat < 0) = 0;
    if opts.smoothingWindow > 1
        procMat = movmean(procMat, opts.smoothingWindow, 2);
    end
    normMat = normalizeRowsByMax(procMat);
end

function normMat = normalizeRowsByMax(procMat)
    maxValue = max(procMat, [], 2) + eps;
    normMat = procMat ./ maxValue;
end

function validPixelMask = buildValidPixelMask(expProcMat, simProcMat, rows, cols, opts)
    validPixelMask = max(expProcMat, [], 2) > opts.minPixelSignal & ...
        max(simProcMat, [], 2) > opts.minPixelSignal;
    if ~isempty(opts.roi)
        roi = round(double(opts.roi(:)'));
        if numel(roi) ~= 4
            error('precies:pixelwise_inversion:InvalidRoi', 'roi must be [rowStart rowEnd colStart colEnd].');
        end
        roi(1) = max(1, roi(1));
        roi(2) = min(rows, roi(2));
        roi(3) = max(1, roi(3));
        roi(4) = min(cols, roi(4));
        roiMask = false(rows, cols);
        roiMask(roi(1):roi(2), roi(3):roi(4)) = true;
        validPixelMask = validPixelMask & roiMask(:);
    end
end

function [inversionPixelMask, selectionInfo] = selectPixelsForInversion(initialMetrics, validPixelMask, opts)
    validPixelMask = logical(validPixelMask(:));
    if ~opts.autoSelectBadPixels
        inversionPixelMask = validPixelMask;
        selectionInfo = buildSelectionInfo(validPixelMask, inversionPixelMask, "all valid pixels");
        return;
    end

    ncc = initialMetrics.nccMap(:);
    nrmse = initialMetrics.nrmseMap(:);
    peakShift = max(abs(initialMetrics.deltaLambdaMQWMap(:)), abs(initialMetrics.deltaLambdaGaNMap(:)));
    fwhmShift = max(abs(initialMetrics.deltaFwhmMQWMap(:)), abs(initialMetrics.deltaFwhmGaNMap(:)));
    ratioError = abs(initialMetrics.ratioErrorMap(:));

    nccBad = validPixelMask & isfinite(ncc) & ncc < opts.badPixelNccThreshold;
    nrmseBad = validPixelMask & isfinite(nrmse) & nrmse > opts.badPixelNRMSEThreshold;
    peakBad = validPixelMask & isfinite(peakShift) & peakShift > opts.badPixelPeakShiftThresholdNm;
    fwhmBad = validPixelMask & isfinite(fwhmShift) & fwhmShift > opts.badPixelFwhmShiftThresholdNm;
    ratioBad = validPixelMask & isfinite(ratioError) & ratioError > opts.badPixelRatioErrorThresholdPercent;
    inversionPixelMask = nccBad | nrmseBad | peakBad | fwhmBad | ratioBad;

    validCount = nnz(validPixelMask);
    if validCount > 0
        maxSelected = max(opts.badPixelMinCount, ceil(opts.badPixelMaxFraction * validCount));
        if nnz(inversionPixelMask) > maxSelected
            score = buildBadPixelScore(ncc, nrmse, peakShift, fwhmShift, ratioError, validPixelMask, opts);
            validIdx = find(validPixelMask);
            [~, order] = sort(score(validIdx), 'descend', 'MissingPlacement', 'last');
            keepIdx = validIdx(order(1:maxSelected));
            cappedMask = false(size(validPixelMask));
            cappedMask(keepIdx) = true;
            inversionPixelMask = cappedMask;
        elseif nnz(inversionPixelMask) < opts.badPixelMinCount
            score = buildBadPixelScore(ncc, nrmse, peakShift, fwhmShift, ratioError, validPixelMask, opts);
            validIdx = find(validPixelMask);
            [~, order] = sort(score(validIdx), 'descend', 'MissingPlacement', 'last');
            keepCount = min(opts.badPixelMinCount, numel(validIdx));
            fallbackMask = false(size(validPixelMask));
            fallbackMask(validIdx(order(1:keepCount))) = true;
            inversionPixelMask = fallbackMask;
        end
    end

    selectionInfo = buildSelectionInfo(validPixelMask, inversionPixelMask, "automatic poor-fit pixels");
    selectionInfo.thresholds = struct( ...
        'NCC', opts.badPixelNccThreshold, ...
        'nRMSE', opts.badPixelNRMSEThreshold, ...
        'PeakShiftNm', opts.badPixelPeakShiftThresholdNm, ...
        'FwhmShiftNm', opts.badPixelFwhmShiftThresholdNm, ...
        'RatioErrorPercent', opts.badPixelRatioErrorThresholdPercent, ...
        'MaxFraction', opts.badPixelMaxFraction, ...
        'MinCount', opts.badPixelMinCount);
end

function score = buildBadPixelScore(ncc, nrmse, peakShift, fwhmShift, ratioError, validPixelMask, opts)
    score = -inf(size(validPixelMask));
    score(validPixelMask) = 0;
    score = score + max(0, (opts.badPixelNccThreshold - finiteVector(ncc)) ./ max(opts.badPixelNccThreshold, eps));
    score = score + max(0, finiteVector(nrmse) ./ max(opts.badPixelNRMSEThreshold, eps));
    score = score + max(0, finiteVector(peakShift) ./ max(opts.badPixelPeakShiftThresholdNm, eps));
    score = score + 0.5 * max(0, finiteVector(fwhmShift) ./ max(opts.badPixelFwhmShiftThresholdNm, eps));
    score = score + 0.5 * max(0, finiteVector(ratioError) ./ max(opts.badPixelRatioErrorThresholdPercent, eps));
    score(~validPixelMask) = -inf;
end

function values = finiteVector(values)
    values = double(values(:));
    values(~isfinite(values)) = 0;
end

function selectionInfo = buildSelectionInfo(validPixelMask, inversionPixelMask, modeText)
    validCount = nnz(validPixelMask);
    selectedCount = nnz(inversionPixelMask);
    selectionInfo = struct();
    selectionInfo.mode = char(modeText);
    selectionInfo.validPixelCount = validCount;
    selectionInfo.selectedPixelCount = selectedCount;
    selectionInfo.selectedFraction = selectedCount / max(validCount, 1);
end

function bounds = parameterBounds(opts)
    b = opts.parameterBounds;
    bounds = [
        b.ganGain
        b.mqwGain
        b.backgroundGain
        b.ganShiftNm
        b.mqwShiftNm
        b.ganWidthScale
        b.mqwWidthScale];
end

function p0 = initialPixelParameters(wl, targetProc, baselineProc, opts, bounds)
    targetNorm = normalizeSpectrum(targetProc);
    baselineNorm = normalizeSpectrum(baselineProc);
    ganGain = boundedRatio(windowArea(wl, targetProc, opts.ganWindow), ...
        windowArea(wl, baselineProc, opts.ganWindow), bounds(1, :));
    mqwGain = boundedRatio(windowArea(wl, targetProc, opts.mqwWindow), ...
        windowArea(wl, baselineProc, opts.mqwWindow), bounds(2, :));

    [lambdaTargetGaN, fwhmTargetGaN] = getPeakAndFWHM(wl, targetNorm, opts.ganWindow, opts);
    [lambdaBaseGaN, fwhmBaseGaN] = getPeakAndFWHM(wl, baselineNorm, opts.ganWindow, opts);
    [lambdaTargetMQW, fwhmTargetMQW] = getPeakAndFWHM(wl, targetNorm, opts.mqwWindow, opts);
    [lambdaBaseMQW, fwhmBaseMQW] = getPeakAndFWHM(wl, baselineNorm, opts.mqwWindow, opts);

    ganShift = finiteOrDefault(lambdaTargetGaN - lambdaBaseGaN, 0);
    mqwShift = finiteOrDefault(lambdaTargetMQW - lambdaBaseMQW, 0);
    ganWidthScale = finiteOrDefault(fwhmTargetGaN / max(fwhmBaseGaN, eps), 1);
    mqwWidthScale = finiteOrDefault(fwhmTargetMQW / max(fwhmBaseMQW, eps), 1);

    p0 = [
        ganGain
        mqwGain
        1
        ganShift
        mqwShift
        ganWidthScale
        mqwWidthScale];
    p0 = min(max(p0, bounds(:, 1) + 1e-6), bounds(:, 2) - 1e-6);
end

function ratio = boundedRatio(numerator, denominator, bounds)
    if denominator <= eps || numerator <= 0
        ratio = 1;
    else
        ratio = numerator / denominator;
    end
    ratio = min(max(ratio, bounds(1)), bounds(2));
end

function value = finiteOrDefault(value, defaultValue)
    if ~isfinite(value)
        value = defaultValue;
    end
end

function areaValue = windowArea(wl, spectrum, window)
    mask = wl >= window(1) & wl <= window(2);
    if nnz(mask) < 2
        areaValue = 0;
        return;
    end
    areaValue = trapz(wl(mask), max(double(spectrum(mask)), 0));
end

function u = paramsToUnconstrained(params, bounds)
    fraction = (params(:) - bounds(:, 1)) ./ max(bounds(:, 2) - bounds(:, 1), eps);
    fraction = min(max(fraction, 1e-5), 1 - 1e-5);
    u = log(fraction ./ (1 - fraction));
end

function params = unconstrainedToParams(u, bounds)
    fraction = 1 ./ (1 + exp(-double(u(:))));
    params = bounds(:, 1) + fraction .* (bounds(:, 2) - bounds(:, 1));
end

function objectiveValue = pixelObjective(params, p0, wl, targetNorm, targetProc, baselineProc, opts)
    fittedProc = synthesizeCorrectedSpectrum(wl, baselineProc, params, opts);
    fittedNorm = normalizeSpectrum(fittedProc);
    weights = 0.20 + 1.20 * sqrt(max(targetNorm(:), 0));
    weights = weights / max(mean(weights), eps);
    diff = targetNorm(:) - fittedNorm(:);
    rmseValue = sqrt(mean(weights .* diff .^ 2, 'omitnan'));
    nccValue = calculateNcc(targetNorm, fittedNorm);
    objectiveValue = rmseValue + opts.objectiveNccWeight * max(0, 1 - nccValue);

    [lambdaTargetMQW, fwhmTargetMQW] = getPeakAndFWHM(wl, targetNorm, opts.mqwWindow, opts);
    [lambdaFitMQW, fwhmFitMQW] = getPeakAndFWHM(wl, fittedNorm, opts.mqwWindow, opts);
    [lambdaTargetGaN, fwhmTargetGaN] = getPeakAndFWHM(wl, targetNorm, opts.ganWindow, opts);
    [lambdaFitGaN, fwhmFitGaN] = getPeakAndFWHM(wl, fittedNorm, opts.ganWindow, opts);

    objectiveValue = objectiveValue + peakPenalty(lambdaFitMQW, lambdaTargetMQW, opts) + ...
        peakPenalty(lambdaFitGaN, lambdaTargetGaN, opts) + ...
        fwhmPenalty(fwhmFitMQW, fwhmTargetMQW, opts) + ...
        fwhmPenalty(fwhmFitGaN, fwhmTargetGaN, opts);

    ratioTarget = calculateIntensityRatio(wl, targetProc, opts.ganWindow, opts.mqwWindow);
    ratioFit = calculateIntensityRatio(wl, fittedProc, opts.ganWindow, opts.mqwWindow);
    if isfinite(ratioTarget) && isfinite(ratioFit)
        ratioErrPercent = (ratioFit - ratioTarget) / (ratioTarget + eps) * 100;
        objectiveValue = objectiveValue + opts.objectiveRatioWeight * ...
            (ratioErrPercent / max(opts.ratioTolerancePercent, eps)) ^ 2;
    end

    scale = [1.5; 1.5; 0.8; 10; 10; 0.35; 0.35];
    objectiveValue = objectiveValue + opts.regularizationWeight * ...
        mean(((params(:) - p0(:)) ./ scale) .^ 2);
end

function penaltyValue = peakPenalty(lambdaFit, lambdaTarget, opts)
    penaltyValue = 0;
    if isfinite(lambdaFit) && isfinite(lambdaTarget)
        penaltyValue = opts.objectivePeakWeight * ...
            ((lambdaFit - lambdaTarget) / max(opts.peakToleranceNm, eps)) ^ 2;
    end
end

function penaltyValue = fwhmPenalty(fwhmFit, fwhmTarget, opts)
    penaltyValue = 0;
    if isfinite(fwhmFit) && isfinite(fwhmTarget)
        penaltyValue = opts.objectiveFwhmWeight * ...
            ((fwhmFit - fwhmTarget) / max(opts.fwhmToleranceNm, eps)) ^ 2;
    end
end

function fittedProc = synthesizeCorrectedSpectrum(wl, baselineProc, params, opts)
    baselineProc = max(double(baselineProc(:)), 0);
    ganWindow = buildSoftWindow(wl, opts.ganWindow, opts.softWindowEdgeNm);
    mqwWindow = buildSoftWindow(wl, opts.mqwWindow, opts.softWindowEdgeNm);
    occupiedWindow = min(1, ganWindow + mqwWindow);
    backgroundWindow = max(0, 1 - occupiedWindow);

    ganComponent = baselineProc .* ganWindow;
    mqwComponent = baselineProc .* mqwWindow;
    backgroundComponent = baselineProc .* backgroundWindow;

    ganCenter = componentCenterNm(wl, ganComponent, mean(opts.ganWindow));
    mqwCenter = componentCenterNm(wl, mqwComponent, mean(opts.mqwWindow));

    ganComponent = shiftAndScaleComponent(wl, ganComponent, params(4), params(6), ganCenter);
    mqwComponent = shiftAndScaleComponent(wl, mqwComponent, params(5), params(7), mqwCenter);

    fittedProc = params(1) * ganComponent + params(2) * mqwComponent + params(3) * backgroundComponent;
    fittedProc = max(fittedProc, 0);
end

function weight = buildSoftWindow(wl, window, edgeWidthNm)
    wl = double(wl(:));
    edgeWidthNm = max(edgeWidthNm, eps);
    rise = 1 ./ (1 + exp(-(wl - window(1)) / edgeWidthNm));
    fall = 1 ./ (1 + exp((wl - window(2)) / edgeWidthNm));
    weight = rise .* fall;
    if max(weight) > 0
        weight = weight / max(weight);
    end
end

function centerNm = componentCenterNm(wl, component, defaultCenterNm)
    component = max(double(component(:)), 0);
    if sum(component) <= eps
        centerNm = defaultCenterNm;
    else
        centerNm = sum(wl(:) .* component) / sum(component);
    end
end

function componentOut = shiftAndScaleComponent(wl, component, shiftNm, widthScale, centerNm)
    widthScale = max(widthScale, 0.05);
    queryAxis = centerNm + (wl(:) - centerNm - shiftNm) ./ widthScale;
    componentOut = interp1(wl(:), double(component(:)), queryAxis, 'linear', 0);
    componentOut = max(componentOut, 0);
end

function result = calculateMetricsFromMatrices(wl, expProcMat, expNormMat, simProcMat, simNormMat, rows, cols, opts, validPixelMask)
    nPix = rows * cols;
    nccVec = nan(nPix, 1);
    nrmseVec = nan(nPix, 1);
    dLambdaMQW = nan(nPix, 1);
    dLambdaGaN = nan(nPix, 1);
    dFwhmMQW = nan(nPix, 1);
    dFwhmGaN = nan(nPix, 1);
    ratioErr = nan(nPix, 1);
    residualRms = nan(nPix, 1);
    residualArea = nan(nPix, 1);

    diffMat = expNormMat - simNormMat;
    validIdx = find(validPixelMask(:)');
    for p = validIdx
        expSpecNorm = expNormMat(p, :)';
        simSpecNorm = simNormMat(p, :)';
        expSpecProc = expProcMat(p, :)';
        simSpecProc = simProcMat(p, :)';

        nccVec(p) = calculateNcc(expSpecNorm, simSpecNorm);
        rmseValue = sqrt(mean((expSpecNorm - simSpecNorm) .^ 2, 'omitnan'));
        expRange = max(expSpecNorm) - min(expSpecNorm) + eps;
        nrmseVec(p) = rmseValue / expRange;
        residualRms(p) = rmseValue;
        residualArea(p) = trapz(wl, expSpecNorm - simSpecNorm);

        [lambdaExpMQW, fwhmExpMQW] = getPeakAndFWHM(wl, expSpecNorm, opts.mqwWindow, opts);
        [lambdaSimMQW, fwhmSimMQW] = getPeakAndFWHM(wl, simSpecNorm, opts.mqwWindow, opts);
        [lambdaExpGaN, fwhmExpGaN] = getPeakAndFWHM(wl, expSpecNorm, opts.ganWindow, opts);
        [lambdaSimGaN, fwhmSimGaN] = getPeakAndFWHM(wl, simSpecNorm, opts.ganWindow, opts);

        dLambdaMQW(p) = lambdaSimMQW - lambdaExpMQW;
        dLambdaGaN(p) = lambdaSimGaN - lambdaExpGaN;
        dFwhmMQW(p) = fwhmSimMQW - fwhmExpMQW;
        dFwhmGaN(p) = fwhmSimGaN - fwhmExpGaN;

        ratioExp = calculateIntensityRatio(wl, expSpecProc, opts.ganWindow, opts.mqwWindow);
        ratioSim = calculateIntensityRatio(wl, simSpecProc, opts.ganWindow, opts.mqwWindow);
        ratioErr(p) = (ratioSim - ratioExp) / (ratioExp + eps) * 100;
    end

    meanResidual = mean(diffMat(validPixelMask, :), 1, 'omitnan')';
    result = struct();
    result.wavelength = wl(:);
    result.expNormCube = reshape(single(expNormMat), rows, cols, []);
    result.simNormCube = reshape(single(simNormMat), rows, cols, []);
    result.expProcCube = reshape(single(expProcMat), rows, cols, []);
    result.simProcCube = reshape(single(simProcMat), rows, cols, []);
    result.nccMap = reshape(nccVec, rows, cols);
    result.nrmseMap = reshape(nrmseVec, rows, cols);
    result.deltaLambdaMQWMap = reshape(dLambdaMQW, rows, cols);
    result.deltaLambdaGaNMap = reshape(dLambdaGaN, rows, cols);
    result.deltaFwhmMQWMap = reshape(dFwhmMQW, rows, cols);
    result.deltaFwhmGaNMap = reshape(dFwhmGaN, rows, cols);
    result.ratioErrorMap = reshape(ratioErr, rows, cols);
    result.residualRmsMap = reshape(residualRms, rows, cols);
    result.signedResidualAreaMap = reshape(residualArea, rows, cols);
    result.meanResidualSpectrum = meanResidual;
    result.summary = summarizeMetrics(nccVec, nrmseVec, dLambdaMQW, dLambdaGaN, dFwhmMQW, dFwhmGaN, ratioErr);
    if opts.storeResidualCube
        result.residualNormCube = reshape(single(diffMat), rows, cols, []);
    end
end

function nccValue = calculateNcc(referenceSpectrum, testSpectrum)
    referenceZero = referenceSpectrum(:) - mean(referenceSpectrum(:), 'omitnan');
    testZero = testSpectrum(:) - mean(testSpectrum(:), 'omitnan');
    numerator = sum(referenceZero .* testZero, 'omitnan');
    denominator = sqrt(sum(referenceZero .^ 2, 'omitnan') * sum(testZero .^ 2, 'omitnan')) + eps;
    nccValue = numerator / denominator;
    if ~isfinite(nccValue)
        nccValue = NaN;
    end
end

function summary = summarizeMetrics(nccVec, nrmseVec, dLambdaMQW, dLambdaGaN, dFwhmMQW, dFwhmGaN, ratioErr)
    summary = struct();
    summary.averageNCC = mean(nccVec, 'omitnan');
    summary.stdNCC = std(nccVec, 0, 'omitnan');
    summary.averageNRMSE = mean(nrmseVec, 'omitnan');
    summary.stdNRMSE = std(nrmseVec, 0, 'omitnan');
    summary.meanAbsDeltaLambdaMQW = mean(abs(dLambdaMQW), 'omitnan');
    summary.stdAbsDeltaLambdaMQW = std(abs(dLambdaMQW), 0, 'omitnan');
    summary.meanAbsDeltaLambdaGaN = mean(abs(dLambdaGaN), 'omitnan');
    summary.stdAbsDeltaLambdaGaN = std(abs(dLambdaGaN), 0, 'omitnan');
    summary.meanAbsDeltaFwhmMQW = mean(abs(dFwhmMQW), 'omitnan');
    summary.stdAbsDeltaFwhmMQW = std(abs(dFwhmMQW), 0, 'omitnan');
    summary.meanAbsDeltaFwhmGaN = mean(abs(dFwhmGaN), 'omitnan');
    summary.stdAbsDeltaFwhmGaN = std(abs(dFwhmGaN), 0, 'omitnan');
    summary.meanAbsRatioError = mean(abs(ratioErr), 'omitnan');
    summary.stdAbsRatioError = std(abs(ratioErr), 0, 'omitnan');
end

function T = buildSummaryTable(result)
    voltage = result.voltage;
    avgNCC = result.summary.averageNCC;
    stdNCC = result.summary.stdNCC;
    avgNRMSE = result.summary.averageNRMSE;
    stdNRMSE = result.summary.stdNRMSE;
    dLamMQW = result.summary.meanAbsDeltaLambdaMQW;
    sdLamMQW = result.summary.stdAbsDeltaLambdaMQW;
    dLamGaN = result.summary.meanAbsDeltaLambdaGaN;
    sdLamGaN = result.summary.stdAbsDeltaLambdaGaN;
    dFWHMMQW = result.summary.meanAbsDeltaFwhmMQW;
    sdFWHMMQW = result.summary.stdAbsDeltaFwhmMQW;
    dFWHMGaN = result.summary.meanAbsDeltaFwhmGaN;
    sdFWHMGaN = result.summary.stdAbsDeltaFwhmGaN;
    ratioErr = result.summary.meanAbsRatioError;
    sdRatioErr = result.summary.stdAbsRatioError;
    T = table(voltage, avgNCC, stdNCC, avgNRMSE, stdNRMSE, ...
        dLamMQW, sdLamMQW, dLamGaN, sdLamGaN, ...
        dFWHMMQW, sdFWHMMQW, dFWHMGaN, sdFWHMGaN, ratioErr, sdRatioErr, ...
        'VariableNames', {'Voltage_kV', 'Average_NCC', 'SD_NCC', ...
        'nRMSE', 'SD_nRMSE', ...
        'Abs_DeltaLambda_MQW_nm', 'SD_Abs_DeltaLambda_MQW_nm', ...
        'Abs_DeltaLambda_GaN_nm', 'SD_Abs_DeltaLambda_GaN_nm', ...
        'Abs_DeltaFWHM_MQW_nm', 'SD_Abs_DeltaFWHM_MQW_nm', ...
        'Abs_DeltaFWHM_GaN_nm', 'SD_Abs_DeltaFWHM_GaN_nm', ...
        'Abs_Error_IGaN_over_IMQW_percent', 'SD_Abs_Error_IGaN_over_IMQW_percent'});
end

function parameterMaps = buildParameterMaps(parameterValues, objectiveMap, exitFlagMap, iterationMap, rows, cols)
    parameterMaps = struct();
    parameterMaps.ganGainMap = reshape(parameterValues(:, 1), rows, cols);
    parameterMaps.mqwGainMap = reshape(parameterValues(:, 2), rows, cols);
    parameterMaps.backgroundGainMap = reshape(parameterValues(:, 3), rows, cols);
    parameterMaps.ganShiftNmMap = reshape(parameterValues(:, 4), rows, cols);
    parameterMaps.mqwShiftNmMap = reshape(parameterValues(:, 5), rows, cols);
    parameterMaps.ganWidthScaleMap = reshape(parameterValues(:, 6), rows, cols);
    parameterMaps.mqwWidthScaleMap = reshape(parameterValues(:, 7), rows, cols);
    parameterMaps.objectiveMap = reshape(objectiveMap, rows, cols);
    parameterMaps.exitFlagMap = reshape(exitFlagMap, rows, cols);
    parameterMaps.iterationMap = reshape(iterationMap, rows, cols);
end

function data = buildSpectraData(wavelengthNm, procMat, rows, cols)
    data = struct();
    data.wavelengthAxis = wavelengthNm(:);
    data.wavelengthIntensityMaps = reshape(single(procMat), rows, cols, []);
    data.totalSpectra = cell(rows, cols);
    for idx = 1:rows * cols
        [r, c] = ind2sub([rows, cols], idx);
        data.totalSpectra{r, c} = [wavelengthNm(:), procMat(idx, :)'];
    end
end

function writeInversionOutputs(result, opts)
    if ~exist(opts.outputDirectory, 'dir')
        mkdir(opts.outputDirectory);
    end
    writetable(result.summaryTable, fullfile(opts.outputDirectory, 'PixelwiseInversionMetrics.csv'));
    save(fullfile(opts.outputDirectory, 'PixelwiseInversionResult.mat'), 'result', '-v7.3');
end

function [lambdaPeak, fwhmValue] = getPeakAndFWHM(wl, spec, window, opts)
    idx = wl >= window(1) & wl <= window(2);
    x = wl(idx);
    y = double(spec(idx));
    if numel(x) < 5 || all(~isfinite(y))
        lambdaPeak = NaN;
        fwhmValue = NaN;
        return;
    end

    y(~isfinite(y)) = 0;
    y = y - min(y);
    if max(y) < opts.minPeakRelativeHeight
        lambdaPeak = NaN;
        fwhmValue = NaN;
        return;
    end

    [peakVal, peakIdx] = max(y);
    lambdaPeak = x(peakIdx);
    if peakIdx > 1 && peakIdx < numel(y)
        x3 = x(peakIdx-1:peakIdx+1);
        y3 = y(peakIdx-1:peakIdx+1);
        polyCoeff = polyfit(x3, y3, 2);
        if isfinite(polyCoeff(1)) && polyCoeff(1) < 0
            lambdaInterp = -polyCoeff(2) / (2 * polyCoeff(1));
            if lambdaInterp >= x3(1) && lambdaInterp <= x3(end)
                lambdaPeak = lambdaInterp;
            end
        end
    end

    halfMax = 0.5 * peakVal;
    leftIdx = find(y(1:peakIdx) <= halfMax, 1, 'last');
    rightRelIdx = find(y(peakIdx:end) <= halfMax, 1, 'first');
    if isempty(leftIdx) || isempty(rightRelIdx) || leftIdx == peakIdx
        fwhmValue = NaN;
        return;
    end
    rightIdx = peakIdx + rightRelIdx - 1;
    leftX = interpHalfMax(x(leftIdx:leftIdx+1), y(leftIdx:leftIdx+1), halfMax);
    rightX = interpHalfMax(x(rightIdx-1:rightIdx), y(rightIdx-1:rightIdx), halfMax);
    fwhmValue = rightX - leftX;
end

function xHalf = interpHalfMax(xPair, yPair, halfMax)
    if abs(diff(yPair)) <= eps
        xHalf = mean(xPair);
    else
        xHalf = interp1(yPair, xPair, halfMax, 'linear', 'extrap');
    end
end

function ratio = calculateIntensityRatio(wl, spec, ganWindow, mqwWindow)
    ganArea = windowArea(wl, spec, ganWindow);
    mqwArea = windowArea(wl, spec, mqwWindow);
    ratio = ganArea / (mqwArea + eps);
end

function normalizedSpectrum = normalizeSpectrum(values)
    normalizedSpectrum = double(values(:));
    normalizedSpectrum(~isfinite(normalizedSpectrum)) = 0;
    normalizedSpectrum = max(normalizedSpectrum, 0);
    peakValue = max(normalizedSpectrum);
    if peakValue > 0
        normalizedSpectrum = normalizedSpectrum / peakValue;
    end
end

function merged = mergeOptions(defaults, overrides)
    merged = defaults;
    if ~isstruct(overrides)
        return;
    end
    fields = fieldnames(overrides);
    for idx = 1:numel(fields)
        fieldName = fields{idx};
        if strcmp(fieldName, 'parameterBounds') && isstruct(overrides.parameterBounds)
            merged.parameterBounds = mergeOptions(merged.parameterBounds, overrides.parameterBounds);
        else
            merged.(fieldName) = overrides.(fieldName);
        end
    end
end

function value = getFieldOrDefault(data, fieldName, defaultValue)
    if isstruct(data) && isfield(data, fieldName) && ~isempty(data.(fieldName))
        value = data.(fieldName);
    else
        value = defaultValue;
    end
end

function ok = selfTest()
    rng(8);
    wl = linspace(300, 800, 256)';
    rows = 2;
    cols = 2;
    simCube = zeros(rows, cols, numel(wl));
    expCube = zeros(rows, cols, numel(wl));
    for pixelIdx = 1:rows * cols
        ganShift = -3 + pixelIdx;
        mqwShift = 3 - 0.8 * pixelIdx;
        gan = 0.18 * gaussianProfile(wl, 388, 9);
        mqw = gaussianProfile(wl, 515, 24);
        base = gan + mqw + 0.02 * gaussianProfile(wl, 455, 55);
        target = 1.18 * shiftProfile(wl, gan, ganShift) + ...
            0.92 * shiftProfile(wl, mqw, mqwShift) + 0.02 * gaussianProfile(wl, 455, 55);
        target = target + 0.002 * rand(size(target));
        [r, c] = ind2sub([rows, cols], pixelIdx);
        simCube(r, c, :) = base;
        expCube(r, c, :) = target;
    end

    expData = struct('wavelengthAxis', wl, 'wavelengthIntensityMaps', expCube);
    simData = struct('wavelengthAxis', wl, 'wavelengthIntensityMaps', simCube);
    opts = struct('maxIterations', 35, 'maxFunctionEvaluations', 160, 'smoothingWindow', 3);
    result = runInversion(expData, simData, opts);
    assert(isfield(result, 'nccMap') && all(isfinite(result.nccMap(:))), 'NCC map missing or invalid.');
    assert(isfield(result, 'inversionPixelMask') && any(result.inversionPixelMask(:)), ...
        'Automatic poor-pixel selection did not select any pixels in self-test.');
    assert(result.summary.averageNCC >= result.initialMetrics.summary.averageNCC - 1e-6, ...
        'Pixel-wise inversion did not improve average NCC in self-test.');
    assert(result.summary.averageNRMSE <= result.initialMetrics.summary.averageNRMSE + 1e-6, ...
        'Pixel-wise inversion did not improve average nRMSE in self-test.');
    ok = true;
    fprintf('precies.pixelwise_inversion selftest passed. NCC %.4f -> %.4f, nRMSE %.4f -> %.4f\n', ...
        result.initialMetrics.summary.averageNCC, result.summary.averageNCC, ...
        result.initialMetrics.summary.averageNRMSE, result.summary.averageNRMSE);
end

function y = gaussianProfile(wl, centerNm, fwhmNm)
    sigma = fwhmNm / (2 * sqrt(2 * log(2)));
    y = exp(-0.5 * ((wl - centerNm) / sigma) .^ 2);
end

function y = shiftProfile(wl, values, shiftNm)
    y = interp1(wl, values, wl - shiftNm, 'linear', 0);
end
