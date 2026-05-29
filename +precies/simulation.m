function varargout = simulation(action, varargin)
    switch lower(string(action))
        case "electrondepthsweep"
            [varargout{1:nargout}] = electronDepthSweep(varargin{:});
        case "simulateoptions"
            [varargout{1:nargout}] = simulateOptions(varargin{:});
        otherwise
            error('precies:simulation:InvalidAction', 'Unsupported action: %s', action);
    end
end
function [positions, intensities, wavelengths, sourceCount] = appendLightSourcesToBuffer( ...
    positions, intensities, wavelengths, sourceCount, clSources, clIntensities, clWavelengths, ...
    plSources, plIntensities, plWavelengths)

    numCl = size(clSources, 1);
    numPl = size(plSources, 1);
    numNew = numCl + numPl;
    if numNew < 1
        return;
    end

    [positions, intensities, wavelengths] = ensureLightBufferCapacity( ...
        positions, intensities, wavelengths, sourceCount + numNew, max(numNew, 1024));
    writeRange = (sourceCount + 1):(sourceCount + numNew);
    positions(writeRange, :) = [clSources; plSources];
    intensities(writeRange) = [clIntensities; plIntensities];
    wavelengths(writeRange) = [clWavelengths; plWavelengths];
    sourceCount = sourceCount + numNew;
end

function [positions, intensities, wavelengths] = ensureLightBufferCapacity(positions, intensities, wavelengths, requiredCount, growthHint)
    currentCapacity = size(positions, 1);
    if requiredCount <= currentCapacity
        return;
    end

    if nargin < 5 || isempty(growthHint)
        growthHint = max(1024, currentCapacity);
    end
    newCapacity = max([requiredCount, currentCapacity + growthHint, ceil(max(currentCapacity, 1) * 1.5)]);
    newPositions = zeros(newCapacity, 3, 'like', positions);
    newIntensities = zeros(newCapacity, 1, 'like', intensities);
    newWavelengths = zeros(newCapacity, 1, 'like', wavelengths);
    if currentCapacity > 0
        newPositions(1:currentCapacity, :) = positions;
        newIntensities(1:currentCapacity) = intensities;
        newWavelengths(1:currentCapacity) = wavelengths;
    end
    positions = newPositions;
    intensities = newIntensities;
    wavelengths = newWavelengths;
end

function matrixOut = growMatrixRows(matrixIn, requiredRows)
    currentRows = size(matrixIn, 1);
    if requiredRows <= currentRows
        matrixOut = matrixIn;
        return;
    end

    newRows = max([requiredRows, currentRows + max(1, currentRows), ceil(max(currentRows, 1) * 1.5)]);
    matrixOut = zeros(newRows, size(matrixIn, 2), 'like', matrixIn);
    if currentRows > 0
        matrixOut(1:currentRows, :) = matrixIn;
    end
end

function [wavelengthBuffer, intensityBuffer, sampleCount] = appendSpectrumSamplesToBuffer( ...
    wavelengthBuffer, intensityBuffer, sampleCount, newWavelengths, newIntensities)

    newWavelengths = double(newWavelengths(:));
    newIntensities = double(newIntensities(:));
    nNew = min(numel(newWavelengths), numel(newIntensities));
    if nNew < 1
        return;
    end

    requiredCount = sampleCount + nNew;
    if requiredCount > numel(wavelengthBuffer)
        newCapacity = max([requiredCount, numel(wavelengthBuffer) + max(nNew, 128), ceil(max(numel(wavelengthBuffer), 1) * 1.5)]);
        newWavelengthBuffer = zeros(newCapacity, 1, 'like', wavelengthBuffer);
        newIntensityBuffer = zeros(newCapacity, 1, 'like', intensityBuffer);
        if sampleCount > 0
            newWavelengthBuffer(1:sampleCount) = wavelengthBuffer(1:sampleCount);
            newIntensityBuffer(1:sampleCount) = intensityBuffer(1:sampleCount);
        end
        wavelengthBuffer = newWavelengthBuffer;
        intensityBuffer = newIntensityBuffer;
    end

    writeRange = sampleCount + 1:requiredCount;
    wavelengthBuffer(writeRange) = newWavelengths(1:nNew);
    intensityBuffer(writeRange) = newIntensities(1:nNew);
    sampleCount = requiredCount;
end

function simulationData = simulateOptions(options)
    opts = normalizeForwardSimulationOptions(options);
    if ~isempty(opts.rngSeed)
        rng(opts.rngSeed);
    end

    params = opts.params;
    params.modelToggles = mergeModelToggles(getDefaultModelToggles(), ...
        getFieldOrDefault(params, 'modelToggles', struct()));
    params.structuralParameters = mergeStructs( ...
        getFieldOrDefault(params, 'structuralParameters', struct()), opts.structuralParameters);
    params = applyStructuralParametersToSimulationParams(params);

    Eg_GaN = 3.3032;
    Eg_InN = 0.6086;
    entryConfig = precies.config(opts.entryMode, params);
    layer_params = entryConfig.layerParams;
    [layers, layerBoundaries] = Layers(layer_params, params);
    [beam_positions, pixelsX, pixelsY] = generatePixelGrid(params.scanX, params.scanY, params.totalPixels);
    params.pixelsX = pixelsX;
    params.pixelsY = pixelsY;

    total_thickness = (layerBoundaries(end - 1) - layerBoundaries(1)) * 1e9;
    depth_planes_nm = 0:params.depth_step:total_thickness;
    if numel(depth_planes_nm) < 2
        depth_planes_nm = [0, max(params.depth_step, total_thickness)];
    end
    depth_planes_bins = [depth_planes_nm(1:end-1); depth_planes_nm(2:end)]';

    Vpits = opts.Vpits;
    if isempty(Vpits) && opts.defaultVpitCount > 0
        Vpits = generateVpits(params.scanX, params.scanY, opts.defaultVpitDensity, opts.defaultVpitDepthNm);
        vpitGeometrySource = 'generated-default';
    elseif isempty(Vpits)
        vpitGeometrySource = 'none';
    else
        vpitGeometrySource = 'user-supplied';
    end
    Vpits = applyStructuralParametersToVpits(Vpits, params.structuralParameters);
    Vpits = precies.vpit('buildGeometries', Vpits, layers, layerBoundaries);
    nominalSurfaceZ_nm = layerBoundaries(end - 1) * 1e9;

    totalIntensityMap = zeros(pixelsY, pixelsX);
    depthIntensityMap = zeros(pixelsY, pixelsX, size(depth_planes_bins, 1));
    layerIntensityMap = zeros(pixelsY, pixelsX, size(layers, 1));
    totalSpectra = cell(pixelsY, pixelsX);
    depthSpectra = cell(pixelsY, pixelsX, size(depth_planes_bins, 1));
    layerSpectra = cell(pixelsY, pixelsX, size(layers, 1));
    sourceCapacity = max(1024, size(beam_positions, 1) * 256);
    allLightPositions = zeros(sourceCapacity, 3);
    allLightIntensities = zeros(sourceCapacity, 1);
    allLightWavelengths = zeros(sourceCapacity, 1);
    sourceCount = 0;

    for idx = 1:size(beam_positions, 1)
        CL_sources = run_electron_simulation(layer_params, params, Eg_GaN, Eg_InN, beam_positions(idx, :), Vpits);
        if isempty(CL_sources)
            continue;
        end
        result = mainMonteCarloSim(params, layer_params, Eg_GaN, Eg_InN, ...
            CL_sources, depth_planes_bins, idx, Vpits, nominalSurfaceZ_nm);

        [rowIdx, colIdx] = ind2sub([pixelsY, pixelsX], idx);
        totalIntensityMap(rowIdx, colIdx) = result.totalIntensity;
        totalSpectra{rowIdx, colIdx} = result.totalSpectrum;
        depthIntensityMap(rowIdx, colIdx, :) = reshape(result.depthIntensity, 1, 1, []);
        layerIntensityMap(rowIdx, colIdx, :) = reshape(result.layerIntensity, 1, 1, []);
        for binIdx = 1:size(depth_planes_bins, 1)
            depthSpectra{rowIdx, colIdx, binIdx} = result.depthSpectra{binIdx};
        end
        for layerIdx = 1:size(layers, 1)
            layerSpectra{rowIdx, colIdx, layerIdx} = result.layerSpectra{layerIdx};
        end

        [allLightPositions, allLightIntensities, allLightWavelengths, sourceCount] = appendLightSourcesToBuffer( ...
            allLightPositions, allLightIntensities, allLightWavelengths, sourceCount, ...
            result.CLSources, result.CLIntensities, result.CLWavelengths, ...
            result.PLSources, result.PLIntensities, result.PLWavelengths);

        if opts.progressEveryPixels > 0 && (mod(idx, opts.progressEveryPixels) == 0 || idx == size(beam_positions, 1))
            fprintf('simulateOptions: %d/%d pixels\n', idx, size(beam_positions, 1));
        end
    end
    allLightPositions = allLightPositions(1:sourceCount, :);
    allLightIntensities = allLightIntensities(1:sourceCount);
    allLightWavelengths = allLightWavelengths(1:sourceCount);

    layerCenters = (layerBoundaries(1:end-1) + layerBoundaries(2:end)) / 2;
    layerThicknesses = diff(layerBoundaries);
    simulationData = struct( ...
        'totalIntensityMap', totalIntensityMap, ...
        'depthIntensityMap', depthIntensityMap, ...
        'layerIntensityMap', layerIntensityMap, ...
        'totalSpectra', {totalSpectra}, ...
        'depthSpectra', {depthSpectra}, ...
        'layerSpectra', {layerSpectra}, ...
        'depthBins', depth_planes_bins, ...
        'layerNames', {layers(:, 1)'}, ...
        'positions', allLightPositions, ...
        'intensities', allLightIntensities, ...
        'wavelengths', allLightWavelengths, ...
        'Vpits', {Vpits}, ...
        'vpitGeometrySource', vpitGeometrySource, ...
        'vpitCount', numel(Vpits), ...
        'nominalSurfaceZ_nm', nominalSurfaceZ_nm, ...
        'isImportedTiff', false, ...
        'spectralRows', pixelsY, ...
        'spectralCols', pixelsX, ...
        'layerBoundaries', layerBoundaries, ...
        'layerCenters', layerCenters, ...
        'layerThicknesses', layerThicknesses, ...
        'simulationParams', params, ...
        'modelToggles', params.modelToggles, ...
        'structuralParameters', params.structuralParameters);
end

function opts = normalizeForwardSimulationOptions(options)
    if nargin < 1 || isempty(options)
        options = struct();
    end
    if ~isstruct(options)
        error('precies:simulation:InvalidOptions', 'simulateOptions expects a struct of options.');
    end

    defaultParams = struct( ...
        'InteTime', 0.3, ...
        'numRays', 64, ...
        'electronEnergy', 5e3, ...
        'beamCurrent', 100e-12, ...
        'depth_step', 25, ...
        'mqw_pairs', 10, ...
        'mqw_barriers', 9, ...
        'scanX', 120, ...
        'scanY', 120, ...
        'totalPixels', 4, ...
        'gridPoints3D', 24, ...
        'barrierThick', 8, ...
        'wellThick', 5, ...
        'prestrainedPeriods', 10, ...
        'prestrainedThickLayerNm', 9, ...
        'prestrainedThinLayerNm', 2, ...
        'prestrainedTotalThicknessNm', 120, ...
        'elementaryCharge', 1.602e-19, ...
        'numElectrons', 16);

    opts = struct();
    opts.entryMode = char(string(getOptionOrDefault(options, {'entryMode', 'mode'}, 'wavelength')));
    opts.params = mergeStructs(defaultParams, getOptionOrDefault(options, {'params'}, struct()));
    directParamFields = fieldnames(defaultParams);
    for fieldIdx = 1:numel(directParamFields)
        fieldName = directParamFields{fieldIdx};
        if isfield(options, fieldName) && ~isempty(options.(fieldName))
            opts.params.(fieldName) = options.(fieldName);
        end
    end
    if isfield(options, 'voltageKeV') && ~isempty(options.voltageKeV)
        opts.params.electronEnergy = double(options.voltageKeV) * 1e3;
    elseif isfield(options, 'voltage_kV') && ~isempty(options.voltage_kV)
        opts.params.electronEnergy = double(options.voltage_kV) * 1e3;
    end
    opts.params.modelToggles = mergeModelToggles(getDefaultModelToggles(), ...
        getOptionOrDefault(options, {'modelToggles'}, getFieldOrDefault(opts.params, 'modelToggles', struct())));
    opts.structuralParameters = getOptionOrDefault(options, {'structuralParameters', 'structuralParams'}, ...
        getFieldOrDefault(opts.params, 'structuralParameters', struct()));

    opts.Vpits = getOptionOrDefault(options, {'Vpits', 'vpits'}, cell(0, 1));
    opts.defaultVpitCount = max(0, round(double(getOptionOrDefault(options, {'defaultVpitCount'}, 0))));
    opts.defaultVpitDensity = double(getOptionOrDefault(options, {'defaultVpitDensity'}, 10));
    opts.defaultVpitDepthNm = double(getOptionOrDefault(options, {'defaultVpitDepthNm'}, 120));
    opts.rngSeed = getOptionOrDefault(options, {'rngSeed', 'seed'}, []);
    opts.progressEveryPixels = max(0, round(double(getOptionOrDefault(options, {'progressEveryPixels'}, 0))));

    numericFields = {'numRays', 'electronEnergy', 'beamCurrent', 'depth_step', 'mqw_pairs', 'mqw_barriers', ...
        'scanX', 'scanY', 'totalPixels', 'barrierThick', 'wellThick', ...
        'prestrainedPeriods', 'prestrainedThickLayerNm', 'prestrainedThinLayerNm', ...
        'prestrainedTotalThicknessNm', 'numElectrons'};
    for fieldIdx = 1:numel(numericFields)
        fieldName = numericFields{fieldIdx};
        opts.params.(fieldName) = double(opts.params.(fieldName));
    end
    opts.params.numRays = max(4, round(opts.params.numRays));
    opts.params.mqw_pairs = max(1, round(opts.params.mqw_pairs));
    opts.params.mqw_barriers = max(0, round(opts.params.mqw_barriers));
    opts.params.prestrainedPeriods = max(1, round(opts.params.prestrainedPeriods));
    opts.params.prestrainedThickLayerNm = max(0.1, opts.params.prestrainedThickLayerNm);
    opts.params.prestrainedThinLayerNm = max(0.1, opts.params.prestrainedThinLayerNm);
    opts.params.prestrainedTotalThicknessNm = max(0.2, opts.params.prestrainedTotalThicknessNm);
    opts.params.totalPixels = max(1, round(opts.params.totalPixels));
    opts.params.numElectrons = max(1, round(opts.params.numElectrons));
    opts.params.depth_step = max(1, opts.params.depth_step);
end

function errorMsg = buildDetailedErrorMessage(ME)
    errorMsg = sprintf('Error Message:\n%s', ME.message);

    try
        if ~isempty(ME.cause)
            errorMsg = sprintf('%s\n\nCaused By:', errorMsg);
            for causeIdx = 1:numel(ME.cause)
                causeReport = getReport(ME.cause{causeIdx}, 'extended', 'hyperlinks', 'off');
                errorMsg = sprintf('%s\n[%d]\n%s', errorMsg, causeIdx, causeReport);
            end
        end

        if ~isempty(ME.stack)
            errorMsg = sprintf('%s\n\nCall Stack:', errorMsg);
            for stackIdx = 1:length(ME.stack)
                stack = ME.stack(stackIdx);
                if isfield(stack, 'file')
                    errorMsg = sprintf('%s\n[%d] %s (Line %d)', ...
                        errorMsg, stackIdx, stack.file, stack.line);
                else
                    errorMsg = sprintf('%s\n[%d] [Internal Function] %s', ...
                        errorMsg, stackIdx, stack.name);
                end
            end
        end
    catch
    end
end

function resultsTable = electronDepthSweep(options)
    opts = normalizeElectronDepthOptions(options);

    params = struct(...
        'mqw_pairs', opts.mqwPairs, ...
        'mqw_barriers', opts.mqwBarriers, ...
        'barrierThick', opts.barrierThick, ...
        'wellThick', opts.wellThick, ...
        'prestrainedPeriods', opts.prestrainedPeriods, ...
        'prestrainedThickLayerNm', opts.prestrainedThickLayerNm, ...
        'prestrainedThinLayerNm', opts.prestrainedThinLayerNm, ...
        'prestrainedTotalThicknessNm', opts.prestrainedTotalThicknessNm);
    entryConfig = precies.config(opts.entryMode, params);
    layer_params = entryConfig.layerParams;

    Eg_GaN = 3.3032;
    Eg_InN = 0.6086;
    start_pos_0 = reshape(opts.startPositionNm, 1, []);
    if numel(start_pos_0) ~= 3
        error('precies:simulation:InvalidStartPosition', 'startPositionNm must contain exactly three values: [x, y, z].');
    end
    start_pos_0 = start_pos_0 * 1e-9;
    Vpits = opts.Vpits;

    numVoltages = numel(opts.voltagesKeV);
    resultsTable = table('Size', [numVoltages, 7], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'Voltage_keV', 'MeanDepth_nm', 'MedianDepth_nm', 'P95Depth_nm', 'MaxDepth_nm', 'StdDepth_nm', 'NumSamples'});

    if ~isempty(opts.rngSeed)
        rng(opts.rngSeed);
    end

    for voltageIdx = 1:numVoltages
        params.electronEnergy = opts.voltagesKeV(voltageIdx) * 1e3;
        penetrationDepthCells = cell(opts.numRepeats, 1);

        for repeatIdx = 1:opts.numRepeats
            [~, electronStats] = simulateElectronExcitations( ...
                layer_params, params, Eg_GaN, Eg_InN, start_pos_0, Vpits, opts.numElectrons);
            penetrationDepthCells{repeatIdx} = electronStats.penetrationDepths_nm(:);
        end
        allPenetrationDepths_nm = vertcat(penetrationDepthCells{:});

        resultsTable.Voltage_keV(voltageIdx) = opts.voltagesKeV(voltageIdx);
        resultsTable.NumSamples(voltageIdx) = numel(allPenetrationDepths_nm);

        if isempty(allPenetrationDepths_nm)
            resultsTable.MeanDepth_nm(voltageIdx) = NaN;
            resultsTable.MedianDepth_nm(voltageIdx) = NaN;
            resultsTable.P95Depth_nm(voltageIdx) = NaN;
            resultsTable.MaxDepth_nm(voltageIdx) = NaN;
            resultsTable.StdDepth_nm(voltageIdx) = NaN;
        else
            resultsTable.MeanDepth_nm(voltageIdx) = mean(allPenetrationDepths_nm, 'omitnan');
            resultsTable.MedianDepth_nm(voltageIdx) = median(allPenetrationDepths_nm, 'omitnan');
            resultsTable.P95Depth_nm(voltageIdx) = prctile(allPenetrationDepths_nm, 95);
            resultsTable.MaxDepth_nm(voltageIdx) = max(allPenetrationDepths_nm);
            resultsTable.StdDepth_nm(voltageIdx) = std(allPenetrationDepths_nm, 0, 'omitnan');
        end
    end
end

function opts = normalizeElectronDepthOptions(options)
    if nargin < 1 || isempty(options)
        options = struct();
    end
    if ~isstruct(options)
        error('precies:simulation:InvalidOptions', 'electronDepthSweep expects a struct of options.');
    end

    opts = struct();
    opts.entryMode = getOptionOrDefault(options, {'entryMode', 'mode'}, 'wavelength');
    opts.voltagesKeV = reshape(getOptionOrDefault(options, {'voltagesKeV', 'voltages'}, 1:10), 1, []);
    opts.numElectrons = getOptionOrDefault(options, {'numElectrons', 'num_electrons'}, 100);
    opts.numRepeats = getOptionOrDefault(options, {'numRepeats', 'num_repeats'}, 3);
    opts.mqwPairs = getOptionOrDefault(options, {'mqwPairs', 'mqw_pairs'}, 10);
    opts.mqwBarriers = getOptionOrDefault(options, {'mqwBarriers', 'mqw_barriers', 'barrierCount'}, opts.mqwPairs - 1);
    opts.barrierThick = getOptionOrDefault(options, {'barrierThick', 'barrier_nm'}, 8);
    opts.wellThick = getOptionOrDefault(options, {'wellThick', 'well_nm'}, 5);
    opts.prestrainedPeriods = getOptionOrDefault(options, {'prestrainedPeriods', 'prelayerPeriods'}, 10);
    opts.prestrainedThickLayerNm = getOptionOrDefault(options, {'prestrainedThickLayerNm', 'prelayerThickNm'}, 9);
    opts.prestrainedThinLayerNm = getOptionOrDefault(options, {'prestrainedThinLayerNm', 'prelayerThinNm'}, 2);
    opts.prestrainedTotalThicknessNm = getOptionOrDefault(options, {'prestrainedTotalThicknessNm', 'prelayerTotalThicknessNm'}, 120);
    opts.startPositionNm = getOptionOrDefault(options, {'startPositionNm', 'startPosition_nm'}, [0, 0, 0]);
    opts.rngSeed = getOptionOrDefault(options, {'rngSeed', 'seed'}, []);

    if isfield(options, 'Vpits') && ~isempty(options.Vpits)
        opts.Vpits = options.Vpits;
    else
        opts.Vpits = cell(0, 1);
    end

    opts.voltagesKeV = double(opts.voltagesKeV);
    opts.numElectrons = max(1, round(double(opts.numElectrons)));
    opts.numRepeats = max(1, round(double(opts.numRepeats)));
    opts.mqwPairs = max(1, round(double(opts.mqwPairs)));
    opts.mqwBarriers = max(0, round(double(opts.mqwBarriers)));
    opts.barrierThick = double(opts.barrierThick);
    opts.wellThick = double(opts.wellThick);
    opts.prestrainedPeriods = max(1, round(double(opts.prestrainedPeriods)));
    opts.prestrainedThickLayerNm = max(0.1, double(opts.prestrainedThickLayerNm));
    opts.prestrainedThinLayerNm = max(0.1, double(opts.prestrainedThinLayerNm));
    opts.prestrainedTotalThicknessNm = max(0.2, double(opts.prestrainedTotalThicknessNm));
    opts.startPositionNm = double(opts.startPositionNm);

    if any(~isfinite(opts.voltagesKeV)) || isempty(opts.voltagesKeV)
        error('precies:simulation:InvalidVoltages', 'voltagesKeV must be a non-empty numeric vector.');
    end
end

function value = getOptionOrDefault(options, fieldNames, defaultValue)
    value = defaultValue;
    for fieldIdx = 1:numel(fieldNames)
        fieldName = fieldNames{fieldIdx};
        if isfield(options, fieldName) && ~isempty(options.(fieldName))
            value = options.(fieldName);
            return;
        end
    end
end

function R = fresnelReflectance(n1_func, n2_func, lambda, cosTheta1, pol)
    n1 = n1_func(lambda);
    n2 = n2_func(lambda);
    sinTheta1 = sqrt(1 - cosTheta1^2);
    sinTheta2 = (n1/n2)*sinTheta1;
    
    if sinTheta2 > 1
        R = 1; 
        return;
    end
    
    cosTheta2 = sqrt(1 - sinTheta2^2);
    if pol == 0 
        Rs = ((n2*cosTheta1 - n1*cosTheta2)/(n2*cosTheta1 + n1*cosTheta2))^2;
        R = Rs;
    else 
        Rp = ((n1*cosTheta2 - n2*cosTheta1)/(n1*cosTheta2 + n2*cosTheta1))^2;
        R = Rp;
    end
end

function [layers, layerBoundaries] = Layers(layer_params,params)
mqw_pairs = max(1, round(params.mqw_pairs));
mqw_barriers = resolveMqwBarrierCount(params, mqw_pairs);
mqwTemplateCount = nnz(strcmp(layer_params(:, 1), 'MQW-Barrier'));
prestrainedPeriods = resolvePrestrainedPeriodCount(params);
prestrainedTemplateCount = nnz(strcmp(layer_params(:, 1), 'Prestrained'));
maxLayerRows = size(layer_params, 1) + ...
    mqwTemplateCount * max(0, mqw_pairs + mqw_barriers) + ...
    prestrainedTemplateCount * max(0, 2 * prestrainedPeriods - 1);
layers = cell(maxLayerRows, size(layer_params, 2));
layerWriteIdx = 0;
skip_mqw_well = false;
for i = 1:size(layer_params,1)
    if skip_mqw_well
        skip_mqw_well = false;
        continue; 
    end
    
    if strcmp(layer_params{i,1}, 'Prestrained') && prestrainedPeriods > 1
        preRows = buildPrestrainedSuperlatticeRows(layer_params(i,:), params, prestrainedPeriods);
        for preIdx = 1:size(preRows, 1)
            layerWriteIdx = layerWriteIdx + 1;
            layers(layerWriteIdx,:) = preRows(preIdx,:);
        end
    elseif strcmp(layer_params{i,1}, 'MQW-Barrier')
        barrierRow = layer_params(i,:);
        wellRow = layer_params(i + 1,:);
        if mqw_barriers >= mqw_pairs
            layerWriteIdx = layerWriteIdx + 1;
            layers(layerWriteIdx,:) = barrierRow;
            remainingBarriers = mqw_barriers - 1;
            for j = 1:mqw_pairs
                layerWriteIdx = layerWriteIdx + 1;
                layers(layerWriteIdx,:) = wellRow;
                if remainingBarriers > 0
                    layerWriteIdx = layerWriteIdx + 1;
                    layers(layerWriteIdx,:) = barrierRow;
                    remainingBarriers = remainingBarriers - 1;
                end
            end
        else
            layerWriteIdx = layerWriteIdx + 1;
            layers(layerWriteIdx,:) = wellRow;
            remainingWells = mqw_pairs - 1;
            remainingBarriers = mqw_barriers;
            while remainingWells > 0 || remainingBarriers > 0
                if remainingBarriers > 0
                    layerWriteIdx = layerWriteIdx + 1;
                    layers(layerWriteIdx,:) = barrierRow;
                    remainingBarriers = remainingBarriers - 1;
                end
                if remainingWells > 0
                    layerWriteIdx = layerWriteIdx + 1;
                    layers(layerWriteIdx,:) = wellRow;
                    remainingWells = remainingWells - 1;
                end
            end
        end
        
        skip_mqw_well = true; 
    else
        layerWriteIdx = layerWriteIdx + 1;
        layers(layerWriteIdx,:) = layer_params(i,:);
    end
end

layers = layers(1:layerWriteIdx, :);
layerBoundaries = zeros(1, size(layers,1)+1); 
layerBoundaries(1) = -layers{1,3}*1e-9;

for i = 1:size(layers,1)
    thickness = layers{i,3}*1e-9;
    layerBoundaries(i+1) = layerBoundaries(i) + thickness;
end
layerBoundaries(end+1) = Inf; 
end

function mqw_barriers = resolveMqwBarrierCount(params, mqw_pairs)
    if isstruct(params) && isfield(params, 'mqw_barriers') && ~isempty(params.mqw_barriers)
        mqw_barriers = params.mqw_barriers;
    elseif isstruct(params) && isfield(params, 'mqwBarrierCount') && ~isempty(params.mqwBarrierCount)
        mqw_barriers = params.mqwBarrierCount;
    else
        mqw_barriers = mqw_pairs - 1;
    end
    mqw_barriers = max(0, round(double(mqw_barriers)));
end

function prestrainedPeriods = resolvePrestrainedPeriodCount(params)
    prestrainedPeriods = resolvePrestrainedOption(params, ...
        {'prestrainedPeriods', 'prelayerPeriods', 'prePeriods'}, 1);
    prestrainedPeriods = max(1, round(double(prestrainedPeriods)));
end

function preRows = buildPrestrainedSuperlatticeRows(baseRow, params, prestrainedPeriods)
    thickLayerNm = resolvePrestrainedOption(params, ...
        {'prestrainedThickLayerNm', 'prestrainedThickNm', 'prelayerThickNm'}, 9);
    thinLayerNm = resolvePrestrainedOption(params, ...
        {'prestrainedThinLayerNm', 'prestrainedThinNm', 'prelayerThinNm'}, 2);
    totalThicknessNm = resolvePrestrainedOption(params, ...
        {'prestrainedTotalThicknessNm', 'prelayerTotalThicknessNm'}, ...
        prestrainedPeriods * (thickLayerNm + thinLayerNm));

    thickLayerNm = max(0.1, double(thickLayerNm));
    thinLayerNm = max(0.1, double(thinLayerNm));
    totalThicknessNm = max(0.2, double(totalThicknessNm));
    scaleFactor = totalThicknessNm / max(prestrainedPeriods * (thickLayerNm + thinLayerNm), eps);
    thickLayerNm = thickLayerNm * scaleFactor;
    thinLayerNm = thinLayerNm * scaleFactor;

    baseIn = double(baseRow{1, 2});
    thickIn = resolvePrestrainedOption(params, ...
        {'prestrainedThickInComposition', 'prelayerThickInComposition'}, ...
        max(0.01, min(baseIn, 0.85 * baseIn)));
    thinInDefault = baseIn;
    if thinLayerNm > 0
        thinInDefault = (baseIn * (thickLayerNm + thinLayerNm) - thickIn * thickLayerNm) / thinLayerNm;
    end
    thinInDefault = max(baseIn + 0.005, thinInDefault);
    thinIn = resolvePrestrainedOption(params, ...
        {'prestrainedThinInComposition', 'prelayerThinInComposition'}, thinInDefault);
    thickIn = clamp(double(thickIn), 0, 0.30);
    thinIn = clamp(double(thinIn), 0, 0.30);

    preRows = cell(2 * prestrainedPeriods, size(baseRow, 2));
    for periodIdx = 1:prestrainedPeriods
        thickRow = baseRow;
        thickRow{1, 1} = sprintf('Prestrained-Thick-%02d', periodIdx);
        thickRow{1, 2} = thickIn;
        thickRow{1, 3} = thickLayerNm;
        thinRow = baseRow;
        thinRow{1, 1} = sprintf('Prestrained-Thin-%02d', periodIdx);
        thinRow{1, 2} = thinIn;
        thinRow{1, 3} = thinLayerNm;

        preRows(2 * periodIdx - 1, :) = thickRow;
        preRows(2 * periodIdx, :) = thinRow;
    end
end

function value = resolvePrestrainedOption(params, fieldNames, defaultValue)
    value = defaultValue;
    if ~isstruct(params)
        return;
    end

    profile = getFieldOrDefault(params, 'calibrationProfile', struct());
    if isstruct(profile)
        for fieldIdx = 1:numel(fieldNames)
            fieldName = fieldNames{fieldIdx};
            if isfield(profile, fieldName) && ~isempty(profile.(fieldName))
                value = profile.(fieldName);
                return;
            end
        end
    end

    for fieldIdx = 1:numel(fieldNames)
        fieldName = fieldNames{fieldIdx};
        if isfield(params, fieldName) && ~isempty(params.(fieldName))
            value = params.(fieldName);
            return;
        end
    end
end

function plotRefractiveIndex(material_data)
    figure('Name','Refractive Index Analysis','Position',[200 200 1000 600]);
    hold on;
    style_config = {
        'Substrate',       'LineStyle', '-',  'Color', [0.7 0.7 0.7], 'LineWidth', 2.5;
        'p/n-GaN',         'LineStyle', '-',  'Color', [0.2 0.6 0.8], 'LineWidth', 2.5;
        'Prestrained',     'LineStyle', '--', 'Color', [0.9 0.4 0.1], 'LineWidth', 2.5;
        'MQW-Barrier',     'LineStyle', '-',  'Color', [0.1 0.8 0.2], 'LineWidth', 2.5;
        'MQW-Well',        'LineStyle', ':',  'Color', [0.8 0.1 0.8], 'LineWidth', 2.5;
        'p-EBL',           'LineStyle', '-',  'Color', [0.9 0.9 0.2], 'LineWidth', 2.5;
    };

    legend_entries = cell(numel(material_data), 1);
    legendCount = 0;
    for idx = 1:length(material_data)
        layer = material_data(idx);

        style_idx = find(strcmp(style_config(:,1), layer.Group));
        if isempty(style_idx), continue; end

        In_percent = sprintf('%.0f%%', layer.InComposition*100);
        leg_text = sprintf('%s (In=%s)', strrep(layer.Group,'-',' '), In_percent);
        
        plot(200:5:800, layer.RefractiveIndex,'LineStyle',...
            style_config{style_idx,3}, 'Color',style_config{style_idx,5},'LineWidth',style_config{style_idx,7},...
            'DisplayName', leg_text);
        
        legendCount = legendCount + 1;
        legend_entries{legendCount} = leg_text;
    end
    legend_entries = legend_entries(1:legendCount);

    title('Comprehensive Refractive Index Profile Analysis', 'FontSize',14);
    xlabel('Wavelength (nm)', 'FontSize',12);
    ylabel('Refractive Index (n)', 'FontSize',12);
    grid on;
    xlim([400 800]);
    legend(legend_entries, 'Location','eastoutside',...
        'FontSize',9, 'EdgeColor','none');
    set(gca, 'FontSize',11, 'LineWidth',1.2);
end

function plotAbsorption(material_data)
    figure('Name','Absorption Analysis','Position',[200 200 1000 600]);
    hold on;
    
    style_config = {
        'Substrate',       'LineStyle', '-',  'Color', [0.7 0.7 0.7], 'LineWidth', 2.5;
        'p/n-GaN',         'LineStyle', '-',  'Color', [0.2 0.6 0.8], 'LineWidth', 2.5;
        'Prestrained',     'LineStyle', '--', 'Color', [0.9 0.4 0.1], 'LineWidth', 2.5;
        'MQW-Barrier',     'LineStyle', '-',  'Color', [0.1 0.8 0.2], 'LineWidth', 2.5;
        'MQW-Well',        'LineStyle', ':',  'Color', [0.8 0.1 0.8], 'LineWidth', 2.5;
        'p-EBL',           'LineStyle', '-',  'Color', [0.9 0.9 0.2], 'LineWidth', 2.5;
    };

    legend_entries = cell(numel(material_data), 1);
    legendCount = 0;
    for idx = 1:length(material_data)
        layer = material_data(idx);
        
        style_idx = find(strcmp(style_config(:,1), layer.Group));
        if isempty(style_idx), continue; end
        
        In_percent = sprintf('%.0f%%', layer.InComposition*100);
        leg_text = sprintf('%s (In=%s)', strrep(layer.Group,'-',' '), In_percent);
        
        semilogy(200:5:800, layer.Absorption*1e-3,'LineStyle',...
            style_config{style_idx,3}, 'Color',style_config{style_idx,5},'LineWidth',style_config{style_idx,7},...
            'DisplayName', leg_text);
        
        legendCount = legendCount + 1;
        legend_entries{legendCount} = leg_text;
    end
    legend_entries = legend_entries(1:legendCount);
    
    title('Depth-Resolved Absorption Characteristics', 'FontSize',14);
    xlabel('Wavelength (nm)', 'FontSize',12);
    ylabel('Absorption Coefficient (mm^{-1})', 'FontSize',12);
    grid on;
    ylim([1e0 1e4]);
    xlim([400 800]);
    set(gca, 'YScale','log', 'FontSize',11, 'LineWidth',1.2);
    legend(legend_entries, 'Location','eastoutside',...
        'FontSize',9, 'EdgeColor','none');
end

function CL_sources = run_electron_simulation(layer_params, params, Eg_GaN, Eg_InN, start_pos_0, Vpits)
    numElectrons = max(1, round(double(getFieldOrDefault(params, 'numElectrons', 100))));
    [EXC_positions, electronStats] = simulateElectronExcitations( ...
        layer_params, params, Eg_GaN, Eg_InN, start_pos_0, Vpits, numElectrons);
    if ~isempty(EXC_positions) && ~isModelToggleEnabled(params, 'enablePEInteractionVolume', true)
        meanExcitationZ = mean(EXC_positions(:, 3), 'omitnan');
        EXC_positions(:, 1) = start_pos_0(1);
        EXC_positions(:, 2) = start_pos_0(2);
        EXC_positions(:, 3) = meanExcitationZ;
    end
    carrierTransport = getCarrierTransportDefaults(params);
    CL_sources = generate_CL_sources( ...
        EXC_positions, layer_params, params, ...
        carrierTransport.defaultLateralDiffusion_m, carrierTransport.defaultVerticalDiffusion_m, ...
        Eg_GaN, Eg_InN, Vpits, electronStats.nominalSurfaceZ_nm);
    if isempty(CL_sources)
        return;
    end
end

function [EXC_positions, electronStats] = simulateElectronExcitations(layer_params, params, Eg_GaN, Eg_InN, start_pos_0, Vpits, numElectrons)
    if nargin < 7 || isempty(numElectrons)
        numElectrons = 100;
    end

    E0 = params.electronEnergy;
    E_threshold = 50;
    [layers, layerBoundaries] = Layers(layer_params, params);
    nominalSurfaceZ = layerBoundaries(end - 1) * 1e9;
    z_min = layerBoundaries(1);
    layerProps = buildElectronLayerProperties(layers, Eg_GaN, Eg_InN);
    [launchPos, entryExcitationScale, entrySurfaceZ_m] = buildElectronLaunchEnsemble( ...
        start_pos_0, numElectrons, params, layers, layerBoundaries, Vpits, nominalSurfaceZ);
    nominalSurfaceZ_m = nominalSurfaceZ * 1e-9;
    pos = launchPos;
    dir = repmat([0, 0, -1], numElectrons, 1);
    E = E0 * ones(numElectrons, 1);
    accumulated_energy = zeros(numElectrons, 1);
    deepestZ = pos(:, 3);
    active = true(numElectrons, 1);
    useGpuElectrons = canUseGpuElectronBatching(numElectrons);
    if useGpuElectrons
        pos = gpuArray(pos);
        dir = gpuArray(dir);
        E = gpuArray(E);
        accumulated_energy = gpuArray(accumulated_energy);
        deepestZ = gpuArray(deepestZ);
        active = gpuArray(active);
        entryExcitationScale = gpuArray(entryExcitationScale);
        entrySurfaceZ_m = gpuArray(entrySurfaceZ_m);
    end
    maxIterations = 2000;
    excitationChunks = cell(maxIterations, 1);
    excitationChunkCount = 0;
    iterationCount = 0;

    while hasActiveRays(active) && iterationCount < maxIterations
        iterationCount = iterationCount + 1;
        activeIdx = gather(find(active & E > E_threshold));
        if isempty(activeIdx)
            break;
        end

        start_pos = pos(activeIdx, :);
        currentDir = dir(activeIdx, :);
        currentLayer = findLayerIndicesFast(gather(start_pos(:, 3)), layerBoundaries, size(layers, 1));
        validMask = currentLayer >= 1 & currentLayer <= size(layers, 1);
        if ~all(validMask)
            active(activeIdx(~validMask)) = false;
            activeIdx = activeIdx(validMask);
            start_pos = start_pos(validMask, :);
            currentDir = currentDir(validMask, :);
            currentLayer = currentLayer(validMask);
            if isempty(activeIdx)
                continue;
            end
        end

        J = layerProps.J(currentLayer);
        Z = layerProps.Z(currentLayer);
        A = layerProps.A(currentLayer);
        rho = layerProps.rho(currentLayer);
        Eg = layerProps.Eg(currentLayer);
        currentE = E(activeIdx);
        J = castNumericState(J, currentE);
        Z = castNumericState(Z, currentE);
        A = castNumericState(A, currentE);
        rho = castNumericState(rho, currentE);
        Eg = castNumericState(Eg, currentE);

        sigma = MottCrossSectionVectorized(currentE, Z);
        Lambda = A ./ (rho * 1e3 * 6.022e23 .* sigma * 1e-4);
        s = -Lambda .* log(randomLike(currentE, [numel(activeIdx), 1]));
        dEds = (-(rho * 1e-3 .* Z) ./ (J .* A)) .* ...
            (1e4 ./ (0.303 * sqrt(J ./ currentE) + 1.16 * sqrt(currentE ./ J) + 0.147 * (currentE ./ J)));
        delta_E = dEds .* s * 1e9;
        delta_E = max(delta_E, -(currentE - E_threshold));
        end_pos = start_pos + currentDir .* s;
        currentEntrySurfaceZ = entrySurfaceZ_m(activeIdx);

        escapeMask = end_pos(:, 3) > currentEntrySurfaceZ;
        if any(escapeMask)
            delta_z = currentEntrySurfaceZ(escapeMask) - start_pos(escapeMask, 3);
            dir_z = currentDir(escapeMask, 3);
            dir_z(abs(dir_z) < eps) = -eps;
            t_escape = delta_z ./ dir_z;
            end_pos(escapeMask, :) = start_pos(escapeMask, :) + currentDir(escapeMask, :) .* t_escape;
        end

        deepestZ(activeIdx) = min(deepestZ(activeIdx), end_pos(:, 3));
        accumulated_before = accumulated_energy(activeIdx);
        depositedEnergy = abs(delta_E) .* entryExcitationScale(activeIdx);
        accumulated_after = accumulated_before + depositedEnergy;
        num_excitons = floor(accumulated_after ./ (3 .* Eg));
        num_excitons_cpu = gather(num_excitons);
        if any(num_excitons_cpu > 0)
            excitationChunkCount = excitationChunkCount + 1;
            excitationChunks{excitationChunkCount} = buildExcitationChunk( ...
                gather(start_pos), gather(end_pos), gather(accumulated_before), gather(depositedEnergy), ...
                gather(Eg), num_excitons_cpu, gather(entryExcitationScale(activeIdx)));
        end
        accumulated_energy(activeIdx) = accumulated_after - num_excitons .* (3 .* Eg);

        pos(activeIdx, :) = end_pos;
        updatedE = currentE + delta_E;
        E(activeIdx) = updatedE;

        deactivateMask = escapeMask | updatedE <= E_threshold;
        active(activeIdx(deactivateMask)) = false;

        survivorMask = ~deactivateMask;
        if any(survivorMask)
            survivorIdx = activeIdx(survivorMask);
            theta = polarScatteringAngleVectorized(updatedE(survivorMask));
            phi = 2 * pi * randomLike(updatedE(survivorMask), [sum(survivorMask), 1]);
            dir(survivorIdx, :) = rotateVectorBatch(currentDir(survivorMask, :), theta, phi);
        end
    end

    if excitationChunkCount == 0
        EXC_positions = zeros(0, 4);
    else
        EXC_positions = vertcat(excitationChunks{1:excitationChunkCount});
        EXC_positions(:, 3) = max(min(EXC_positions(:, 3), nominalSurfaceZ_m), z_min);
    end

    penetrationDepths_nm = gather(max(0, (entrySurfaceZ_m - deepestZ) * 1e9));

    electronStats = struct( ...
        'penetrationDepths_nm', penetrationDepths_nm, ...
        'surfaceZ_nm', mean(gather(entrySurfaceZ_m)) * 1e9, ...
        'entrySurfaceZ_nm', gather(entrySurfaceZ_m) * 1e9, ...
        'meanEntryExcitationScale', mean(gather(entryExcitationScale)), ...
        'nominalSurfaceZ_nm', nominalSurfaceZ, ...
        'layerBoundaries_nm', layerBoundaries * 1e9, ...
        'layerNames', {layers(:, 1)});
end

function [launchPos, entryExcitationScale, entrySurfaceZ_m] = buildElectronLaunchEnsemble(start_pos_0, numElectrons, params, layers, layerBoundaries, Vpits, nominalSurfaceZ_nm)
    defaults = getElectronEntryDefaults(params);
    baseXY_nm = start_pos_0(1:2) * 1e9;
    launchXY_nm = baseXY_nm + defaults.beamSpotSigma_nm * randn(numElectrons, 2);
    launchPos = zeros(numElectrons, 3);
    launchPos(:, 1:2) = launchXY_nm * 1e-9;
    entryExcitationScale = ones(numElectrons, 1);
    entrySurfaceZ_m = nominalSurfaceZ_nm * 1e-9 * ones(numElectrons, 1);
    incidentDir = [0, 0, -1];
    materialProbeOffset_nm = 0.25;

    for electronIdx = 1:numElectrons
        [inVpit, surfaceZ_nm] = checkVpit( ...
            launchXY_nm(electronIdx, 1), launchXY_nm(electronIdx, 2), nominalSurfaceZ_nm, Vpits, nominalSurfaceZ_nm);
        launchPos(electronIdx, 3) = surfaceZ_nm * 1e-9;
        entrySurfaceZ_m(electronIdx) = surfaceZ_nm * 1e-9;
        if ~inVpit
            continue;
        end

        probePoint_nm = [launchXY_nm(electronIdx, :), max(layerBoundaries(1) * 1e9, surfaceZ_nm - materialProbeOffset_nm)];
        entryState = precies.vpit('resolvePointState', probePoint_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ_nm);
        incidenceCos = resolveElectronEntryIncidenceCos(entryState, incidentDir);
        entryProfile = resolveVpitEntryProfile(launchXY_nm(electronIdx, :), surfaceZ_nm, Vpits, nominalSurfaceZ_nm);
        depthFraction = entryProfile.depthFraction;

        sidewallScale = defaults.sidewallYieldFloor + ...
            (1 - defaults.sidewallYieldFloor) * incidenceCos ^ defaults.sidewallYieldExponent;
        apexScale = 1 - defaults.apexSuppressionStrength * depthFraction ^ defaults.apexSuppressionExponent;
        rimEnhancement = computeVpitRimEntryEnhancement(entryProfile, defaults);
        centerPenalty = computeVpitCoreEntryPenalty(entryProfile, defaults);
        centerRefill = computeVpitCenterRefillEnhancement(entryProfile, entryState, defaults);
        if strcmp(entryState.materialKind, 'semipolar_shell')
            materialScale = defaults.semipolarEntryScale;
        else
            materialScale = 1;
        end

        entryExcitationScale(electronIdx) = max( ...
            defaults.minEntryScale, ...
            min(defaults.maxEntryScale, sidewallScale * apexScale * materialScale * rimEnhancement * centerPenalty * centerRefill));
    end
end

function defaults = getElectronEntryDefaults(params)
    scanX_nm = 500;
    scanY_nm = 500;
    totalPixels = 225;
    electronEnergy_eV = 5e3;
    if isstruct(params)
        if isfield(params, 'scanX') && isfinite(params.scanX)
            scanX_nm = params.scanX;
        end
        if isfield(params, 'scanY') && isfinite(params.scanY)
            scanY_nm = params.scanY;
        end
        if isfield(params, 'totalPixels') && isfinite(params.totalPixels)
            totalPixels = params.totalPixels;
        elseif isfield(params, 'pixelNum') && isfinite(params.pixelNum)
            totalPixels = params.pixelNum;
        end
        if isfield(params, 'electronEnergy') && isfinite(params.electronEnergy)
            electronEnergy_eV = params.electronEnergy;
        end
    end

    equivalentPixelPitch_nm = sqrt(max(scanX_nm * scanY_nm / max(totalPixels, 1), 1));
    electronEnergy_keV = max(1, electronEnergy_eV / 1e3);
    beamSpotSigma_nm = clamp(0.32 * equivalentPixelPitch_nm + 1.8 * sqrt(electronEnergy_keV), 4, 18);
    defaults = struct( ...
        'beamSpotSigma_nm', beamSpotSigma_nm, ...
        'sidewallYieldFloor', 0.66, ...
        'sidewallYieldExponent', 1.35, ...
        'apexSuppressionStrength', 0.11, ...
        'apexSuppressionExponent', 1.28, ...
        'semipolarEntryScale', 1.34, ...
        'minEntryScale', 0.54, ...
        'maxEntryScale', 1.95, ...
        'rimEnhancementStrength', 0.58, ...
        'rimEnhancementCenterRatio', 0.82, ...
        'rimEnhancementWidthRatio', 0.18, ...
        'rimDepthCenter', 0.26, ...
        'rimDepthWidth', 0.20, ...
        'centerPenaltyStrength', 0.06, ...
        'centerPenaltyWidthRatio', 0.18, ...
        'centerPenaltyDepthWeight', 0.22, ...
        'centerRefillStrength', 0.46, ...
        'centerRefillWidthRatio', 0.42, ...
        'centerRefillDepthCenter', 0.50, ...
        'centerRefillDepthWidth', 0.28, ...
        'centerRefillSemipolarBonus', 0.25);
end

function incidenceCos = resolveElectronEntryIncidenceCos(state, incidentDir)
    incidenceCos = 1;
    if isfield(state, 'closestFacetNormal3D') && numel(state.closestFacetNormal3D) == 3 && all(isfinite(state.closestFacetNormal3D))
        normal = state.closestFacetNormal3D(:)';
        if norm(normal) > eps
            normal = normal / norm(normal);
            incidenceCos = abs(dot(incidentDir, normal));
        end
    end
    incidenceCos = max(0.2, min(1, incidenceCos));
end

function entryProfile = resolveVpitEntryProfile(sampleXY_nm, surfaceZ_nm, Vpits, nominalSurfaceZ_nm)
    entryProfile = struct( ...
        'isInside', false, ...
        'vpitIndex', 0, ...
        'depthFraction', 0, ...
        'radialRatio', 1, ...
        'radialDistance_nm', inf);
    depth_nm = max(0, nominalSurfaceZ_nm - surfaceZ_nm);
    bestConsistency = inf;
    for pitIdx = 1:numel(Vpits)
        pit = Vpits{pitIdx};
        if isempty(pit) || ~isstruct(pit) || ~isfield(pit, 'depth') || pit.depth <= 0
            continue;
        end

        dx = sampleXY_nm(1) - pit.center(1);
        dy = sampleXY_nm(2) - pit.center(2);
        radialDistance_nm = hypot(dx, dy);
        orientationDeg = 0;
        if isfield(pit, 'orientationDeg') && ~isempty(pit.orientationDeg) && isfinite(pit.orientationDeg)
            orientationDeg = mod(pit.orientationDeg, 60);
        end
        facetNormalsDeg = orientationDeg + (0:60:300);
        hexDistance_nm = max(dx .* cosd(facetNormalsDeg) + dy .* sind(facetNormalsDeg));
        if hexDistance_nm > pit.topRadius + eps || depth_nm > pit.depth + eps
            continue;
        end

        radialRatio = clamp(hexDistance_nm / max(pit.topRadius, eps), 0, 1);
        depthFraction = clamp(depth_nm / max(pit.depth, eps), 0, 1);
        surfaceConsistency = abs((1 - radialRatio) - depthFraction);
        if surfaceConsistency > bestConsistency
            continue;
        end

        bestConsistency = surfaceConsistency;
        entryProfile.isInside = true;
        entryProfile.vpitIndex = pitIdx;
        entryProfile.depthFraction = depthFraction;
        entryProfile.radialRatio = radialRatio;
        entryProfile.radialDistance_nm = radialDistance_nm;
    end
end

function rimEnhancement = computeVpitRimEntryEnhancement(entryProfile, defaults)
    rimEnhancement = 1;
    if ~isstruct(entryProfile) || ~getFieldOrDefault(entryProfile, 'isInside', false)
        return;
    end

    radialRatio = clamp(getFieldOrDefault(entryProfile, 'radialRatio', 1), 0, 1);
    depthFraction = clamp(getFieldOrDefault(entryProfile, 'depthFraction', 0), 0, 1);
    rimKernel = exp(-0.5 * ...
        ((radialRatio - defaults.rimEnhancementCenterRatio) / max(defaults.rimEnhancementWidthRatio, 0.05)) ^ 2);
    depthKernel = exp(-0.5 * ...
        ((depthFraction - defaults.rimDepthCenter) / max(defaults.rimDepthWidth, 0.05)) ^ 2);
    rimEnhancement = 1 + defaults.rimEnhancementStrength * rimKernel * (0.65 + 0.35 * depthKernel);
end

function centerPenalty = computeVpitCoreEntryPenalty(entryProfile, defaults)
    centerPenalty = 1;
    if ~isstruct(entryProfile) || ~getFieldOrDefault(entryProfile, 'isInside', false)
        return;
    end

    radialRatio = clamp(getFieldOrDefault(entryProfile, 'radialRatio', 1), 0, 1);
    depthFraction = clamp(getFieldOrDefault(entryProfile, 'depthFraction', 0), 0, 1);
    centerKernel = exp(-0.5 * (radialRatio / max(defaults.centerPenaltyWidthRatio, 0.05)) ^ 2);
    depthWeight = 0.45 + defaults.centerPenaltyDepthWeight * depthFraction;
    centerPenalty = 1 - defaults.centerPenaltyStrength * centerKernel * depthWeight;
    centerPenalty = max(0.88, min(1.0, centerPenalty));
end

function centerRefill = computeVpitCenterRefillEnhancement(entryProfile, entryState, defaults)
    centerRefill = 1;
    if ~isstruct(entryProfile) || ~getFieldOrDefault(entryProfile, 'isInside', false)
        return;
    end

    radialRatio = clamp(getFieldOrDefault(entryProfile, 'radialRatio', 1), 0, 1);
    depthFraction = clamp(getFieldOrDefault(entryProfile, 'depthFraction', 0), 0, 1);
    centerKernel = exp(-0.5 * (radialRatio / max(defaults.centerRefillWidthRatio, 0.08)) ^ 2);
    depthKernel = exp(-0.5 * ...
        ((depthFraction - defaults.centerRefillDepthCenter) / max(defaults.centerRefillDepthWidth, 0.08)) ^ 2);

    materialBonus = 1;
    if isstruct(entryState) && strcmp(getFieldOrDefault(entryState, 'materialKind', ''), 'semipolar_shell')
        materialBonus = materialBonus + defaults.centerRefillSemipolarBonus;
    elseif isstruct(entryState) && getFieldOrDefault(entryState, 'isMQW', false)
        materialBonus = materialBonus + 0.5 * defaults.centerRefillSemipolarBonus;
    end

    centerRefill = 1 + defaults.centerRefillStrength * centerKernel * depthKernel * materialBonus;
    centerRefill = min(centerRefill, 1.45);
end

function [J_compound, Z_bar, A_bar, rho] = calculate_compound_parameters(x)
    Z_elements = [49, 31, 7];
    A_elements = [114.82, 69.72, 14.007];
    C = [x, 1-x, 1];
    C = C / sum(C);
    J = zeros(1,3);
    for i = 1:3
        if Z_elements(i) >= 13
            J(i) = 9.76*Z_elements(i) + 58.5*Z_elements(i)^(-0.19);
        else
            J(i) = 11.5*Z_elements(i);
        end
    end
    J_compound = exp((sum((C.*Z_elements./A_elements).*log(J)))/sum(C.*Z_elements./A_elements));
    Z_bar = sum(C.*Z_elements./A_elements) / sum(C./A_elements);
    A_bar = 1 / sum(C./A_elements);
    rho = x*6810 + (1-x)*6150;
end

function layerProps = buildElectronLayerProperties(layers, Eg_GaN, Eg_InN)
    numLayers = size(layers, 1);
    layerProps = struct( ...
        'J', zeros(numLayers, 1), ...
        'Z', zeros(numLayers, 1), ...
        'A', zeros(numLayers, 1), ...
        'rho', zeros(numLayers, 1), ...
        'Eg', zeros(numLayers, 1));
    for layerIdx = 1:numLayers
        x = layers{layerIdx, 2};
        [layerProps.J(layerIdx), layerProps.Z(layerIdx), layerProps.A(layerIdx), layerProps.rho(layerIdx)] = ...
            calculate_compound_parameters(x);
        layerProps.Eg(layerIdx) = get_bandgap(layerIdx, layers, Eg_GaN, Eg_InN);
    end
end

function sigma = MottCrossSectionVectorized(E, Z)
    E_kev = E ./ 1e3;
    alpha = (3.4e-3) .* (Z .^ 0.67) ./ E_kev;
    beta = 26.42 ./ (Z .^ 1.42);
    lambdaValue = 1.162 + (1.28e-2) .* Z;
    term1 = (5.21e-21) .* (Z .^ 2) ./ (E_kev .^ 2);
    term2 = 4 * pi .* lambdaValue ./ (alpha .* (alpha + 1));
    term3 = (1 - exp(-beta .* sqrt(E_kev))) .* ((E_kev + 511) ./ (E_kev + 1022)) .^ 2;
    sigma = term1 .* term2 .* term3;
end

function CL_light_sources = generate_CL_sources(EXC_positions,layer_params,params, L_diff_xy, L_diff_z, Eg_GaN, Eg_InN, Vpits, nominalSurfaceZ)
    num_excitons = size(EXC_positions,1);
    vpitDefaults = precies.vpit('getDefaults');
    carrierTransport = getCarrierTransportDefaults(params);
    lateralSigma_nm = max(L_diff_xy * 1e9, carrierTransport.defaultLateralDiffusion_nm);
    verticalSigma_nm = max(L_diff_z * 1e9, carrierTransport.defaultVerticalDiffusion_nm);

    [layers, layerBoundaries] = Layers(layer_params,params);
    valid_sources = zeros(max(1, 2 * num_excitons), 9);
    validSourceCount = 0;

    for i = 1:num_excitons
        startPos_nm = EXC_positions(i, 1:3) * 1e9;
        sourceWeight = 1;
        if size(EXC_positions, 2) >= 4 && isfinite(EXC_positions(i, 4))
            sourceWeight = EXC_positions(i, 4);
        end
        startState = precies.vpit('resolvePointState', startPos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
        if ~isCarrierStateAccessible(startState)
            [valid_sources, validSourceCount] = appendVpitCavitySidewallEmissionSource( ...
                valid_sources, validSourceCount, startPos_nm, startState, sourceWeight, ...
                layers, layerBoundaries, Vpits, nominalSurfaceZ, Eg_GaN, Eg_InN, carrierTransport);
            continue;
        end

        directFraction = resolveShallowDirectRecombinationFraction( ...
            startPos_nm, startState, nominalSurfaceZ, carrierTransport, Vpits);
        if directFraction > 1e-4
            directBandgap_eV = getCarrierStateBandgap(startState, layers, Eg_GaN, Eg_InN);
            directSourceType = classifyShallowDirectEmissionSource(startState);
            directRadiativeYield = computeCarrierRadiativeYield( ...
                startPos_nm, startState, startPos_nm, startState, Vpits, carrierTransport);
            if directSourceType > 0 && directRadiativeYield > 0.02
                [valid_sources, validSourceCount] = appendCarrierEmissionSource( ...
                    valid_sources, validSourceCount, startPos_nm, directBandgap_eV, ...
                    0, directSourceType, 0, sourceWeight * directFraction, directRadiativeYield);
            end
        end

        transportedSourceWeight = sourceWeight * max(0, 1 - directFraction);
        if transportedSourceWeight <= 1e-4
            continue;
        end

        pathResult = propagateCarrierDiffusionPath( ...
            startPos_nm, lateralSigma_nm, verticalSigma_nm, layers, layerBoundaries, ...
            Vpits, nominalSurfaceZ, Eg_GaN, Eg_InN, carrierTransport);
        if ~pathResult.isValid
            continue;
        end

        candidatePos_nm = pathResult.finalPos_nm;
        candidateState = pathResult.finalState;
        if strcmp(candidateState.materialKind, 'cavity')
            continue;
        end

        originPos_nm = candidatePos_nm;
        originState = candidateState;
        finalPos_nm = candidatePos_nm;
        finalState = candidateState;
        shouldProjectToCPlane = false;
        if strcmp(candidateState.materialKind, 'semipolar_shell')
            directRecombinationProbability = resolveSemipolarDirectRecombinationProbability( ...
                candidateState, carrierTransport);
            shouldProjectToCPlane = rand() > directRecombinationProbability;
        elseif strcmp(candidateState.materialKind, 'bulk') && candidateState.isMQW && ...
                isfinite(candidateState.closestFacetDistance_nm) && ...
                candidateState.closestFacetDistance_nm <= vpitDefaults.carrierCaptureDistance_nm
            captureProb = exp(-candidateState.closestFacetDistance_nm / vpitDefaults.carrierCaptureDistance_nm);
            shouldProjectToCPlane = rand() < captureProb;
        end

        if shouldProjectToCPlane
            finalPos_nm = precies.vpit('projectCarrierSource', candidatePos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ, candidateState);
            projectedState = precies.vpit('resolvePointState', finalPos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
            if ~strcmp(projectedState.materialKind, 'bulk') || ~projectedState.isMQW
                continue;
            end
            finalState = projectedState;
        end

        emissionBandgap_eV = getCarrierStateBandgap(finalState, layers, Eg_GaN, Eg_InN);
        emissionHotShift_eV = computeCarrierEmissionHotShift( ...
            pathResult.remainingEnergy_eV, pathResult.pathLength_nm, finalState, carrierTransport);
        emissionSourceType = classifyCarrierEmissionSource(finalState);
        radiativeYield = computeCarrierRadiativeYield( ...
            originPos_nm, originState, finalPos_nm, finalState, Vpits, carrierTransport);
        if radiativeYield <= 0.02
            continue;
        end

        [valid_sources, validSourceCount] = appendCarrierEmissionSource( ...
            valid_sources, validSourceCount, finalPos_nm, emissionBandgap_eV, ...
            emissionHotShift_eV, emissionSourceType, pathResult.pathLength_nm, ...
            transportedSourceWeight, radiativeYield);
    end
    CL_light_sources = valid_sources(1:validSourceCount, :);
end

function [valid_sources, validSourceCount] = appendCarrierEmissionSource( ...
    valid_sources, validSourceCount, sourcePos_nm, emissionBandgap_eV, emissionHotShift_eV, ...
    emissionSourceType, pathLength_nm, sourceWeight, radiativeYield)

    if emissionSourceType <= 0 || sourceWeight <= 0 || radiativeYield <= 0
        return;
    end
    if validSourceCount + 1 > size(valid_sources, 1)
        valid_sources = growMatrixRows(valid_sources, validSourceCount + 1);
    end

    validSourceCount = validSourceCount + 1;
    valid_sources(validSourceCount, :) = [ ...
        sourcePos_nm * 1e-9, ...
        emissionBandgap_eV, ...
        emissionHotShift_eV, ...
        emissionSourceType, ...
        pathLength_nm, ...
        sourceWeight, ...
        radiativeYield];
end

function [valid_sources, validSourceCount] = appendVpitCavitySidewallEmissionSource( ...
    valid_sources, validSourceCount, startPos_nm, startState, sourceWeight, ...
    layers, layerBoundaries, Vpits, nominalSurfaceZ, Eg_GaN, Eg_InN, transportDefaults)

    if ~transportDefaults.enableVpitShortwaveEnhancement || ...
            ~strcmp(getFieldOrDefault(startState, 'materialKind', ''), 'cavity') || ...
            ~isfield(startState, 'vpitIndex') || startState.vpitIndex < 1 || startState.vpitIndex > numel(Vpits)
        return;
    end

    vpit = Vpits{startState.vpitIndex};
    if isempty(vpit) || ~isstruct(vpit) || ~isfield(vpit, 'surfaceZ_nm') || ...
            ~isfield(vpit, 'depth') || ~isfinite(vpit.depth) || vpit.depth <= 0
        return;
    end

    projectedPos_nm = precies.vpit('projectCarrierSource', ...
        startPos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ, startState);
    projectedState = precies.vpit('resolvePointState', projectedPos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
    if ~isCarrierStateAccessible(projectedState)
        return;
    end

    depthBelowSurface_nm = vpit.surfaceZ_nm - startPos_nm(3);
    if depthBelowSurface_nm < -1 || depthBelowSurface_nm > vpit.depth + 2
        return;
    end

    lowVoltageWeight = clamp( ...
        (transportDefaults.shallowDirectVoltagePivot_keV - transportDefaults.electronEnergy_keV) / ...
        max(transportDefaults.shallowDirectVoltageWidth_keV, eps), 0, 1);
    if lowVoltageWeight <= 0
        return;
    end

    facetDistance_nm = getFieldOrDefault(startState, 'closestFacetDistance_nm', inf);
    shellThickness_nm = getFieldOrDefault(startState, 'layerShellThickness_nm', 0);
    facetWidth_nm = max(transportDefaults.vpitFacetShortwaveWidth_nm, 1.5 * max(shellThickness_nm, 1));
    if isfinite(facetDistance_nm)
        facetKernel = exp(-max(facetDistance_nm, 0) / max(facetWidth_nm, eps));
    else
        facetKernel = 0.35;
    end

    depthNorm = clamp(depthBelowSurface_nm / max(vpit.depth, eps), 0, 1);
    depthKernel = 0.35 + 0.65 * exp(-0.5 * ...
        ((depthNorm - transportDefaults.vpitCavitySidewallDepthCenterNorm) / ...
        max(transportDefaults.vpitCavitySidewallDepthWidthNorm, 0.05)) ^ 2);
    channelScale = min(1, resolveVpitShortwaveChannelCoupling( ...
        projectedPos_nm, projectedState, transportDefaults, Vpits) / 2.35);
    sidewallFraction = transportDefaults.vpitCavitySidewallDirectScale * ...
        lowVoltageWeight * (0.35 + 0.65 * facetKernel) * depthKernel * (0.50 + 0.50 * channelScale);
    sidewallFraction = clamp(sidewallFraction, 0, transportDefaults.vpitCavitySidewallMaxFraction);
    if sidewallFraction <= 1e-4
        return;
    end

    highBandgapShare = clamp(transportDefaults.vpitHighBandgapSidewallChannelShare * ...
        (0.55 + 0.45 * lowVoltageWeight) * (0.70 + 0.30 * depthKernel), 0, 0.92);
    [highBandgapPos_nm, highBandgapState, hasHighBandgapSite, highBandgapLocalBandgap_eV] = ...
        resolveVpitGuidedHighBandgapEmissionSite( ...
        startPos_nm, projectedPos_nm, projectedState, vpit, layers, layerBoundaries, ...
        nominalSurfaceZ, Eg_GaN, Eg_InN, transportDefaults);
    if hasHighBandgapSite
        highBandgapFraction = clamp(sidewallFraction * highBandgapShare * ...
            transportDefaults.vpitHighBandgapSidewallChannelScale, 0, transportDefaults.vpitHighBandgapSidewallMaxFraction);
        if highBandgapFraction > 1e-4
            highBandgapRadiativeYield = max(0.24, computeCarrierRadiativeYield( ...
                startPos_nm, startState, highBandgapPos_nm, highBandgapState, Vpits, transportDefaults));
            highBandgapRadiativeYield = highBandgapRadiativeYield * ...
                computeHighBandgapSidewallRadiativeRetention(highBandgapState, transportDefaults);
            [valid_sources, validSourceCount] = appendCarrierEmissionSource( ...
                valid_sources, validSourceCount, highBandgapPos_nm, highBandgapLocalBandgap_eV, ...
                0, 1, norm(highBandgapPos_nm - startPos_nm), ...
                sourceWeight * highBandgapFraction, highBandgapRadiativeYield);
        end
        sidewallFraction = sidewallFraction * (1 - 0.82 * highBandgapShare);
    end

    sourceType = classifyShallowDirectEmissionSource(projectedState);
    if sourceType <= 0
        sourceType = classifyCarrierEmissionSource(projectedState);
    end
    if ~ismember(sourceType, [1, 2, 3, 6])
        sourceType = 6;
    end

    emissionBandgap_eV = NaN;
    if sourceType ~= 6
        emissionBandgap_eV = getCarrierStateBandgap(projectedState, layers, Eg_GaN, Eg_InN);
        if ~isfinite(emissionBandgap_eV)
            return;
        end
    end
    emissionHotShift_eV = min(0.035, 0.010 + 0.020 * lowVoltageWeight * depthKernel);
    radiativeYield = max(0.08, computeCarrierRadiativeYield( ...
        startPos_nm, startState, projectedPos_nm, projectedState, Vpits, transportDefaults));

    [valid_sources, validSourceCount] = appendCarrierEmissionSource( ...
        valid_sources, validSourceCount, projectedPos_nm, emissionBandgap_eV, ...
        emissionHotShift_eV, sourceType, norm(projectedPos_nm - startPos_nm), ...
        sourceWeight * sidewallFraction, radiativeYield);
end

function [targetPos_nm, targetState, isValid, localBandgap_eV] = resolveVpitGuidedHighBandgapEmissionSite( ...
    startPos_nm, projectedPos_nm, projectedState, vpit, layers, layerBoundaries, ...
    nominalSurfaceZ, Eg_GaN, Eg_InN, transportDefaults)

    targetPos_nm = projectedPos_nm;
    targetState = projectedState;
    isValid = isCarrierStateAccessible(targetState);
    localBandgap_eV = NaN;

    if isValid && shouldUseSidewallHighBandgapSite(targetState)
        localBandgap_eV = resolveEffectiveHighBandgapChannelEnergy( ...
            targetState, layers, Eg_GaN, Eg_InN, transportDefaults);
        return;
    end

    nGaNLayerIdx = findLayerIndexByName(layers, 'n-GaN');
    if nGaNLayerIdx < 1 || nGaNLayerIdx + 1 > numel(layerBoundaries)
        return;
    end

    nBottom_nm = layerBoundaries(nGaNLayerIdx) * 1e9;
    nTop_nm = layerBoundaries(nGaNLayerIdx + 1) * 1e9;
    if ~(isfinite(nBottom_nm) && isfinite(nTop_nm)) || nTop_nm <= nBottom_nm
        return;
    end

    nTopDepth_nm = nominalSurfaceZ - nTop_nm;
    entryDepth_nm = max(0, nominalSurfaceZ - startPos_nm(3));
    reachesNganAdjacentRegion = max(entryDepth_nm, vpit.depth) + ...
        transportDefaults.vpitHighBandgapNganAdjacencyMargin_nm >= nTopDepth_nm;
    if ~reachesNganAdjacentRegion
        return;
    end

    guidedOffset_nm = clamp(transportDefaults.vpitHighBandgapGuidedDepthBelowNTop_nm, 2, 35);
    zTarget_nm = min(nTop_nm - guidedOffset_nm, nTop_nm - 1);
    zTarget_nm = max(nBottom_nm + 1, zTarget_nm);
    targetXY_nm = 0.62 * startPos_nm(1:2) + 0.38 * vpit.center(:)';
    lateralJitter_nm = transportDefaults.vpitHighBandgapGuidedLateralJitter_nm;
    if isfinite(lateralJitter_nm) && lateralJitter_nm > 0
        targetXY_nm = targetXY_nm + lateralJitter_nm * randn(1, 2);
    end
    targetPos_nm = [targetXY_nm, zTarget_nm];

    targetState = precies.vpit('resolvePointState', targetPos_nm, layers, layerBoundaries, {vpit}, nominalSurfaceZ);
    if isCarrierStateAccessible(targetState) && contains(string(targetState.layerName), 'n-GaN')
        isValid = true;
        localBandgap_eV = resolveEffectiveHighBandgapChannelEnergy( ...
            targetState, layers, Eg_GaN, Eg_InN, transportDefaults);
    end
end

function localBandgap_eV = resolveEffectiveHighBandgapChannelEnergy(state, layers, Eg_GaN, Eg_InN, transportDefaults)
    baseBandgap_eV = getCarrierStateBandgap(state, layers, Eg_GaN, Eg_InN);
    if ~isfinite(baseBandgap_eV)
        localBandgap_eV = NaN;
        return;
    end

    layerName = string(getFieldOrDefault(state, 'layerName', ''));
    materialKind = string(getFieldOrDefault(state, 'materialKind', ''));
    inDepletionShift_eV = transportDefaults.vpitInCompositionDepletionEnergyShift_eV;
    qcseShift_eV = transportDefaults.vpitQcseScreeningEnergyShift_eV;
    strainShift_eV = transportDefaults.vpitStrainRelaxationEnergyShift_eV;
    if strcmp(materialKind, "semipolar_shell")
        sidewallFactor = 1.0;
    elseif contains(layerName, 'MQW-Well') || contains(layerName, 'MQW-Barrier')
        sidewallFactor = 0.82;
    elseif contains(layerName, 'Prestrained')
        sidewallFactor = 0.62;
    elseif contains(layerName, 'GaN')
        sidewallFactor = 0.25;
    else
        sidewallFactor = 0.45;
    end

    effectiveShift_eV = sidewallFactor * (inDepletionShift_eV + qcseShift_eV + strainShift_eV);
    localBandgap_eV = clamp(baseBandgap_eV + effectiveShift_eV, 2.95, Eg_GaN + 0.09);
end

function retention = computeHighBandgapSidewallRadiativeRetention(state, transportDefaults)
    layerName = string(getFieldOrDefault(state, 'layerName', ''));
    materialKind = string(getFieldOrDefault(state, 'materialKind', ''));
    nonradiativeScale = clamp(transportDefaults.vpitSidewallNonradiativeScale, 0, 3);
    if nonradiativeScale <= 0
        retention = 1;
        return;
    end

    if strcmp(materialKind, "semipolar_shell")
        sensitivity = 0.42;
    elseif contains(layerName, 'MQW-Well') || contains(layerName, 'MQW-Barrier')
        sensitivity = 0.56;
    elseif contains(layerName, 'Prestrained')
        sensitivity = 0.50;
    elseif contains(layerName, 'p-GaN') || contains(layerName, 'p-EBL')
        sensitivity = 0.38;
    else
        sensitivity = 0.26;
    end
    retention = clamp(1 - 0.30 * nonradiativeScale * sensitivity, 0.18, 1.0);
end

function tf = shouldUseSidewallHighBandgapSite(state)
    if ~isCarrierStateAccessible(state)
        tf = false;
        return;
    end

    materialKind = string(getFieldOrDefault(state, 'materialKind', ''));
    layerName = string(getFieldOrDefault(state, 'layerName', ''));
    tf = strcmp(materialKind, "semipolar_shell") || ...
        contains(layerName, 'MQW-Barrier') || ...
        contains(layerName, 'Prestrained') || ...
        contains(layerName, 'p-GaN') || ...
        contains(layerName, 'p-EBL') || ...
        contains(layerName, 'GaN');
end

function layerIdx = findLayerIndexByName(layers, layerNamePattern)
    layerIdx = 0;
    if isempty(layers)
        return;
    end
    names = string(layers(:, 1));
    matchIdx = find(contains(names, string(layerNamePattern), 'IgnoreCase', true), 1, 'first');
    if ~isempty(matchIdx)
        layerIdx = matchIdx;
    end
end

function selectedSources = ensureVpitShortwaveSourcesRepresented(selectedSources, allSources, numRays, params, Vpits, nominalSurfaceZ)
    if isempty(selectedSources) || isempty(allSources) || size(allSources, 2) < 6 || isempty(Vpits) || ...
            ~isModelToggleEnabled(params, 'enableVpitTransportCorrection', true)
        return;
    end

    selectedSources = selectedSources(1:min(size(selectedSources, 1), max(numRays, 1)), :);
    sourceTypesAll = round(allSources(:, 6));
    shortwaveMask = ismember(sourceTypesAll, [1, 2, 3, 6]) & ...
        classifySourcesNearVpits(allSources, Vpits, nominalSurfaceZ);
    if ~any(shortwaveMask)
        return;
    end

    sourceTypesSelected = round(selectedSources(:, 6));
    currentShortwaveCount = nnz(ismember(sourceTypesSelected, [1, 2, 3, 6]) & ...
        classifySourcesNearVpits(selectedSources, Vpits, nominalSurfaceZ));
    electronEnergy_keV = max(1, getFieldOrDefault(params, 'electronEnergy', 5e3) / 1e3);
    lowVoltageWeight = clamp((7.0 - electronEnergy_keV) / 3.0, 0, 1);
    targetShortwaveFraction = 0.14 + 0.26 * lowVoltageWeight;
    targetShortwaveCount = min(nnz(shortwaveMask), max(1, ceil(targetShortwaveFraction * max(numRays, 1))));
    missingCount = targetShortwaveCount - currentShortwaveCount;
    if missingCount <= 0
        return;
    end

    shortwaveSources = allSources(shortwaveMask, :);
    sourceWeights = ones(size(shortwaveSources, 1), 1);
    if size(shortwaveSources, 2) >= 8
        sourceWeights = max(0, shortwaveSources(:, 8));
    end
    [~, order] = sort(sourceWeights, 'descend');
    shortwaveSources = shortwaveSources(order, :);
    addCount = min(missingCount, size(shortwaveSources, 1));

    longwaveMask = ~ismember(sourceTypesSelected, [1, 2, 3, 6]);
    replaceIdx = find(longwaveMask, addCount, 'last');
    if numel(replaceIdx) < addCount
        fallbackIdx = setdiff((1:size(selectedSources, 1))', replaceIdx(:), 'stable');
        replaceIdx = [replaceIdx(:); fallbackIdx(1:min(addCount - numel(replaceIdx), numel(fallbackIdx)))];
    end
    if isempty(replaceIdx)
        return;
    end
    replaceCount = min(numel(replaceIdx), addCount);
    selectedSources(replaceIdx(1:replaceCount), :) = shortwaveSources(1:replaceCount, :);
end

function nearMask = classifySourcesNearVpits(sourceData, Vpits, nominalSurfaceZ)
    nearMask = false(size(sourceData, 1), 1);
    if isempty(sourceData) || isempty(Vpits)
        return;
    end

    sourcePos_nm = sourceData(:, 1:3) * 1e9;
    for pitIdx = 1:numel(Vpits)
        pit = Vpits{pitIdx};
        if isempty(pit) || ~isstruct(pit) || ~isfield(pit, 'center') || ~isfield(pit, 'topRadius') || ...
                ~isfield(pit, 'depth') || ~isfinite(pit.depth) || pit.depth <= 0
            continue;
        end
        surfaceZ_nm = nominalSurfaceZ;
        if isfield(pit, 'surfaceZ_nm') && isfinite(pit.surfaceZ_nm)
            surfaceZ_nm = pit.surfaceZ_nm;
        end
        radialDistance_nm = hypot(sourcePos_nm(:, 1) - pit.center(1), sourcePos_nm(:, 2) - pit.center(2));
        depthBelowSurface_nm = surfaceZ_nm - sourcePos_nm(:, 3);
        highBandgapDepthAllowance_nm = max(35, 155 - min(pit.depth, 155));
        nearMask = nearMask | (radialDistance_nm <= pit.topRadius + 35 & ...
            depthBelowSurface_nm >= -5 & depthBelowSurface_nm <= pit.depth + highBandgapDepthAllowance_nm);
    end
end

function transportDefaults = getCarrierTransportDefaults(params)
    if nargin < 1
        params = struct();
    end
    electronEnergy_keV = max(1, getFieldOrDefault(params, 'electronEnergy', 5e3) / 1e3);
    transportDefaults = struct( ...
        'electronEnergy_keV', electronEnergy_keV, ...
        'defaultLateralDiffusion_nm', 60, ...
        'defaultVerticalDiffusion_nm', 12, ...
        'defaultLateralDiffusion_m', 60e-9, ...
        'defaultVerticalDiffusion_m', 12e-9, ...
        'pathStep_nm', 2.5, ...
        'maxPathSteps', 120, ...
        'minWalkBeforeRecombine_nm', 4, ...
        'maxLateralStep_nm', 4.0, ...
        'maxVerticalStep_nm', 1.8, ...
        'driftStrength', 0.85, ...
        'facetCaptureDriftWeight', 0.60, ...
        'hotCarrierRetentionFactor', 0.11, ...
        'hotCarrierRelaxationLength_nm', 70, ...
        'maxEmissionHotShift_eV', 0.026, ...
        'minCarrierEnergy_eV', 0.01, ...
        'baseTransportEnergy_eV', 0.18, ...
        'mqwBonusEnergy_eV', 0.06, ...
        'semipolarBonusEnergy_eV', 0.16, ...
        'semipolarBarrierReduction', 0.44, ...
        'semipolarDirectRecombinationProbability', 0.48, ...
        'semipolarDirectCoreBonus', 0.20, ...
        'thermalAssist_eV', 0.030, ...
        'tunnelDeficitScale_eV', 0.040, ...
        'tunnelDecayLength_nm', 1.4, ...
        'postTunnelEnergyFraction', 0.65, ...
        'barrierDissipationFactor', 1.0, ...
        'lowVoltageShallowDirectScale', 1.0, ...
        'shallowDirectVoltagePivot_keV', 7.0, ...
        'shallowDirectVoltageWidth_keV', 3.0, ...
        'shallowDirectDepth_nm', 58, ...
        'shallowDirectMaxFraction', 0.82, ...
        'pTypeSurfaceDirectBonus', 0.28, ...
        'nonVpitShortwaveDirectFloor', 0.055, ...
        'vpitShortwaveDirectBoost', 2.65, ...
        'vpitFacetShortwaveWidth_nm', 18, ...
        'vpitShortwaveDepthCenterNorm', 0.42, ...
        'vpitShortwaveDepthWidthNorm', 0.34, ...
        'vpitCavitySidewallDirectScale', 0.74, ...
        'vpitCavitySidewallMaxFraction', 0.62, ...
        'vpitCavitySidewallDepthCenterNorm', 0.38, ...
        'vpitCavitySidewallDepthWidthNorm', 0.36, ...
        'vpitHighBandgapSidewallChannelShare', 0.74, ...
        'vpitHighBandgapSidewallChannelScale', 1.80, ...
        'vpitHighBandgapSidewallMaxFraction', 0.78, ...
        'vpitHighBandgapGuidedDepthBelowNTop_nm', 8, ...
        'vpitHighBandgapGuidedLateralJitter_nm', 4, ...
        'vpitHighBandgapNganAdjacencyMargin_nm', 28, ...
        'vpitInCompositionDepletionEnergyShift_eV', 0.030, ...
        'vpitQcseScreeningEnergyShift_eV', 0.020, ...
        'vpitStrainRelaxationEnergyShift_eV', 0.018, ...
        'vpitSidewallNonradiativeScale', 0.38, ...
        'enableVpitShortwaveEnhancement', isModelToggleEnabled(params, 'enableVpitTransportCorrection', true), ...
        'surfaceQuenchDepth_nm', 5, ...
        'surfaceQuenchStrength', 0.95, ...
        'vpitCenterQuenchWidthRatio', 0.18, ...
        'vpitCenterQuenchDepthNorm', 0.50, ...
        'vpitCenterQuenchDepthWidth', 0.24, ...
        'vpitCenterBaseQuench', 0.08, ...
        'vpitDefectBaseQuench', 0.045, ...
        'vpitDefectHaloCenterRatio', 0.78, ...
        'vpitDefectHaloWidthRatio', 0.18, ...
        'vpitDefectCoreWidthRatio', 0.12, ...
        'vpitDefectDepthNorm', 0.42, ...
        'vpitDefectDepthWidth', 0.24, ...
        'vpitCenterRefillBoost', 0.56, ...
        'vpitCenterRefillWidthRatio', 0.48, ...
        'vpitCenterRefillDepthNorm', 0.50, ...
        'vpitCenterRefillDepthWidth', 0.30, ...
        'vpitCenterRefillCollectionBonus', 0.32);

    structuralParams = getStructuralParameterStruct(params);
    diffusionLength_nm = getFieldOrDefault(structuralParams, 'carrierDiffusionLengthNm', NaN);
    if isfinite(diffusionLength_nm) && diffusionLength_nm > 0
        lateralScale = diffusionLength_nm / max(transportDefaults.defaultLateralDiffusion_nm, eps);
        verticalScale = sqrt(max(lateralScale, 0.05));
        transportDefaults.defaultLateralDiffusion_nm = diffusionLength_nm;
        transportDefaults.defaultLateralDiffusion_m = diffusionLength_nm * 1e-9;
        transportDefaults.defaultVerticalDiffusion_nm = clamp(transportDefaults.defaultVerticalDiffusion_nm * verticalScale, 3, 40);
        transportDefaults.defaultVerticalDiffusion_m = transportDefaults.defaultVerticalDiffusion_nm * 1e-9;
    end

    quenchScale = getFieldOrDefault(structuralParams, 'vpitQuenchingScale', ...
        getFieldOrDefault(structuralParams, 'vpitQuenchStrength', 1));
    if isfinite(quenchScale)
        quenchScale = clamp(quenchScale, 0, 4);
        transportDefaults.vpitCenterBaseQuench = clamp(transportDefaults.vpitCenterBaseQuench * quenchScale, 0, 0.95);
        transportDefaults.vpitDefectBaseQuench = clamp(transportDefaults.vpitDefectBaseQuench * quenchScale, 0, 0.95);
    end

    semipolarScale = getFieldOrDefault(structuralParams, 'semipolarCaptureScale', 1);
    if isfinite(semipolarScale)
        semipolarScale = clamp(semipolarScale, 0.1, 3);
        transportDefaults.semipolarDirectRecombinationProbability = clamp( ...
            transportDefaults.semipolarDirectRecombinationProbability * semipolarScale, 0.02, 0.95);
        transportDefaults.facetCaptureDriftWeight = clamp( ...
            transportDefaults.facetCaptureDriftWeight * sqrt(semipolarScale), 0.05, 1.5);
    end

    refillScale = getFieldOrDefault(structuralParams, 'vpitCenterRefillScale', 1);
    if isfinite(refillScale)
        refillScale = clamp(refillScale, 0, 3);
        transportDefaults.vpitCenterRefillBoost = clamp(transportDefaults.vpitCenterRefillBoost * refillScale, 0, 2.5);
        transportDefaults.vpitCenterRefillCollectionBonus = clamp( ...
            transportDefaults.vpitCenterRefillCollectionBonus * refillScale, 0, 2.0);
    end

    directScale = getFieldOrDefault(structuralParams, 'lowVoltageShallowDirectScale', ...
        getFieldOrDefault(structuralParams, 'shortwaveDirectRecombinationScale', 1));
    if isfinite(directScale)
        transportDefaults.lowVoltageShallowDirectScale = clamp(directScale, 0, 3);
    end

    vpitShortwaveScale = getFieldOrDefault(structuralParams, 'vpitShortwaveChannelScale', 1);
    if isfinite(vpitShortwaveScale)
        transportDefaults.vpitShortwaveDirectBoost = clamp( ...
            transportDefaults.vpitShortwaveDirectBoost * vpitShortwaveScale, 0, 5);
        transportDefaults.vpitCavitySidewallDirectScale = clamp( ...
            transportDefaults.vpitCavitySidewallDirectScale * sqrt(max(vpitShortwaveScale, 0)), 0, 2.5);
        transportDefaults.vpitHighBandgapSidewallChannelScale = clamp( ...
            transportDefaults.vpitHighBandgapSidewallChannelScale * vpitShortwaveScale, 0, 6);
    end

    inDepletionScale = getFieldOrDefault(structuralParams, 'vpitInCompositionGradientScale', ...
        getFieldOrDefault(structuralParams, 'vpitInDepletionScale', 1));
    if isfinite(inDepletionScale)
        transportDefaults.vpitInCompositionDepletionEnergyShift_eV = clamp( ...
            transportDefaults.vpitInCompositionDepletionEnergyShift_eV * inDepletionScale, 0, 0.11);
    end

    qcseScale = getFieldOrDefault(structuralParams, 'vpitQcseScreeningScale', 1);
    if isfinite(qcseScale)
        transportDefaults.vpitQcseScreeningEnergyShift_eV = clamp( ...
            transportDefaults.vpitQcseScreeningEnergyShift_eV * qcseScale, 0, 0.08);
    end

    strainScale = getFieldOrDefault(structuralParams, 'vpitStrainRelaxationScale', 1);
    if isfinite(strainScale)
        transportDefaults.vpitStrainRelaxationEnergyShift_eV = clamp( ...
            transportDefaults.vpitStrainRelaxationEnergyShift_eV * strainScale, 0, 0.08);
    end

    sidewallNonradiativeScale = getFieldOrDefault(structuralParams, 'vpitSidewallNonradiativeScale', ...
        getFieldOrDefault(structuralParams, 'dopingNonradiativeScale', 1));
    if isfinite(sidewallNonradiativeScale)
        transportDefaults.vpitSidewallNonradiativeScale = clamp(sidewallNonradiativeScale, 0, 3);
    end
end

function directFraction = resolveShallowDirectRecombinationFraction(position_nm, state, nominalSurfaceZ_nm, transportDefaults, Vpits)
    directFraction = 0;
    if ~isCarrierStateAccessible(state)
        return;
    end

    sourceTypeCode = classifyShallowDirectEmissionSource(state);
    if ~ismember(sourceTypeCode, [1, 2, 3, 6])
        return;
    end

    depthBelowSurface_nm = max(0, nominalSurfaceZ_nm - position_nm(3));
    lowVoltageWeight = clamp( ...
        (transportDefaults.shallowDirectVoltagePivot_keV - transportDefaults.electronEnergy_keV) / ...
        max(transportDefaults.shallowDirectVoltageWidth_keV, eps), 0, 1);
    depthWeight = exp(-depthBelowSurface_nm / max(transportDefaults.shallowDirectDepth_nm, eps));

    baseByType = [0.52, 0.12, 0.16, 0.0, 0.0, 0.88];
    typeBase = baseByType(sourceTypeCode);
    bandgapWeight = 0.80;
    if isfield(state, 'bandgap_eV') && isfinite(state.bandgap_eV)
        bandgapWeight = 0.65 + 0.35 * clamp((state.bandgap_eV - 2.65) / 0.65, 0, 1);
    end
    vpitCoupling = resolveVpitShortwaveChannelCoupling(position_nm, state, transportDefaults, Vpits);

    directFraction = transportDefaults.lowVoltageShallowDirectScale * typeBase * ...
        (0.18 + 0.82 * lowVoltageWeight) * (0.28 + 0.72 * depthWeight) * ...
        bandgapWeight * vpitCoupling;

    layerName = string(getFieldOrDefault(state, 'layerName', ""));
    if sourceTypeCode == 6 && contains(layerName, 'p-GaN')
        pSurfaceWeight = exp(-depthBelowSurface_nm / 18);
        directFraction = directFraction + ...
            transportDefaults.lowVoltageShallowDirectScale * transportDefaults.pTypeSurfaceDirectBonus * ...
            (0.25 + 0.75 * lowVoltageWeight) * pSurfaceWeight * vpitCoupling;
    end

    directFraction = clamp(directFraction, 0, transportDefaults.shallowDirectMaxFraction);
end

function sourceTypeCode = classifyShallowDirectEmissionSource(state)
    sourceTypeCode = 0;
    layerName = string(getFieldOrDefault(state, 'layerName', ""));
    if contains(layerName, 'p-EBL') || contains(layerName, 'p-GaN')
        sourceTypeCode = 6;
    elseif contains(layerName, 'MQW-Barrier')
        sourceTypeCode = 3;
    elseif contains(layerName, 'Prestrained')
        sourceTypeCode = 2;
    elseif contains(layerName, 'GaN')
        sourceTypeCode = 1;
    end
end

function coupling = resolveVpitShortwaveChannelCoupling(position_nm, state, transportDefaults, Vpits)
    coupling = transportDefaults.nonVpitShortwaveDirectFloor;
    if ~transportDefaults.enableVpitShortwaveEnhancement || isempty(Vpits) || ...
            ~isfield(state, 'vpitIndex') || state.vpitIndex < 1 || state.vpitIndex > numel(Vpits)
        return;
    end

    vpit = Vpits{state.vpitIndex};
    if isempty(vpit) || ~isstruct(vpit) || ~isfield(vpit, 'depth') || ~isfield(vpit, 'surfaceZ_nm') || ...
            ~isfinite(vpit.depth) || vpit.depth <= 0
        return;
    end

    depthBelowSurface_nm = max(0, vpit.surfaceZ_nm - position_nm(3));
    if depthBelowSurface_nm > vpit.depth + 2
        return;
    end

    facetDistance_nm = getFieldOrDefault(state, 'closestFacetDistance_nm', inf);
    shellThickness_nm = getFieldOrDefault(state, 'layerShellThickness_nm', 0);
    if strcmp(getFieldOrDefault(state, 'materialKind', ''), 'semipolar_shell')
        facetWeight = 1;
    elseif isfinite(facetDistance_nm)
        facetWidth_nm = max(transportDefaults.vpitFacetShortwaveWidth_nm, 1.5 * max(shellThickness_nm, 1));
        facetWeight = exp(-max(facetDistance_nm, 0) / facetWidth_nm);
    else
        facetWeight = 0.25;
    end

    depthNorm = clamp(depthBelowSurface_nm / max(vpit.depth, eps), 0, 1);
    depthKernel = 0.35 + 0.65 * exp(-0.5 * ...
        ((depthNorm - transportDefaults.vpitShortwaveDepthCenterNorm) / ...
        max(transportDefaults.vpitShortwaveDepthWidthNorm, 0.05)) ^ 2);

    apertureBoost = 1.0;
    currentApothem_nm = currentVpitApothemAtZ(vpit, position_nm(3));
    if isfinite(currentApothem_nm) && currentApothem_nm > 1
        radialDistance_nm = norm(position_nm(1:2) - vpit.center(:)');
        radialRatio = radialDistance_nm / max(currentApothem_nm, 1);
        apertureBoost = 0.72 + 0.28 * exp(-0.5 * ((radialRatio - 0.72) / 0.30) ^ 2);
    end

    coupling = transportDefaults.nonVpitShortwaveDirectFloor + ...
        transportDefaults.vpitShortwaveDirectBoost * (0.40 + 0.60 * facetWeight) * depthKernel * apertureBoost;
    coupling = clamp(coupling, transportDefaults.nonVpitShortwaveDirectFloor, 3.20);
end

function probability = resolveSemipolarDirectRecombinationProbability(state, transportDefaults)
    probability = transportDefaults.semipolarDirectRecombinationProbability;
    if isfield(state, 'closestFacetDistance_nm') && isfinite(state.closestFacetDistance_nm) && ...
            isfield(state, 'layerShellThickness_nm') && isfinite(state.layerShellThickness_nm) && state.layerShellThickness_nm > 0
        coreProximity = 1 - clamp(state.closestFacetDistance_nm / max(state.layerShellThickness_nm, eps), 0, 1);
        probability = probability + transportDefaults.semipolarDirectCoreBonus * coreProximity;
    end
    probability = clamp(probability, 0.18, 0.78);
end

function radiativeYield = computeCarrierRadiativeYield(originPos_nm, originState, finalPos_nm, finalState, Vpits, transportDefaults)
    radiativeYield = 1;
    if ~isCarrierStateAccessible(finalState)
        radiativeYield = 0;
        return;
    end

    radiativeYield = radiativeYield * computeVpitCenterRadiativeRetention( ...
        originPos_nm, originState, finalPos_nm, finalState, Vpits, transportDefaults);
    radiativeYield = radiativeYield * computeVpitDefectFieldRetention( ...
        originPos_nm, originState, finalPos_nm, finalState, Vpits, transportDefaults);
    radiativeYield = radiativeYield * computeVpitCenterRefillBoost( ...
        originPos_nm, originState, finalPos_nm, finalState, Vpits, transportDefaults);
    radiativeYield = max(0.04, min(1.0, radiativeYield));
end

function retentionFactor = computeVpitCenterRadiativeRetention(originPos_nm, originState, finalPos_nm, finalState, Vpits, transportDefaults)
    retentionFactor = 1;
    [referencePos_nm, referenceState] = selectVpitReferenceState(originPos_nm, originState, finalPos_nm, finalState);
    if ~isfield(referenceState, 'vpitIndex') || referenceState.vpitIndex < 1 || referenceState.vpitIndex > numel(Vpits)
        return;
    end

    vpit = Vpits{referenceState.vpitIndex};
    if isempty(vpit) || ~isstruct(vpit) || ~isfield(vpit, 'surfaceZ_nm') || ~isfield(vpit, 'depth') || ...
            ~isfinite(vpit.surfaceZ_nm) || ~isfinite(vpit.depth) || vpit.depth <= 0
        return;
    end

    depthBelowSurface_nm = vpit.surfaceZ_nm - referencePos_nm(3);
    if depthBelowSurface_nm < 0 || depthBelowSurface_nm > vpit.depth
        return;
    end

    currentApothem_nm = currentVpitApothemAtZ(vpit, referencePos_nm(3));
    radialDistance_nm = norm(referencePos_nm(1:2) - vpit.center(:)');
    radialRatio = radialDistance_nm / max(currentApothem_nm, 5);
    centerKernel = exp(-0.5 * (radialRatio / transportDefaults.vpitCenterQuenchWidthRatio) ^ 2);

    depthNorm = clamp(depthBelowSurface_nm / max(vpit.depth, eps), 0, 1);
    depthKernel = 0.74 + 0.26 * exp(-0.5 * ...
        ((depthNorm - transportDefaults.vpitCenterQuenchDepthNorm) / transportDefaults.vpitCenterQuenchDepthWidth) ^ 2);

    materialSensitivity = getCarrierRadiativeQuenchSensitivity(originState, finalState);
    quenchStrength = transportDefaults.vpitCenterBaseQuench * materialSensitivity * centerKernel * depthKernel;
    retentionFactor = 1 - clamp(quenchStrength, 0, 0.92);
end

function retentionFactor = computeVpitDefectFieldRetention(originPos_nm, originState, finalPos_nm, finalState, Vpits, transportDefaults)
    retentionFactor = 1;
    [referencePos_nm, referenceState] = selectVpitReferenceState(originPos_nm, originState, finalPos_nm, finalState);
    if ~isfield(referenceState, 'vpitIndex') || referenceState.vpitIndex < 1 || referenceState.vpitIndex > numel(Vpits)
        return;
    end

    vpit = Vpits{referenceState.vpitIndex};
    if isempty(vpit) || ~isstruct(vpit) || ~isfield(vpit, 'surfaceZ_nm') || ~isfield(vpit, 'depth') || ...
            ~isfinite(vpit.surfaceZ_nm) || ~isfinite(vpit.depth) || vpit.depth <= 0
        return;
    end

    depthBelowSurface_nm = vpit.surfaceZ_nm - referencePos_nm(3);
    if depthBelowSurface_nm < 0 || depthBelowSurface_nm > vpit.depth
        return;
    end

    currentApothem_nm = currentVpitApothemAtZ(vpit, referencePos_nm(3));
    radialDistance_nm = norm(referencePos_nm(1:2) - vpit.center(:)');
    radialRatio = radialDistance_nm / max(currentApothem_nm, 5);
    depthNorm = clamp(depthBelowSurface_nm / max(vpit.depth, eps), 0, 1);

    haloKernel = exp(-0.5 * ...
        ((radialRatio - transportDefaults.vpitDefectHaloCenterRatio) / max(transportDefaults.vpitDefectHaloWidthRatio, 0.05)) ^ 2);
    coreKernel = exp(-0.5 * ...
        (radialRatio / max(transportDefaults.vpitDefectCoreWidthRatio, 0.05)) ^ 2);
    depthKernel = 0.60 + 0.40 * exp(-0.5 * ...
        ((depthNorm - transportDefaults.vpitDefectDepthNorm) / max(transportDefaults.vpitDefectDepthWidth, 0.05)) ^ 2);
    defectKernel = 0.86 * haloKernel + 0.14 * coreKernel;

    if strcmp(referenceState.materialKind, 'semipolar_shell') || strcmp(finalState.materialKind, 'semipolar_shell')
        defectKernel = defectKernel * 0.82;
    end

    defectSensitivity = getCarrierDefectFieldSensitivity(originState, finalState);
    quenchStrength = transportDefaults.vpitDefectBaseQuench * defectSensitivity * defectKernel * depthKernel;
    retentionFactor = 1 - clamp(quenchStrength, 0, 0.42);
end

function boostFactor = computeVpitCenterRefillBoost(originPos_nm, originState, finalPos_nm, finalState, Vpits, transportDefaults)
    boostFactor = 1;
    [~, referenceState] = selectVpitReferenceState(originPos_nm, originState, finalPos_nm, finalState);
    if ~isfield(referenceState, 'vpitIndex') || referenceState.vpitIndex < 1 || referenceState.vpitIndex > numel(Vpits)
        return;
    end

    vpit = Vpits{referenceState.vpitIndex};
    if isempty(vpit) || ~isstruct(vpit) || ~isfield(vpit, 'surfaceZ_nm') || ~isfield(vpit, 'depth') || ...
            ~isfinite(vpit.surfaceZ_nm) || ~isfinite(vpit.depth) || vpit.depth <= 0
        return;
    end

    depthBelowSurface_nm = vpit.surfaceZ_nm - finalPos_nm(3);
    if depthBelowSurface_nm < 0 || depthBelowSurface_nm > vpit.depth
        return;
    end

    currentApothem_nm = currentVpitApothemAtZ(vpit, finalPos_nm(3));
    radialDistance_nm = norm(finalPos_nm(1:2) - vpit.center(:)');
    radialRatio = radialDistance_nm / max(currentApothem_nm, 5);
    depthNorm = clamp(depthBelowSurface_nm / max(vpit.depth, eps), 0, 1);

    centerKernel = exp(-0.5 * (radialRatio / max(transportDefaults.vpitCenterRefillWidthRatio, 0.08)) ^ 2);
    depthKernel = exp(-0.5 * ...
        ((depthNorm - transportDefaults.vpitCenterRefillDepthNorm) / max(transportDefaults.vpitCenterRefillDepthWidth, 0.08)) ^ 2);

    emissionType = classifyCarrierEmissionSource(finalState);
    switch emissionType
        case 4
            materialBoost = 1.00;
        case 5
            materialBoost = 1.05;
        case 3
            materialBoost = 0.62;
        otherwise
            materialBoost = 0.28;
    end

    collectionBonus = 0;
    if strcmp(getFieldOrDefault(finalState, 'materialKind', ''), 'bulk') && getFieldOrDefault(finalState, 'isMQW', false)
        collectionBonus = transportDefaults.vpitCenterRefillCollectionBonus;
    elseif strcmp(getFieldOrDefault(originState, 'materialKind', ''), 'semipolar_shell')
        collectionBonus = 0.6 * transportDefaults.vpitCenterRefillCollectionBonus;
    end

    boostStrength = (transportDefaults.vpitCenterRefillBoost + collectionBonus) * materialBoost;
    boostFactor = 1 + boostStrength * centerKernel * depthKernel;
    boostFactor = min(boostFactor, 1.52);
end

function [referencePos_nm, referenceState] = selectVpitReferenceState(originPos_nm, originState, finalPos_nm, finalState)
    referencePos_nm = finalPos_nm;
    referenceState = finalState;
    if isfield(originState, 'vpitIndex') && originState.vpitIndex > 0
        referencePos_nm = originPos_nm;
        referenceState = originState;
    end
end

function currentApothem_nm = currentVpitApothemAtZ(vpit, z_nm)
    if ~isfield(vpit, 'topRadius') || ~isfield(vpit, 'surfaceZ_nm') || ~isfield(vpit, 'depth') || ...
            ~isfinite(vpit.topRadius) || ~isfinite(vpit.surfaceZ_nm) || ~isfinite(vpit.depth) || vpit.depth <= 0
        currentApothem_nm = 5;
        return;
    end

    depthBelowSurface_nm = vpit.surfaceZ_nm - z_nm;
    currentApothem_nm = vpit.topRadius * max(0.06, 1 - depthBelowSurface_nm / max(vpit.depth, eps));
    currentApothem_nm = max(currentApothem_nm, 5);
end

function materialSensitivity = getCarrierRadiativeQuenchSensitivity(originState, finalState)
    originType = classifyCarrierEmissionSource(originState);
    finalType = classifyCarrierEmissionSource(finalState);
    effectiveType = originType;
    if effectiveType < 1 || effectiveType > 6
        effectiveType = finalType;
    end

    switch effectiveType
        case 1
            materialSensitivity = 0.68;
        case 2
            materialSensitivity = 0.74;
        case 3
            materialSensitivity = 0.84;
        case 4
            materialSensitivity = 0.48;
        case 5
            materialSensitivity = 0.62;
        case 6
            materialSensitivity = 0.72;
        otherwise
            materialSensitivity = 0.70;
    end

    if originType == 5 && finalType == 4
        materialSensitivity = max(materialSensitivity, 0.68);
    end
end

function defectSensitivity = getCarrierDefectFieldSensitivity(originState, finalState)
    originType = classifyCarrierEmissionSource(originState);
    finalType = classifyCarrierEmissionSource(finalState);
    effectiveType = originType;
    if effectiveType < 1 || effectiveType > 6
        effectiveType = finalType;
    end

    switch effectiveType
        case 1
            defectSensitivity = 0.80;
        case 2
            defectSensitivity = 0.88;
        case 3
            defectSensitivity = 0.82;
        case 4
            defectSensitivity = 0.42;
        case 5
            defectSensitivity = 0.42;
        case 6
            defectSensitivity = 0.72;
        otherwise
            defectSensitivity = 0.72;
    end

    if originType == 5 && finalType == 4
        defectSensitivity = min(defectSensitivity, 0.44);
    end
end

function pathResult = propagateCarrierDiffusionPath(startPos_nm, lateralSigma_nm, verticalSigma_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ, Eg_GaN, Eg_InN, transportDefaults)
    currentPos_nm = startPos_nm;
    currentState = precies.vpit('resolvePointState', currentPos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);

    pathResult = struct( ...
        'isValid', false, ...
        'finalPos_nm', currentPos_nm, ...
        'finalState', currentState, ...
        'remainingEnergy_eV', 0, ...
        'pathLength_nm', 0);

    if ~isCarrierStateAccessible(currentState)
        return;
    end

    remainingEnergy_eV = initializeCarrierTransportEnergy(currentState, transportDefaults);
    pathLength_nm = 0;
    topBoundary_nm = layerBoundaries(end - 1) * 1e9;
    bottomBoundary_nm = layerBoundaries(1) * 1e9;
    numSteps = max(10, min(transportDefaults.maxPathSteps, ...
        ceil(4 * max(lateralSigma_nm, verticalSigma_nm) / transportDefaults.pathStep_nm)));

    for stepIdx = 1:numSteps
        if (topBoundary_nm - currentPos_nm(3)) < transportDefaults.surfaceQuenchDepth_nm
            surfaceQuenchProb = 1 - exp(-transportDefaults.surfaceQuenchStrength * transportDefaults.pathStep_nm / ...
                transportDefaults.surfaceQuenchDepth_nm);
            if rand() < surfaceQuenchProb
                return;
            end
        end

        recombinationLength_nm = getCarrierRecombinationLength(currentState, lateralSigma_nm, verticalSigma_nm, transportDefaults);
        if pathLength_nm >= transportDefaults.minWalkBeforeRecombine_nm
            recombinationProb = 1 - exp(-transportDefaults.pathStep_nm / max(recombinationLength_nm, eps));
            if rand() < recombinationProb
                pathResult.isValid = true;
                pathResult.finalPos_nm = currentPos_nm;
                pathResult.finalState = currentState;
                pathResult.remainingEnergy_eV = max(remainingEnergy_eV, 0);
                pathResult.pathLength_nm = pathLength_nm;
                return;
            end
        end

        [stepVector_nm, stepLength_nm] = proposeCarrierStep( ...
            currentPos_nm, currentState, layers, layerBoundaries, Vpits, nominalSurfaceZ, ...
            Eg_GaN, Eg_InN, lateralSigma_nm, verticalSigma_nm, transportDefaults);

        candidateAccepted = false;
        for attemptIdx = 1:4
            nextPos_nm = currentPos_nm + stepVector_nm;
            if nextPos_nm(3) < bottomBoundary_nm || nextPos_nm(3) > topBoundary_nm
                [stepVector_nm, stepLength_nm] = proposeCarrierStep( ...
                    currentPos_nm, currentState, layers, layerBoundaries, Vpits, nominalSurfaceZ, ...
                    Eg_GaN, Eg_InN, lateralSigma_nm, verticalSigma_nm, transportDefaults);
                continue;
            end

            nextState = precies.vpit('resolvePointState', nextPos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
            if isCarrierStateAccessible(nextState)
                candidateAccepted = true;
                break;
            end

            [stepVector_nm, stepLength_nm] = proposeCarrierStep( ...
                currentPos_nm, currentState, layers, layerBoundaries, Vpits, nominalSurfaceZ, ...
                Eg_GaN, Eg_InN, lateralSigma_nm, verticalSigma_nm, transportDefaults);
        end

        if ~candidateAccepted
            pathResult.isValid = true;
            pathResult.finalPos_nm = currentPos_nm;
            pathResult.finalState = currentState;
            pathResult.remainingEnergy_eV = max(remainingEnergy_eV, 0);
            pathResult.pathLength_nm = pathLength_nm;
            return;
        end

        attenuationLength_nm = getCarrierAttenuationLength(currentState, transportDefaults);
        remainingEnergy_eV = remainingEnergy_eV * exp(-stepLength_nm / max(attenuationLength_nm, eps));

        barrierHeight_eV = computeCarrierBarrierHeight(currentState, nextState, layers, Eg_GaN, Eg_InN, transportDefaults);
        if barrierHeight_eV > 0
            [canPassBarrier, remainingEnergy_eV] = attemptCarrierBarrierCrossing( ...
                remainingEnergy_eV, barrierHeight_eV, stepLength_nm, nextState, transportDefaults);
            if ~canPassBarrier
                pathResult.isValid = true;
                pathResult.finalPos_nm = currentPos_nm;
                pathResult.finalState = currentState;
                pathResult.remainingEnergy_eV = max(remainingEnergy_eV, 0);
                pathResult.pathLength_nm = pathLength_nm;
                return;
            end
        end

        currentPos_nm = nextPos_nm;
        currentState = nextState;
        pathLength_nm = pathLength_nm + stepLength_nm;

        if remainingEnergy_eV <= transportDefaults.minCarrierEnergy_eV
            pathResult.isValid = true;
            pathResult.finalPos_nm = currentPos_nm;
            pathResult.finalState = currentState;
            pathResult.remainingEnergy_eV = max(remainingEnergy_eV, 0);
            pathResult.pathLength_nm = pathLength_nm;
            return;
        end
    end

    pathResult.isValid = true;
    pathResult.finalPos_nm = currentPos_nm;
    pathResult.finalState = currentState;
    pathResult.remainingEnergy_eV = max(remainingEnergy_eV, 0);
    pathResult.pathLength_nm = pathLength_nm;
end

function energy_eV = initializeCarrierTransportEnergy(state, transportDefaults)
    energy_eV = transportDefaults.baseTransportEnergy_eV;
    if state.isMQW
        energy_eV = energy_eV + transportDefaults.mqwBonusEnergy_eV;
    end
    if strcmp(state.materialKind, 'semipolar_shell')
        energy_eV = energy_eV + transportDefaults.semipolarBonusEnergy_eV;
    end
end

function attenuationLength_nm = getCarrierAttenuationLength(state, ~)
    if strcmp(state.materialKind, 'semipolar_shell')
        attenuationLength_nm = 220;
        return;
    end

    layerName = '';
    if isfield(state, 'layerName') && ~isempty(state.layerName)
        layerName = state.layerName;
    end

    if contains(layerName, 'MQW-Well')
        attenuationLength_nm = 95;
    elseif contains(layerName, 'MQW-Barrier')
        attenuationLength_nm = 35;
    elseif contains(layerName, 'Prestrained') || contains(layerName, 'p-EBL')
        attenuationLength_nm = 20;
    elseif contains(layerName, 'GaN')
        attenuationLength_nm = 130;
    else
        attenuationLength_nm = 80;
    end
end

function recombinationLength_nm = getCarrierRecombinationLength(state, lateralSigma_nm, verticalSigma_nm, ~)
    if strcmp(state.materialKind, 'semipolar_shell')
        recombinationLength_nm = max(180, 2.4 * lateralSigma_nm);
        return;
    end

    layerName = '';
    if isfield(state, 'layerName') && ~isempty(state.layerName)
        layerName = state.layerName;
    end

    if contains(layerName, 'MQW-Well')
        recombinationLength_nm = max(22, 0.55 * lateralSigma_nm);
    elseif contains(layerName, 'MQW-Barrier')
        recombinationLength_nm = max(45, 1.35 * lateralSigma_nm);
    elseif contains(layerName, 'Prestrained') || contains(layerName, 'p-EBL')
        recombinationLength_nm = max(14, 0.70 * verticalSigma_nm);
    elseif contains(layerName, 'GaN')
        recombinationLength_nm = max(110, 2.1 * lateralSigma_nm);
    else
        recombinationLength_nm = max(70, 1.5 * lateralSigma_nm);
    end
end

function [stepVector_nm, stepLength_nm] = proposeCarrierStep(currentPos_nm, currentState, layers, layerBoundaries, Vpits, nominalSurfaceZ, Eg_GaN, Eg_InN, lateralSigma_nm, verticalSigma_nm, transportDefaults)
    lateralStep_nm = min(transportDefaults.maxLateralStep_nm, max(0.8, lateralSigma_nm / 20));
    verticalStep_nm = min(transportDefaults.maxVerticalStep_nm, max(0.4, verticalSigma_nm / 10));

    if strcmp(currentState.materialKind, 'semipolar_shell')
        lateralStep_nm = lateralStep_nm * 1.45;
        verticalStep_nm = verticalStep_nm * 0.75;
    end

    randomStep_nm = [ ...
        lateralStep_nm * randn(), ...
        lateralStep_nm * randn(), ...
        verticalStep_nm * randn()];

    driftDir = computeCarrierDriftDirection( ...
        currentPos_nm, currentState, layers, layerBoundaries, Vpits, nominalSurfaceZ, ...
        Eg_GaN, Eg_InN, transportDefaults);
    driftStep_nm = transportDefaults.driftStrength * transportDefaults.pathStep_nm * driftDir;

    stepVector_nm = randomStep_nm + driftStep_nm;
    if norm(stepVector_nm) <= eps
        stepVector_nm = [lateralStep_nm, 0, 0];
    end
    stepLength_nm = norm(stepVector_nm);
end

function driftDir = computeCarrierDriftDirection(currentPos_nm, currentState, layers, layerBoundaries, Vpits, nominalSurfaceZ, Eg_GaN, Eg_InN, transportDefaults)
    driftDir = [0, 0, 0];
    currentBandgap_eV = getCarrierStateBandgap(currentState, layers, Eg_GaN, Eg_InN);

    if currentState.layerIndex > 1
        lowerLayerState = currentState;
        lowerLayerState.layerIndex = currentState.layerIndex - 1;
        lowerLayerState.bandgap_eV = get_bandgap(lowerLayerState.layerIndex, layers, Eg_GaN, Eg_InN);
        lowerDelta = max(0, currentBandgap_eV - lowerLayerState.bandgap_eV);
        driftDir(3) = driftDir(3) - lowerDelta;
    end
    if currentState.layerIndex < size(layers, 1)
        upperLayerState = currentState;
        upperLayerState.layerIndex = currentState.layerIndex + 1;
        upperLayerState.bandgap_eV = get_bandgap(upperLayerState.layerIndex, layers, Eg_GaN, Eg_InN);
        upperDelta = max(0, currentBandgap_eV - upperLayerState.bandgap_eV);
        driftDir(3) = driftDir(3) + upperDelta;
    end

    if currentState.isMQW && isfinite(currentState.closestFacetDistance_nm) && ...
            currentState.closestFacetDistance_nm <= transportDefaults.defaultLateralDiffusion_nm
        shellWeight = 1 - currentState.closestFacetDistance_nm / transportDefaults.defaultLateralDiffusion_nm;
        driftDir = driftDir - transportDefaults.facetCaptureDriftWeight * shellWeight * currentState.closestFacetNormal3D;
    elseif strcmp(currentState.materialKind, 'semipolar_shell')
        shellStateOutside = precies.vpit('resolvePointState', currentPos_nm + 6 * currentState.closestFacetNormal3D, ...
            layers, layerBoundaries, Vpits, nominalSurfaceZ);
        if isCarrierStateAccessible(shellStateOutside)
            driftDir = driftDir + 0.35 * currentState.closestFacetNormal3D;
        end
    end

    if norm(driftDir) > eps
        driftDir = driftDir / norm(driftDir);
    end
end

function barrierHeight_eV = computeCarrierBarrierHeight(currentState, nextState, layers, Eg_GaN, Eg_InN, transportDefaults)
    currentBandgap_eV = getCarrierStateBandgap(currentState, layers, Eg_GaN, Eg_InN);
    nextBandgap_eV = getCarrierStateBandgap(nextState, layers, Eg_GaN, Eg_InN);
    if ~isfinite(currentBandgap_eV) || ~isfinite(nextBandgap_eV)
        barrierHeight_eV = 0;
        return;
    end

    barrierHeight_eV = max(0, nextBandgap_eV - currentBandgap_eV);
    if strcmp(currentState.materialKind, 'semipolar_shell') || strcmp(nextState.materialKind, 'semipolar_shell')
        barrierHeight_eV = barrierHeight_eV * transportDefaults.semipolarBarrierReduction;
    end
end

function bandgap_eV = getCarrierStateBandgap(state, layers, Eg_GaN, Eg_InN)
    if isfield(state, 'bandgap_eV') && isfinite(state.bandgap_eV)
        bandgap_eV = state.bandgap_eV;
    elseif isfield(state, 'layerIndex') && state.layerIndex >= 1 && state.layerIndex <= size(layers, 1)
        bandgap_eV = get_bandgap(state.layerIndex, layers, Eg_GaN, Eg_InN);
    else
        bandgap_eV = NaN;
    end
end

function [canPassBarrier, remainingEnergy_eV] = attemptCarrierBarrierCrossing(remainingEnergy_eV, barrierHeight_eV, segmentLength_nm, nextState, transportDefaults)
    if remainingEnergy_eV >= barrierHeight_eV
        remainingEnergy_eV = remainingEnergy_eV - barrierHeight_eV * transportDefaults.barrierDissipationFactor;
        canPassBarrier = true;
        return;
    end

    energyDeficit_eV = barrierHeight_eV - remainingEnergy_eV;
    barrierWidth_nm = max(segmentLength_nm, transportDefaults.pathStep_nm);
    thermalProb = exp(-energyDeficit_eV / transportDefaults.thermalAssist_eV);
    tunnelProb = exp(-energyDeficit_eV / transportDefaults.tunnelDeficitScale_eV) * ...
        exp(-barrierWidth_nm / transportDefaults.tunnelDecayLength_nm);

    if contains(nextState.layerName, 'p-EBL') || contains(nextState.layerName, 'Prestrained')
        thermalProb = thermalProb * 0.15;
        tunnelProb = tunnelProb * 0.35;
    end

    penetrationProb = min(1, thermalProb + tunnelProb);
    canPassBarrier = rand() <= penetrationProb;
    if canPassBarrier
        remainingEnergy_eV = max( ...
            transportDefaults.minCarrierEnergy_eV, ...
            remainingEnergy_eV * transportDefaults.postTunnelEnergyFraction);
    end
end

function hotShift_eV = computeCarrierEmissionHotShift(remainingEnergy_eV, pathLength_nm, finalState, transportDefaults)
    retentionFactor = transportDefaults.hotCarrierRetentionFactor * ...
        exp(-pathLength_nm / transportDefaults.hotCarrierRelaxationLength_nm);

    if strcmp(finalState.materialKind, 'semipolar_shell')
        retentionFactor = retentionFactor * 0.7;
    elseif isfield(finalState, 'isMQW') && finalState.isMQW
        retentionFactor = retentionFactor * 1.0;
    else
        retentionFactor = retentionFactor * 0.5;
    end

    hotShift_eV = min(transportDefaults.maxEmissionHotShift_eV, ...
        max(0, remainingEnergy_eV * retentionFactor));
end

function tf = isCarrierStateAccessible(state)
    tf = ~(strcmp(state.materialKind, 'cavity') || strcmp(state.materialKind, 'air'));
end

function sourceTypeCode = classifyCarrierEmissionSource(state)
    sourceTypeCode = 0;
    if strcmp(state.materialKind, 'semipolar_shell')
        sourceTypeCode = 5;
        return;
    end

    layerName = string(state.layerName);
    if contains(layerName, 'p-EBL') || contains(layerName, 'p-GaN')
        sourceTypeCode = 6;
    elseif contains(layerName, 'MQW-Well')
        sourceTypeCode = 4;
    elseif contains(layerName, 'MQW-Barrier')
        sourceTypeCode = 3;
    elseif contains(layerName, 'Prestrained')
        sourceTypeCode = 2;
    elseif contains(layerName, 'GaN')
        sourceTypeCode = 1;
    end
end

function theta = polarScatteringAngleVectorized(E)
    E_kev = E ./ 1e3;
    a1 = -3.73265; b1 = 0.13947; c1 = -0.0793875; d1 = 0.995866;
    a2 = 0.655582; b2 = 0.00382064; c2 = -0.00809; d2 = 0.00776978;
    logE = log10(E_kev);
    log_alpha = a1 + b1 .* logE + c1 .* (logE .^ 2) + d1 ./ exp(logE);
    alpha = 10 .^ log_alpha;
    beta = a2 + b2 .* sqrt(E_kev) .* log(E_kev) + c2 .* log(E_kev) ./ E_kev + d2 ./ E_kev;
    beta = min(beta, 1);
    R = randomLike(E_kev, size(E_kev));
    cos_theta_beta = 1 - (2 .* alpha .* R) ./ (1 + alpha - R);
    theta_beta = acos(max(-1, min(cos_theta_beta, 1)));
    theta = theta_beta .^ (1 ./ beta);
end

function new_dir = rotateVectorBatch(dir, theta, phi)
    if isempty(dir)
        new_dir = zeros(0, 3, 'like', dir);
        return;
    end

    local_z = normalizeRowVectors(dir);
    reference = zeros(size(local_z), 'like', local_z);
    useZAxis = abs(local_z(:, 3)) < 0.9999;
    reference(useZAxis, 3) = 1;
    reference(~useZAxis, 2) = 1;

    local_x = rowwiseCross(reference, local_z);
    fallbackMask = sqrt(sum(local_x .^ 2, 2)) <= eps;
    if any(gather(fallbackMask))
        fallbackReference = zeros(sum(gather(fallbackMask)), 3, 'like', local_z);
        fallbackReference(:, 1) = 1;
        local_x(fallbackMask, :) = rowwiseCross(fallbackReference, local_z(fallbackMask, :));
    end
    local_x = normalizeRowVectors(local_x);
    local_y = normalizeRowVectors(rowwiseCross(local_z, local_x));

    sinTheta = sin(theta(:));
    cosTheta = cos(theta(:));
    cosPhi = cos(phi(:));
    sinPhi = sin(phi(:));
    localCoeffX = sinTheta .* cosPhi;
    localCoeffY = sinTheta .* sinPhi;

    new_dir = local_x .* localCoeffX + local_y .* localCoeffY + local_z .* cosTheta;
    new_dir = normalizeRowVectors(new_dir);
end

function crossProduct = rowwiseCross(a, b)
    crossProduct = [ ...
        a(:, 2) .* b(:, 3) - a(:, 3) .* b(:, 2), ...
        a(:, 3) .* b(:, 1) - a(:, 1) .* b(:, 3), ...
        a(:, 1) .* b(:, 2) - a(:, 2) .* b(:, 1)];
end

function normalized = normalizeRowVectors(vectors)
    if isempty(vectors)
        normalized = vectors;
        return;
    end

    norms = sqrt(sum(vectors .^ 2, 2));
    norms = max(norms, eps);
    normalized = vectors ./ norms;
end

function excitationChunk = buildExcitationChunk(start_pos, end_pos, accumulated_before, depositedEnergy, Eg, num_excitons, entryScale)
    excitationRows = sum(num_excitons);
    excitationChunk = zeros(excitationRows, 4);
    if excitationRows < 1
        return;
    end

    validIdx = find(num_excitons > 0);
    repeatedIdx = repelem(validIdx(:), num_excitons(validIdx));
    groupCounts = num_excitons(validIdx);
    ordinalCells = arrayfun(@(count) (1:count)', groupCounts, 'UniformOutput', false);
    localOrdinal = vertcat(ordinalCells{:});
    repeatedIdx = repeatedIdx(:);
    repeatedEg = Eg(repeatedIdx);
    repeatedAccumulated = accumulated_before(repeatedIdx);
    repeatedDepositedEnergy = depositedEnergy(repeatedIdx);
    energyLevels = localOrdinal .* (3 .* repeatedEg);
    tValues = (energyLevels - repeatedAccumulated) ./ max(abs(repeatedDepositedEnergy), eps);
    tValues = max(0, min(1, tValues));
    segmentDelta = end_pos(repeatedIdx, :) - start_pos(repeatedIdx, :);
    excitationChunk(:, 1:3) = start_pos(repeatedIdx, :) + bsxfun(@times, tValues(:), segmentDelta);
    excitationChunk(:, 4) = entryScale(repeatedIdx);
end

function Eg = get_bandgap(layer_idx, layers, Eg_GaN, Eg_InN)
    In_comp = layers{layer_idx, 2};
    Eg = Eg_InN*In_comp + Eg_GaN*(1-In_comp) - 1.43*In_comp*(1-In_comp);
end

function layer_idx = find_layer(z, layerBoundaries)
    layer_idx = find(z <= layerBoundaries(2:end - 1), 1);
    if isempty(layer_idx)
        layer_idx = length(layerBoundaries(2:end - 1));
    end
end

function [wavelengths,currentLayer_0,numRays1, emissionPhotonEnergy_eV] = GetWavelength(sourceData,layers, layerBoundaries, Eg_GaN, Eg_InN, params)
    if nargin < 6 || isempty(params)
        params = struct();
    end
    numRays1 = size(sourceData,1);
    wavelengths = zeros(numRays1, 1);
    currentLayer_0 = zeros(numRays1, 1);
    emissionPhotonEnergy_eV = zeros(numRays1, 1);
    for rayIdx = 1:numRays1
        Z = sourceData(rayIdx, 3);
        currentLayer =  find(Z <= layerBoundaries(2:end - 1), 1);
        if isempty(currentLayer)
            currentLayer = size(layers, 1);
        end
        layerName = layers{currentLayer, 1};

        currentLayer_0(rayIdx) = currentLayer;

        if ~strcmp(layerName, 'Substrate')
            spectrumModel = resolveEmissionSpectrumModel(sourceData(rayIdx, :), currentLayer, layers, params);
            lambda = sampleEmissionWavelength(spectrumModel);
        else
            In_comp = layers{currentLayer, 2};
            bandgap = Eg_InN * In_comp + Eg_GaN*(1 - In_comp) - 1.43 * In_comp * (1 - In_comp);
            lambda_center = 1240e-9 / bandgap;
            fwhm = 30e-9;    
            gamma = fwhm / 2;
            u = rand();
            cauchy_sample = gamma * tan(pi*(u - 0.5));
            lambda = lambda_center + cauchy_sample;
            lambda = max(min(lambda, lambda_center + 3*fwhm), lambda_center - 3*fwhm);
        end

        basePhotonEnergy_eV = 1240e-9 / lambda;
        layerBandgap_eV = get_bandgap(currentLayer, layers, Eg_GaN, Eg_InN);

        if size(sourceData, 2) >= 5
            localBandgap_eV = sourceData(rayIdx, 4);
            hotShift_eV = sourceData(rayIdx, 5);
            if isfinite(localBandgap_eV)
                basePhotonEnergy_eV = basePhotonEnergy_eV + (localBandgap_eV - layerBandgap_eV);
            end
            if isfinite(hotShift_eV)
                basePhotonEnergy_eV = basePhotonEnergy_eV + hotShift_eV;
            end
        end

        basePhotonEnergy_eV = max(basePhotonEnergy_eV, 0.05);
        lambda = 1240e-9 / basePhotonEnergy_eV;
        wavelengths(rayIdx) = lambda;   
        emissionPhotonEnergy_eV(rayIdx) = basePhotonEnergy_eV;
    end
end

function spectrumModel = resolveEmissionSpectrumModel(sourceRow, currentLayer, layers, params)
    layerName = string(layers{currentLayer, 1});
    sourceTypeCode = inferEmissionSourceType(sourceRow, currentLayer, layers);
    fallbackSpectrum = layers{currentLayer, 8};

    switch sourceTypeCode
        case 1
            spectrumModel = resolveLayerGroupSpectrum(layers, 'n-GaN', fallbackSpectrum);
        case 2
            spectrumModel = resolveLayerGroupSpectrum(layers, 'Prestrained', fallbackSpectrum);
        case 3
            spectrumModel = resolveLayerGroupSpectrum(layers, 'MQW-Barrier', fallbackSpectrum);
        case 4
            spectrumModel = resolveLayerGroupSpectrum(layers, 'MQW-Well', fallbackSpectrum);
        case 5
            spectrumModel = getCalibrationSpectrum(params, 'semipolarSpectrum', []);
            if isempty(spectrumModel)
                referenceSpectrum = resolveLayerGroupSpectrum(layers, 'MQW-Well', fallbackSpectrum);
                referenceBandgap = get_bandgap(currentLayer, layers, 3.3032, 0.6086);
                localBandgap = referenceBandgap;
                if numel(sourceRow) >= 4 && isfinite(sourceRow(4))
                    localBandgap = sourceRow(4);
                end
                spectrumModel = shiftSpectrumModelByEnergy(referenceSpectrum, max(0, localBandgap - referenceBandgap));
            end
        case 6
            spectrumModel = resolveLayerGroupSpectrum(layers, 'p-EBL', fallbackSpectrum);
        otherwise
            if contains(layerName, 'MQW-Well')
                spectrumModel = resolveLayerGroupSpectrum(layers, 'MQW-Well', fallbackSpectrum);
            elseif contains(layerName, 'MQW-Barrier')
                spectrumModel = resolveLayerGroupSpectrum(layers, 'MQW-Barrier', fallbackSpectrum);
            elseif contains(layerName, 'Prestrained')
                spectrumModel = resolveLayerGroupSpectrum(layers, 'Prestrained', fallbackSpectrum);
            elseif contains(layerName, 'p-EBL') || contains(layerName, 'p-GaN')
                spectrumModel = resolveLayerGroupSpectrum(layers, 'p-EBL', fallbackSpectrum);
            elseif contains(layerName, 'GaN')
                spectrumModel = resolveLayerGroupSpectrum(layers, 'n-GaN', fallbackSpectrum);
            else
                spectrumModel = fallbackSpectrum;
            end
    end

    spectrumModel = applyPixelLocalSpectrumShift(spectrumModel, sourceTypeCode, params);
end

function sourceTypeCode = inferEmissionSourceType(sourceRow, currentLayer, layers)
    sourceTypeCode = 0;
    if numel(sourceRow) >= 6 && isfinite(sourceRow(6))
        candidate = round(sourceRow(6));
        if any(candidate == [1, 2, 3, 4, 5, 6])
            sourceTypeCode = candidate;
            return;
        end
    end

    layerName = string(layers{currentLayer, 1});
    if contains(layerName, 'p-EBL') || contains(layerName, 'p-GaN')
        sourceTypeCode = 6;
    elseif contains(layerName, 'MQW-Well')
        sourceTypeCode = 4;
    elseif contains(layerName, 'MQW-Barrier')
        sourceTypeCode = 3;
    elseif contains(layerName, 'Prestrained')
        sourceTypeCode = 2;
    elseif contains(layerName, 'GaN')
        sourceTypeCode = 1;
    end
end

function spectrumModel = resolveLayerGroupSpectrum(layers, targetLayerName, fallbackSpectrum)
    rowIdx = find(strcmp(layers(:, 1), targetLayerName), 1, 'first');
    if isempty(rowIdx)
        layerNames = string(layers(:, 1));
        rowIdx = find(contains(layerNames, string(targetLayerName), 'IgnoreCase', true), 1, 'first');
    end
    if isempty(rowIdx)
        spectrumModel = fallbackSpectrum;
    else
        spectrumModel = layers{rowIdx, 8};
    end
end

function spectrumModel = getCalibrationSpectrum(params, fieldName, defaultValue)
    spectrumModel = defaultValue;
    if isstruct(params) && isfield(params, 'calibrationProfile') && ...
            isstruct(params.calibrationProfile) && isfield(params.calibrationProfile, fieldName) && ...
            ~isempty(params.calibrationProfile.(fieldName))
        spectrumModel = params.calibrationProfile.(fieldName);
    end
end

function pixelParams = resolvePixelParameterContext(params, pixelIndex)
    pixelParams = params;
    pixelParams.pixelLocalContext = struct( ...
        'nGaNEnergyShift_eV', 0, ...
        'gaNEnergyShift_eV', 0, ...
        'prestrainedEnergyShift_eV', 0, ...
        'barrierEnergyShift_eV', 0, ...
        'wellEnergyShift_eV', 0, ...
        'semipolarEnergyShift_eV', 0, ...
        'pTypeEnergyShift_eV', 0, ...
        'nGaNSourceGain', 1, ...
        'gaNSourceGain', 1, ...
        'prestrainedSourceGain', 1, ...
        'barrierSourceGain', 1, ...
        'wellSourceGain', 1, ...
        'semipolarSourceGain', 1, ...
        'pTypeSourceGain', 1, ...
        'targetFamilyFractions', normalizeFamilyFractionsLocal([1, 1, 1, 1, 1, 1]), ...
        'targetPeakNmByFamily', [376; 391; 480; 540; 510; 388], ...
        'isVpitLikePixel', false, ...
        'vpitSuppressionStrength', 0, ...
        'vpitLongwaveCompetition', 0, ...
        'vpitMqwRecovery', 0, ...
        'highVoltageMqwRecovery', 0, ...
        'totalIntensityGain', 1, ...
        'vpitCollectionGain', 1, ...
        'centerFillGain', 1, ...
        'shellThicknessScale', 1, ...
        'vpitInfluence', 0, ...
        'regionCategory', 1, ...
        'shortwaveToMainRatio', 0, ...
        'shoulderToMainRatio', 0, ...
        'pTypeShare', 0.34, ...
        'wellInComposition', NaN, ...
        'barrierInComposition', NaN, ...
        'semipolarInComposition', NaN);

    if ~isstruct(params) || ~isfield(params, 'localParameterField') || isempty(params.localParameterField) || ...
            ~isfield(params, 'pixelsX') || ~isfield(params, 'pixelsY')
        return;
    end

    localField = params.localParameterField;
    [rowIdx, colIdx] = ind2sub([params.pixelsY, params.pixelsX], pixelIndex);
    if rowIdx < 1 || colIdx < 1 || rowIdx > localField.rows || colIdx > localField.cols
        return;
    end

    globalWellIn = getFieldOrDefault(localField, 'globalWellIn', NaN);
    globalBarrierIn = getFieldOrDefault(localField, 'globalBarrierIn', NaN);
    globalSemipolarIn = getFieldOrDefault(localField, 'globalSemipolarIn', NaN);
    totalIntensityGainMap = getFieldOrDefault(localField, 'totalIntensityGainMap', ones(localField.rows, localField.cols));
    vpitCollectionGainMap = getFieldOrDefault(localField, 'vpitCollectionGainMap', ones(localField.rows, localField.cols));
    centerFillGainMap = getFieldOrDefault(localField, 'centerFillGainMap', ones(localField.rows, localField.cols));
    vpitInfluenceMap = getFieldOrDefault(localField, 'vpitInfluenceMap', zeros(localField.rows, localField.cols));
    regionCategoryMap = getFieldOrDefault(localField, 'regionCategoryMap', ones(localField.rows, localField.cols));
    shortMainRatioMap = getFieldOrDefault(localField, 'shortMainRatioMap', zeros(localField.rows, localField.cols));
    shoulderMainRatioMap = getFieldOrDefault(localField, 'shoulderMainRatioMap', zeros(localField.rows, localField.cols));
    pTypeShareMap = getFieldOrDefault(localField, 'pTypeShareMap', 0.34 * ones(localField.rows, localField.cols));
    localWellIn = localField.wellInMap(rowIdx, colIdx);
    localBarrierIn = localField.barrierInMap(rowIdx, colIdx);
    localSemipolarIn = localField.semipolarInMap(rowIdx, colIdx);
    targetFamilyFractions = normalizeFamilyFractionsLocal([ ...
        localField.nGaNAreaMap(rowIdx, colIdx), ...
        localField.prestrainedAreaMap(rowIdx, colIdx), ...
        localField.barrierAreaMap(rowIdx, colIdx), ...
        localField.wellAreaMap(rowIdx, colIdx), ...
        localField.semipolarAreaMap(rowIdx, colIdx), ...
        localField.pTypeAreaMap(rowIdx, colIdx)]);
    targetPeakNmByFamily = [ ...
        localField.nGaNPeakMap_nm(rowIdx, colIdx); ...
        localField.prestrainedPeakMap_nm(rowIdx, colIdx); ...
        localField.barrierPeakMap_nm(rowIdx, colIdx); ...
        localField.wellPeakMap_nm(rowIdx, colIdx); ...
        localField.semipolarPeakMap_nm(rowIdx, colIdx); ...
        localField.pTypePeakMap_nm(rowIdx, colIdx)];
    centerFillGain = centerFillGainMap(rowIdx, colIdx);
    vpitInfluence = clamp(vpitInfluenceMap(rowIdx, colIdx), 0, 1);
    regionCategory = regionCategoryMap(rowIdx, colIdx);
    shortMainRatio = shortMainRatioMap(rowIdx, colIdx);
    shoulderMainRatio = shoulderMainRatioMap(rowIdx, colIdx);
    if ~isfinite(shortMainRatio)
        shortMainRatio = 0;
    end
    if ~isfinite(shoulderMainRatio)
        shoulderMainRatio = 0;
    end
    shellThicknessScale = localField.shellThicknessScaleMap(rowIdx, colIdx);
    electronEnergy_keV = getFieldOrDefault(params, 'electronEnergy', NaN) / 1e3;
    if ~isfinite(electronEnergy_keV) || electronEnergy_keV <= 0
        electronEnergy_keV = getFieldOrDefault(params, 'electronEnergy_keV', 6);
    end
    if ~isfinite(electronEnergy_keV) || electronEnergy_keV <= 0
        electronEnergy_keV = 6;
    end
    highVoltageMqwRecovery = clamp(1 ./ (1 + exp(-(electronEnergy_keV - 7.2) ./ 0.75)), 0, 1);
    rawVpitSuppressionStrength = ...
        0.32 * vpitInfluence + ...
        0.40 * min(max(shortMainRatio, 0), 1.35) + ...
        0.14 * min(max(shoulderMainRatio, 0), 0.85) + ...
        0.18 * max(0, shellThicknessScale - 1) - ...
        0.06 * max(0, centerFillGain - 1);
    vpitSuppressionStrength = clamp(rawVpitSuppressionStrength * (1 - 0.34 * highVoltageMqwRecovery), ...
        0, 0.90);
    rawVpitLongwaveCompetition = ...
        0.12 * vpitInfluence + ...
        0.26 * min(max(shortMainRatio, 0), 1.35) + ...
        0.10 * min(max(shoulderMainRatio, 0), 0.85) + ...
        0.10 * max(0, shellThicknessScale - 1);
    vpitLongwaveCompetition = clamp(rawVpitLongwaveCompetition * (1 - 0.58 * highVoltageMqwRecovery), ...
        0, 0.52);
    vpitMqwRecovery = clamp(highVoltageMqwRecovery * ...
        (0.42 * vpitInfluence + 0.26 * max(0, centerFillGain - 1) + 0.22 * max(0, shellThicknessScale - 1)), ...
        0, 0.65);

    pixelParams.pixelLocalContext = struct( ...
        'nGaNEnergyShift_eV', energyShiftFromPeakNm(localField.nGaNPeakMap_nm(rowIdx, colIdx), getFieldOrDefault(localField, 'globalNGaNPeakNm', NaN)), ...
        'gaNEnergyShift_eV', energyShiftFromPeakNm(localField.nGaNPeakMap_nm(rowIdx, colIdx), getFieldOrDefault(localField, 'globalNGaNPeakNm', NaN)), ...
        'prestrainedEnergyShift_eV', energyShiftFromPeakNm(localField.prestrainedPeakMap_nm(rowIdx, colIdx), getFieldOrDefault(localField, 'globalPrestrainedPeakNm', NaN)), ...
        'barrierEnergyShift_eV', bandgapFromCompositionLocal(localBarrierIn) - bandgapFromCompositionLocal(globalBarrierIn), ...
        'wellEnergyShift_eV', bandgapFromCompositionLocal(localWellIn) - bandgapFromCompositionLocal(globalWellIn), ...
        'semipolarEnergyShift_eV', bandgapFromCompositionLocal(localSemipolarIn) - bandgapFromCompositionLocal(globalSemipolarIn), ...
        'pTypeEnergyShift_eV', energyShiftFromPeakNm(localField.pTypePeakMap_nm(rowIdx, colIdx), getFieldOrDefault(localField, 'globalPTypePeakNm', NaN)), ...
        'nGaNSourceGain', localField.nGaNGainMap(rowIdx, colIdx), ...
        'gaNSourceGain', localField.nGaNGainMap(rowIdx, colIdx), ...
        'prestrainedSourceGain', localField.prestrainedGainMap(rowIdx, colIdx), ...
        'barrierSourceGain', localField.barrierGainMap(rowIdx, colIdx), ...
        'wellSourceGain', localField.wellGainMap(rowIdx, colIdx), ...
        'semipolarSourceGain', localField.semipolarGainMap(rowIdx, colIdx), ...
        'pTypeSourceGain', localField.pTypeGainMap(rowIdx, colIdx), ...
        'targetFamilyFractions', targetFamilyFractions, ...
        'targetPeakNmByFamily', targetPeakNmByFamily, ...
        'isVpitLikePixel', vpitInfluence > 0.35 || localField.shellThicknessScaleMap(rowIdx, colIdx) > 1.04 || localField.semipolarGainMap(rowIdx, colIdx) > 1.08, ...
        'vpitSuppressionStrength', vpitSuppressionStrength, ...
        'vpitLongwaveCompetition', vpitLongwaveCompetition, ...
        'vpitMqwRecovery', vpitMqwRecovery, ...
        'highVoltageMqwRecovery', highVoltageMqwRecovery, ...
        'totalIntensityGain', totalIntensityGainMap(rowIdx, colIdx), ...
        'vpitCollectionGain', vpitCollectionGainMap(rowIdx, colIdx), ...
        'centerFillGain', centerFillGain, ...
        'shellThicknessScale', shellThicknessScale, ...
        'vpitInfluence', vpitInfluence, ...
        'regionCategory', regionCategory, ...
        'shortwaveToMainRatio', shortMainRatio, ...
        'shoulderToMainRatio', shoulderMainRatio, ...
        'pTypeShare', pTypeShareMap(rowIdx, colIdx), ...
        'wellInComposition', localWellIn, ...
        'barrierInComposition', localBarrierIn, ...
        'semipolarInComposition', localSemipolarIn);
end

function Eg = bandgapFromCompositionLocal(inComposition)
    if ~isfinite(inComposition)
        Eg = NaN;
        return;
    end
    EgGaN = 3.3032;
    EgInN = 0.6086;
    Eg = EgInN * inComposition + EgGaN * (1 - inComposition) - 1.43 * inComposition * (1 - inComposition);
end

function energyShift_eV = energyShiftFromPeakNm(localPeakNm, globalPeakNm)
    if ~isfinite(localPeakNm) || ~isfinite(globalPeakNm) || localPeakNm <= 0 || globalPeakNm <= 0
        energyShift_eV = 0;
        return;
    end
    energyShift_eV = (1240 / localPeakNm - 1240 / globalPeakNm);
end

function fractions = normalizeFamilyFractionsLocal(rawFractions)
    fractions = double(rawFractions(:));
    fractions(~isfinite(fractions) | fractions < 0) = 0;
    if ~any(fractions > 0)
        fractions = [0.04; 0.22; 0.14; 0.36; 0.16; 0.08];
    end
    fractions = fractions / sum(fractions);
end

function entropyNorm = computeNormalizedFamilyEntropy(rawFractions)
    fractions = normalizeFamilyFractionsLocal(rawFractions);
    fractions = fractions(fractions > 0);
    if isempty(fractions)
        entropyNorm = 0;
        return;
    end
    entropyValue = -sum(fractions .* log(fractions));
    entropyNorm = clamp(entropyValue / log(6), 0, 1);
end

function spectrumModel = applyPixelLocalSpectrumShift(spectrumModel, sourceTypeCode, params)
    if ~isstruct(params) || ~isfield(params, 'pixelLocalContext') || isempty(params.pixelLocalContext)
        return;
    end

    localContext = params.pixelLocalContext;
    energyShift_eV = 0;
    switch sourceTypeCode
        case 1
            energyShift_eV = localContext.nGaNEnergyShift_eV;
        case 2
            energyShift_eV = localContext.prestrainedEnergyShift_eV;
        case 3
            energyShift_eV = localContext.barrierEnergyShift_eV;
        case 4
            energyShift_eV = localContext.wellEnergyShift_eV;
        case 5
            energyShift_eV = localContext.semipolarEnergyShift_eV;
        case 6
            energyShift_eV = localContext.pTypeEnergyShift_eV;
    end

    if isfinite(energyShift_eV) && abs(energyShift_eV) > 1e-6
        spectrumModel = shiftSpectrumModelByEnergy(spectrumModel, energyShift_eV);
    end
end

function sourceWeights = applyPixelCalibrationRefinement(sourceWeights, sourceData, params, pixelIndex)
    if isstruct(params) && isfield(params, 'pixelLocalContext') && ~isempty(params.pixelLocalContext)
        localContext = params.pixelLocalContext;
        gainByType = [ ...
            localContext.nGaNSourceGain, ...
            localContext.prestrainedSourceGain, ...
            localContext.barrierSourceGain, ...
            localContext.wellSourceGain, ...
            localContext.semipolarSourceGain * localContext.shellThicknessScale, ...
            localContext.pTypeSourceGain];
        sourceTypeCodes = inferEmissionSourceTypeBatch(sourceData);
        for typeCode = 1:6
            mask = sourceTypeCodes == typeCode;
            sourceWeights(mask) = sourceWeights(mask) * gainByType(typeCode);
        end
    end

    if ~isstruct(params) || ~isfield(params, 'calibrationRefinement') || isempty(params.calibrationRefinement)
        return;
    end

    refinement = params.calibrationRefinement;
    requiredFields = {'nGaNGainMap', 'prestrainedGainMap', 'barrierGainMap', 'wellGainMap', 'semipolarGainMap', 'pTypeGainMap'};
    if ~all(isfield(refinement, requiredFields)) || ~isfield(params, 'pixelsX') || ~isfield(params, 'pixelsY')
        return;
    end

    [rowIdx, colIdx] = ind2sub([params.pixelsY, params.pixelsX], pixelIndex);
    if rowIdx < 1 || colIdx < 1 || rowIdx > size(refinement.nGaNGainMap, 1) || colIdx > size(refinement.nGaNGainMap, 2)
        return;
    end

    gainByType = [ ...
        refinement.nGaNGainMap(rowIdx, colIdx), ...
        refinement.prestrainedGainMap(rowIdx, colIdx), ...
        refinement.barrierGainMap(rowIdx, colIdx), ...
        refinement.wellGainMap(rowIdx, colIdx), ...
        refinement.semipolarGainMap(rowIdx, colIdx), ...
        refinement.pTypeGainMap(rowIdx, colIdx)];
    gainByType(~isfinite(gainByType)) = 1;
    gainByType = max(0.3, min(1.8, gainByType));

    sourceTypeCodes = inferEmissionSourceTypeBatch(sourceData);
    for typeCode = 1:6
        mask = sourceTypeCodes == typeCode;
        sourceWeights(mask) = sourceWeights(mask) * gainByType(typeCode);
    end
end

function sourceWeights = applyPixelOneStepBandMatching(sourceWeights, sourceData, sourceLayerIdx, qeWeights, params)
    if ~isstruct(params) || ~isfield(params, 'pixelLocalContext') || isempty(params.pixelLocalContext)
        return;
    end

    localContext = params.pixelLocalContext;
    targetFractions = getFieldOrDefault(localContext, 'targetFamilyFractions', []);
    if isempty(targetFractions) || numel(targetFractions) < 6
        return;
    end

    sourceTypeCodes = inferEmissionSourceTypeBatch(sourceData);
    effectiveWeights = max(sourceWeights(:), 0) .* max(qeWeights(:), 0);
    familyMass = accumarray(sourceTypeCodes, effectiveWeights, [6, 1], @sum, 0);
    if sum(familyMass) <= eps
        return;
    end

    currentFractions = familyMass / sum(familyMass);
    targetFractions = normalizeFamilyFractionsLocal(targetFractions);

    gainByType = ones(6, 1);
    validMask = currentFractions > 1e-6 & targetFractions > 0;
    gainByType(validMask) = (targetFractions(validMask) ./ currentFractions(validMask)) .^ 0.90;

    zeroCurrentMask = currentFractions <= 1e-6 & targetFractions > 0.08;
    gainByType(zeroCurrentMask) = 1.20;

    gainByType = max(0.45, min(1.85, gainByType));

    if getFieldOrDefault(localContext, 'isVpitLikePixel', false)
        shortwaveIdx = [1, 2, 6];
        currentShortwave = sum(currentFractions(shortwaveIdx));
        targetShortwave = sum(targetFractions(shortwaveIdx));
        if currentShortwave > targetShortwave + 0.05
            shortwaveScale = max(0.86, min(1.0, (targetShortwave / max(currentShortwave, eps)) ^ (0.42 + 0.06 * getFieldOrDefault(localContext, 'vpitSuppressionStrength', 0))));
            gainByType(shortwaveIdx) = gainByType(shortwaveIdx) * shortwaveScale;
            gainByType([4, 5]) = min(1.75, gainByType([4, 5]) / max(shortwaveScale ^ 0.16, 0.90));
        elseif targetShortwave > currentShortwave + 0.03
            shortwaveBoost = min(1.70, (targetShortwave / max(currentShortwave, 1e-6)) ^ 0.48);
            gainByType(shortwaveIdx) = gainByType(shortwaveIdx) * shortwaveBoost;
        end

        vpitLongwaveCompetition = clamp(getFieldOrDefault(localContext, 'vpitLongwaveCompetition', 0), 0, 0.60);
        vpitMqwRecovery = clamp(getFieldOrDefault(localContext, 'vpitMqwRecovery', 0), 0, 0.70);
        if vpitLongwaveCompetition > 0
            gainByType(shortwaveIdx) = gainByType(shortwaveIdx) * (1 + 0.56 * vpitLongwaveCompetition);
            gainByType([4, 5]) = gainByType([4, 5]) * max(0.54, 1 - 0.66 * vpitLongwaveCompetition);
            gainByType(3) = gainByType(3) * max(0.72, 1 - 0.30 * vpitLongwaveCompetition);
        end
        if vpitMqwRecovery > 0
            gainByType(4) = min(2.05, gainByType(4) * (1 + 0.42 * vpitMqwRecovery));
            gainByType(5) = min(2.18, gainByType(5) * (1 + 0.30 * vpitMqwRecovery));
            gainByType(3) = min(1.72, gainByType(3) * (1 + 0.12 * vpitMqwRecovery));
        end

        centerFillGain = clamp(getFieldOrDefault(localContext, 'centerFillGain', 1), 0.85, 1.35);
        if centerFillGain > 1.02
            centerBoost = centerFillGain - 1;
            gainByType(4) = min(2.05, gainByType(4) * (1 + 0.14 * centerBoost + 0.18 * vpitMqwRecovery) * max(0.88, 1 - 0.18 * vpitLongwaveCompetition));
            gainByType(5) = min(2.18, gainByType(5) * (1 + 0.18 * centerBoost + 0.12 * vpitMqwRecovery) * max(0.90, 1 - 0.13 * vpitLongwaveCompetition));
            gainByType(shortwaveIdx) = min(2.20, gainByType(shortwaveIdx) * (1 + 0.16 * centerBoost));
        end
    end

    adjustedMass = familyMass .* gainByType;
    if sum(adjustedMass) > eps
        gainByType = gainByType * (sum(familyMass) / sum(adjustedMass));
    end
    gainByType = max(0.35, min(2.10, gainByType));

    for typeCode = 1:6
        mask = sourceTypeCodes == typeCode;
        sourceWeights(mask) = sourceWeights(mask) * gainByType(typeCode);
    end

    if ~isempty(sourceLayerIdx)
        topLayerMask = sourceLayerIdx <= 3;
        if any(topLayerMask) && getFieldOrDefault(localContext, 'isVpitLikePixel', false)
            sourceWeights(topLayerMask) = (1.00 + 0.06 * max(getFieldOrDefault(localContext, 'centerFillGain', 1) - 1, 0) - ...
                0.01 * getFieldOrDefault(localContext, 'vpitSuppressionStrength', 0)) * sourceWeights(topLayerMask);
        end
    end
end

function sourceWeights = applySampledBandMatching(sourceWeights, wavelengths_m, qeWeights, params)
    if ~isstruct(params) || ~isfield(params, 'pixelLocalContext') || isempty(params.pixelLocalContext)
        return;
    end

    localContext = params.pixelLocalContext;
    targetFractions = getFieldOrDefault(localContext, 'targetFamilyFractions', []);
    targetPeaks = getFieldOrDefault(localContext, 'targetPeakNmByFamily', []);
    if isempty(targetFractions) || isempty(targetPeaks) || numel(targetFractions) < 6 || numel(targetPeaks) < 6
        return;
    end

    wavelengths_nm = double(wavelengths_m(:)) * 1e9;
    weights = max(sourceWeights(:), 0) .* max(qeWeights(:), 0);
    if isempty(wavelengths_nm) || sum(weights) <= eps
        return;
    end

    bandIdx = assignWavelengthBands(wavelengths_nm, targetPeaks(:));
    currentFractions = accumarray(bandIdx, weights, [6, 1], @sum, 0);
    currentFractions = currentFractions / max(sum(currentFractions), eps);
    targetFractions = normalizeFamilyFractionsLocal(targetFractions);

    gainByBand = ones(6, 1);
    validMask = currentFractions > 1e-6 & targetFractions > 0;
    gainByBand(validMask) = (targetFractions(validMask) ./ currentFractions(validMask)) .^ 0.95;
    gainByBand = max(0.18, min(2.40, gainByBand));

    if getFieldOrDefault(localContext, 'isVpitLikePixel', false)
        shortwaveIdx = [1, 2, 6];
        currentShortwave = sum(currentFractions(shortwaveIdx));
        targetShortwave = sum(targetFractions(shortwaveIdx));
        if currentShortwave > targetShortwave + 0.05
            shortwaveScale = max(0.84, min(1.0, (targetShortwave / max(currentShortwave, eps)) ^ (0.44 + 0.06 * getFieldOrDefault(localContext, 'vpitSuppressionStrength', 0))));
            gainByBand(shortwaveIdx) = gainByBand(shortwaveIdx) * shortwaveScale;
        elseif targetShortwave > currentShortwave + 0.03
            shortwaveBoost = min(1.62, (targetShortwave / max(currentShortwave, 1e-6)) ^ 0.46);
            gainByBand(shortwaveIdx) = gainByBand(shortwaveIdx) * shortwaveBoost;
        end

        vpitLongwaveCompetition = clamp(getFieldOrDefault(localContext, 'vpitLongwaveCompetition', 0), 0, 0.60);
        vpitMqwRecovery = clamp(getFieldOrDefault(localContext, 'vpitMqwRecovery', 0), 0, 0.70);
        if vpitLongwaveCompetition > 0
            longwaveIdx = [4, 5];
            gainByBand(shortwaveIdx) = gainByBand(shortwaveIdx) * (1 + 0.48 * vpitLongwaveCompetition);
            gainByBand(longwaveIdx) = gainByBand(longwaveIdx) * max(0.58, 1 - 0.58 * vpitLongwaveCompetition);
            gainByBand(3) = gainByBand(3) * max(0.74, 1 - 0.26 * vpitLongwaveCompetition);
        end
        if vpitMqwRecovery > 0
            longwaveIdx = [4, 5];
            gainByBand(longwaveIdx) = min(2.15, gainByBand(longwaveIdx) * (1 + 0.34 * vpitMqwRecovery));
            gainByBand(3) = min(1.70, gainByBand(3) * (1 + 0.10 * vpitMqwRecovery));
        end

        centerFillGain = clamp(getFieldOrDefault(localContext, 'centerFillGain', 1), 0.85, 1.35);
        if centerFillGain > 1.02
            longwaveIdx = [4, 5];
            centerBoost = centerFillGain - 1;
            gainByBand(longwaveIdx) = min(2.12, gainByBand(longwaveIdx) * (1 + 0.14 * centerBoost + 0.14 * vpitMqwRecovery) * max(0.88, 1 - 0.15 * vpitLongwaveCompetition));
            gainByBand(shortwaveIdx) = min(2.10, gainByBand(shortwaveIdx) * (1 + 0.12 * centerBoost));
        end
    end

    sourceWeights = sourceWeights(:) .* gainByBand(bandIdx);
end

function bandIdx = assignWavelengthBands(wavelengths_nm, targetPeaks_nm)
    targetPeaks_nm = double(targetPeaks_nm(:));
    defaultPeaks_nm = [376; 391; 480; 540; 510; 388];
    invalidMask = ~isfinite(targetPeaks_nm) | targetPeaks_nm <= 0;
    targetPeaks_nm(invalidMask) = defaultPeaks_nm(invalidMask);
    [sortedPeaks, sortOrder] = sort(targetPeaks_nm, 'ascend');
    for peakIdx = 2:numel(sortedPeaks)
        if sortedPeaks(peakIdx) <= sortedPeaks(peakIdx - 1) + 0.25
            sortedPeaks(peakIdx) = sortedPeaks(peakIdx - 1) + 0.25;
        end
    end
    boundaries = [-inf; 0.5 * (sortedPeaks(1:end-1) + sortedPeaks(2:end)); inf];
    sortedBandIdx = discretize(wavelengths_nm, boundaries);
    sortedBandIdx(~isfinite(sortedBandIdx)) = 1;
    bandIdx = sortOrder(sortedBandIdx);
    bandIdx = bandIdx(:);
end

function sourceWeights = applyPixelIntensityField(sourceWeights, sourceData, params)
    if ~isstruct(params) || ~isfield(params, 'pixelLocalContext') || isempty(params.pixelLocalContext)
        return;
    end

    localContext = params.pixelLocalContext;
    totalIntensityGain = clamp(getFieldOrDefault(localContext, 'totalIntensityGain', 1), 0.55, 2.40);
    sourceWeights = sourceWeights(:) * totalIntensityGain;

    if ~getFieldOrDefault(localContext, 'isVpitLikePixel', false)
        return;
    end

    vpitCollectionGain = clamp(getFieldOrDefault(localContext, 'vpitCollectionGain', 1), 0.90, 1.60);
    centerFillGain = clamp(getFieldOrDefault(localContext, 'centerFillGain', 1), 0.85, 1.35);
    vpitLongwaveCompetition = clamp(getFieldOrDefault(localContext, 'vpitLongwaveCompetition', 0), 0, 0.60);
    vpitMqwRecovery = clamp(getFieldOrDefault(localContext, 'vpitMqwRecovery', 0), 0, 0.70);
    sourceTypeCodes = inferEmissionSourceTypeBatch(sourceData);
    wellMask = sourceTypeCodes == 4;
    semipolarMask = sourceTypeCodes == 5;
    barrierMask = sourceTypeCodes == 3;
    upperMask = sourceTypeCodes == 1 | sourceTypeCodes == 2 | sourceTypeCodes == 6;

    if any(wellMask)
        wellCollectionScale = clamp( ...
            (0.96 + 0.32 * max(vpitCollectionGain - 1, 0) + 0.18 * max(centerFillGain - 1, 0) + ...
            0.48 * vpitMqwRecovery) * max(0.58, 1 - 0.52 * vpitLongwaveCompetition), ...
            0.50, 1.72);
        sourceWeights(wellMask) = sourceWeights(wellMask) * wellCollectionScale;
    end
    if any(semipolarMask)
        semipolarScale = clamp( ...
            (1.00 + 0.42 * max(vpitCollectionGain - 1, 0) + 0.24 * max(centerFillGain - 1, 0) + ...
            0.30 * vpitMqwRecovery) * max(0.68, 1 - 0.30 * vpitLongwaveCompetition), ...
            0.64, 1.92);
        sourceWeights(semipolarMask) = sourceWeights(semipolarMask) * semipolarScale;
    end
    if any(barrierMask)
        barrierScale = clamp( ...
            (1.00 + 0.18 * max(vpitCollectionGain - 1, 0) + 0.10 * max(centerFillGain - 1, 0)) * ...
            max(0.72, 1 - 0.22 * vpitLongwaveCompetition), ...
            0.70, 1.36);
        sourceWeights(barrierMask) = sourceWeights(barrierMask) * barrierScale;
    end
    if any(upperMask)
        upperScale = clamp( ...
            1.00 + 0.26 * max(vpitCollectionGain - 1, 0) + ...
            0.20 * max(centerFillGain - 1, 0) + 0.58 * vpitLongwaveCompetition, ...
            0.92, 1.86);
        sourceWeights(upperMask) = sourceWeights(upperMask) * upperScale;
    end
end

function qeValue = resolveLayerEmissionQuantumYield(layerIdx, layers, sourceRow, params)
    qeValue = layers{layerIdx, 9};
    wellQeReference = resolveLayerFamilyQuantumYield(layers, 'MQW-Well', max(qeValue, 0.72));
    structuralParams = getStructuralParameterStruct(params);
    shortwaveIqeRatio = getFieldOrDefault(structuralParams, 'ganMqwIqeRatio', 1);
    pTypeIqeRatio = getFieldOrDefault(structuralParams, 'pTypeMqwIqeRatio', shortwaveIqeRatio);
    if ~isfinite(shortwaveIqeRatio)
        shortwaveIqeRatio = 1;
    end
    if ~isfinite(pTypeIqeRatio)
        pTypeIqeRatio = shortwaveIqeRatio;
    end

    sourceTypeCode = inferEmissionSourceType(sourceRow, layerIdx, layers);
    switch sourceTypeCode
        case 1
            qeValue = qeValue * clamp(shortwaveIqeRatio, 0.10, 5.0);
        case 2
            qeValue = qeValue * sqrt(clamp(shortwaveIqeRatio, 0.10, 2.25));
            qeValue = min(qeValue, 4.0 * wellQeReference);
        case 3
            qeValue = qeValue * sqrt(clamp(shortwaveIqeRatio, 0.10, 5.0));
            qeValue = min(qeValue, 0.45 * wellQeReference);
        case 6
            qeValue = qeValue * clamp(pTypeIqeRatio, 0.10, 5.0);
    end
    qeValue = max(qeValue, 0);
end

function referenceQe = resolveLayerFamilyQuantumYield(layers, layerNamePattern, fallbackValue)
    referenceQe = fallbackValue;
    if isempty(layers)
        return;
    end

    layerNames = string(layers(:, 1));
    mask = contains(layerNames, string(layerNamePattern), 'IgnoreCase', true);
    if ~any(mask)
        return;
    end

    candidateValues = cell2mat(layers(mask, 9));
    candidateValues = candidateValues(isfinite(candidateValues) & candidateValues > 0);
    if ~isempty(candidateValues)
        referenceQe = max(candidateValues);
    end
end

function adaptiveNumRays = resolveAdaptivePixelRayBudget(baseNumRays, params)
    adaptiveNumRays = max(16, round(baseNumRays));
    if ~isstruct(params) || ~isfield(params, 'pixelLocalContext') || isempty(params.pixelLocalContext)
        return;
    end

    localContext = params.pixelLocalContext;
    familyEntropy = computeNormalizedFamilyEntropy(getFieldOrDefault(localContext, 'targetFamilyFractions', []));
    isVpitLike = getFieldOrDefault(localContext, 'isVpitLikePixel', false);
    totalIntensityGain = clamp(getFieldOrDefault(localContext, 'totalIntensityGain', 1), 0.6, 2.4);
    vpitCollectionGain = clamp(getFieldOrDefault(localContext, 'vpitCollectionGain', 1), 0.9, 1.6);
    shellThicknessScale = clamp(getFieldOrDefault(localContext, 'shellThicknessScale', 1), 0.75, 1.35);

    factor = 1 ...
        + 0.12 * familyEntropy ...
        + 0.10 * max(totalIntensityGain - 1, 0) ...
        + 0.10 * max(vpitCollectionGain - 1, 0) ...
        + 0.08 * max(shellThicknessScale - 1, 0);
    if isVpitLike
        factor = factor + 0.14;
    end
    factor = clamp(factor, 1.0, 1.45);

    adaptiveNumRays = round(baseNumRays * factor);
    adaptiveNumRays = max(baseNumRays, adaptiveNumRays);
    adaptiveNumRays = min(adaptiveNumRays, max(baseNumRays + 24, round(1.45 * baseNumRays)));
end

function sourceTypeCodes = inferEmissionSourceTypeBatch(sourceData)
    sourceTypeCodes = zeros(size(sourceData, 1), 1);
    if size(sourceData, 2) >= 6
        sourceTypeCodes = round(sourceData(:, 6));
    end
    invalidMask = ~ismember(sourceTypeCodes, [1, 2, 3, 4, 5, 6]);
    sourceTypeCodes(invalidMask) = 4;
end

function shiftedSpectrum = shiftSpectrumModelByEnergy(baseSpectrum, energyShift_eV)
    if abs(energyShift_eV) <= 1e-9
        shiftedSpectrum = baseSpectrum;
        return;
    end

    components = normalizeSpectrumComponents(baseSpectrum);
    photonEnergy_eV = 1240e-9 ./ components(:, 2);
    shiftedCenters = 1240e-9 ./ max(photonEnergy_eV + energyShift_eV, 0.05);
    modelType = 'gaussian_mixture';
    if isstruct(baseSpectrum) && isfield(baseSpectrum, 'modelType') && ~isempty(baseSpectrum.modelType)
        modelType = char(string(baseSpectrum.modelType));
    elseif size(components, 2) >= 4
        modelType = 'pseudo_voigt_mixture';
    end
    shiftedSpectrum = struct( ...
        'modelType', modelType, ...
        'components', [components(:, 1), shiftedCenters, components(:, 3), ensureSpectrumEtaColumn(components)]);
end

function lambda = sampleEmissionWavelength(spectrumModel)
    components = normalizeSpectrumComponents(spectrumModel);
    componentWeights = components(:, 1);
    cumulativeWeights = cumsum(componentWeights);
    randomValue = rand();
    componentIdx = find(randomValue <= cumulativeWeights, 1, 'first');
    if isempty(componentIdx)
        componentIdx = size(components, 1);
    end

    centerLambda = components(componentIdx, 2);
    fwhmLambda = components(componentIdx, 3);
    etaValue = ensureSpectrumEtaColumn(components);
    etaValue = etaValue(componentIdx);

    if rand() <= etaValue
        sigmaLambda = fwhmLambda / (2 * sqrt(2 * log(2)));
        lambda = normrnd(centerLambda, sigmaLambda);
        lambda = clamp(lambda, centerLambda - 4 * sigmaLambda, centerLambda + 4 * sigmaLambda);
    else
        gammaLambda = max(fwhmLambda / 2, 1e-12);
        lambda = centerLambda + gammaLambda * tan(pi * (rand() - 0.5));
        lambda = clamp(lambda, centerLambda - 8 * gammaLambda, centerLambda + 8 * gammaLambda);
    end
end

function components = normalizeSpectrumComponents(spectrumModel)
    if isstruct(spectrumModel) && isfield(spectrumModel, 'components') && ~isempty(spectrumModel.components)
        components = double(spectrumModel.components);
    elseif isnumeric(spectrumModel) && numel(spectrumModel) >= 6
        rawValues = double(spectrumModel(:)');
        componentCount = floor(numel(rawValues) / 3);
        rawValues = rawValues(1:3*componentCount);
        components = reshape(rawValues, 3, [])';
    else
        components = [1, 500e-9, 40e-9];
    end

    components(:,1) = abs(components(:,1));
    if ~any(components(:,1))
        components(:,1) = 1;
    end
    components(:,1) = components(:,1) / sum(components(:,1));
    components(:,3) = max(abs(components(:,3)), 1e-9);
    if size(components, 2) < 4
        components(:, 4) = 0.72;
    else
        components(:, 4) = ensureSpectrumEtaColumn(components);
    end
    [~, order] = sort(components(:,2), 'ascend');
    components = components(order, :);
end

function etaColumn = ensureSpectrumEtaColumn(components)
    if size(components, 2) >= 4
        etaColumn = max(0.05, min(0.95, double(components(:, 4))));
    else
        etaColumn = 0.72 * ones(size(components, 1), 1);
    end
end

function value = getFieldOrDefault(data, fieldName, defaultValue)
    if isstruct(data) && isfield(data, fieldName) && ~isempty(data.(fieldName))
        value = data.(fieldName);
    else
        value = defaultValue;
    end
end

function value = clamp(value, lowerBound, upperBound)
    value = max(lowerBound, min(upperBound, value));
end

function out = mergeStructs(base, override)
    if nargin < 1 || ~isstruct(base)
        base = struct();
    end
    out = base;
    if nargin < 2 || ~isstruct(override)
        return;
    end
    names = fieldnames(override);
    for idx = 1:numel(names)
        fieldName = names{idx};
        if isempty(override.(fieldName))
            continue;
        end
        if isstruct(override.(fieldName)) && isfield(out, fieldName) && isstruct(out.(fieldName))
            out.(fieldName) = mergeStructs(out.(fieldName), override.(fieldName));
        else
            out.(fieldName) = override.(fieldName);
        end
    end
end

function toggles = getDefaultModelToggles()
    toggles = struct( ...
        'enablePEInteractionVolume', true, ...
        'enablePhotonRefractionReflection', true, ...
        'enableMQWInterfaceOptics', true, ...
        'enableAbsorption', true, ...
        'enableSecondaryExcitation', true, ...
        'enableVpitTransportCorrection', true);
end

function toggles = mergeModelToggles(defaultToggles, overrideToggles)
    toggles = mergeStructs(defaultToggles, overrideToggles);
end

function tf = isModelToggleEnabled(params, fieldName, defaultValue)
    tf = defaultValue;
    if isstruct(params) && isfield(params, 'modelToggles') && isstruct(params.modelToggles) && ...
            isfield(params.modelToggles, fieldName) && ~isempty(params.modelToggles.(fieldName))
        tf = logical(params.modelToggles.(fieldName));
    end
end

function structuralParams = getStructuralParameterStruct(params)
    structuralParams = struct();
    if isstruct(params) && isfield(params, 'structuralParameters') && isstruct(params.structuralParameters)
        structuralParams = params.structuralParameters;
    end
end

function value = getStructuralOrModelParameter(params, fieldName, defaultValue)
    value = defaultValue;
    structuralParams = getStructuralParameterStruct(params);
    if isfield(structuralParams, fieldName) && ~isempty(structuralParams.(fieldName))
        value = structuralParams.(fieldName);
        return;
    end
    if isstruct(params) && isfield(params, fieldName) && ~isempty(params.(fieldName))
        value = params.(fieldName);
    end
end

function params = applyStructuralParametersToSimulationParams(params)
    structuralParams = getStructuralParameterStruct(params);
    if isempty(fieldnames(structuralParams))
        return;
    end

    targetActiveThicknessNm = getFieldOrDefault(structuralParams, 'mqwActiveThicknessNm', NaN);
    if isfinite(targetActiveThicknessNm) && targetActiveThicknessNm > 0 && ...
            isfield(params, 'barrierThick') && isfield(params, 'wellThick') && isfield(params, 'mqw_pairs')
        barrierCount = getFieldOrDefault(params, 'mqw_barriers', max(0, params.mqw_pairs - 1));
        currentThicknessNm = params.mqw_pairs * params.wellThick + barrierCount * params.barrierThick;
        if isfinite(currentThicknessNm) && currentThicknessNm > 0
            thicknessScale = clamp(targetActiveThicknessNm / currentThicknessNm, 0.35, 3.0);
            params.barrierThick = max(0.2, params.barrierThick * thicknessScale);
            params.wellThick = max(0.2, params.wellThick * thicknessScale);
        end
    end

    mqwPairs = getFieldOrDefault(structuralParams, 'mqwPairs', NaN);
    if isfinite(mqwPairs) && mqwPairs >= 1
        params.mqw_pairs = max(1, round(mqwPairs));
        params.mqw_barriers = getFieldOrDefault(structuralParams, 'mqwBarriers', max(0, params.mqw_pairs - 1));
        params.mqw_barriers = max(0, round(params.mqw_barriers));
    end
end

function Vpits = applyStructuralParametersToVpits(Vpits, structuralParams)
    if isempty(Vpits) || ~isstruct(structuralParams)
        return;
    end

    targetDepthNm = getFieldOrDefault(structuralParams, 'vpitDepthNm', NaN);
    radiusScale = getFieldOrDefault(structuralParams, 'vpitRadiusScale', NaN);
    targetRadiusNm = getFieldOrDefault(structuralParams, 'vpitEffectiveRadiusNm', NaN);
    sidewallAngleDeg = getFieldOrDefault(structuralParams, 'vpitSidewallAngleDeg', NaN);

    for idx = 1:numel(Vpits)
        if isempty(Vpits{idx}) || ~isstruct(Vpits{idx})
            continue;
        end
        pit = Vpits{idx};
        if isfinite(targetRadiusNm) && targetRadiusNm > 0
            pit.topRadius = targetRadiusNm;
        elseif isfinite(radiusScale) && radiusScale > 0 && isfield(pit, 'topRadius')
            pit.topRadius = pit.topRadius * radiusScale;
        end
        if isfinite(targetDepthNm) && targetDepthNm > 0
            pit.depth = targetDepthNm;
        elseif isfinite(sidewallAngleDeg) && sidewallAngleDeg > 1 && sidewallAngleDeg < 89 && isfield(pit, 'topRadius')
            pit.depth = pit.topRadius / tand(sidewallAngleDeg);
        end
        Vpits{idx} = pit;
    end
end

function layer_idx = findLayerIndicesFast(z, layerBoundaries, numLayers)
    upperBounds = layerBoundaries(2:end-1);
    layer_idx = sum(z(:) > upperBounds(:)', 2) + 1;
    layer_idx = max(1, min(layer_idx, numLayers + 1));
end

function boundaryZ = inferPlanarBoundaryZ(oldZ, newZ, layerBoundaries)
    lowerZ = min(oldZ, newZ);
    upperZ = max(oldZ, newZ);
    candidateBoundaries = layerBoundaries(2:end-1);
    hitIdx = find(candidateBoundaries >= lowerZ - eps & candidateBoundaries <= upperZ + eps, 1, 'first');
    if isempty(hitIdx)
        boundaryZ = newZ;
    else
        boundaryZ = candidateBoundaries(hitIdx);
    end
end

function state = buildBulkStateFromLayerIndex(layerIdx, layers)
    state = struct( ...
        'materialKind', 'air', ...
        'regionType', 'air', ...
        'layerIndex', layerIdx, ...
        'layerName', 'Air', ...
        'inComposition', NaN, ...
        'layerThickness_nm', NaN, ...
        'bandgap_eV', NaN, ...
        'isMQW', false, ...
        'closestFacetDistance_nm', inf, ...
        'nFunc', @(lambda) 1.0, ...
        'alphaFunc', []);

    if layerIdx >= 1 && layerIdx <= size(layers, 1)
        state.materialKind = 'bulk';
        state.regionType = 'bulk';
        state.layerName = layers{layerIdx, 1};
        state.inComposition = layers{layerIdx, 2};
        state.layerThickness_nm = layers{layerIdx, 3};
        state.isMQW = contains(string(state.layerName), "MQW");
        if numel(layers(layerIdx, :)) >= 5 && ~isempty(layers{layerIdx, 5})
            state.nFunc = layers{layerIdx, 5};
        end
        if numel(layers(layerIdx, :)) >= 6 && ~isempty(layers{layerIdx, 6})
            state.alphaFunc = layers{layerIdx, 6};
        end
    end
end

function sensitiveMask = classifyVpitSensitiveSegments(old_pos_m, new_pos_m, Vpits, nominalSurfaceZ_nm)
    if isempty(Vpits)
        sensitiveMask = false(size(old_pos_m, 1), 1);
        return;
    end

    old_pos_nm = old_pos_m * 1e9;
    new_pos_nm = new_pos_m * 1e9;
    points_nm = 0.5 * (old_pos_nm + new_pos_nm);
    z_min = min(old_pos_nm(:, 3), new_pos_nm(:, 3));
    z_max = max(old_pos_nm(:, 3), new_pos_nm(:, 3));
    sensitiveMask = false(size(points_nm, 1), 1);

    for pitIdx = 1:numel(Vpits)
        pit = Vpits{pitIdx};
        if isempty(pit) || ~isfield(pit, 'center') || ~isfield(pit, 'topRadius') || ~isfield(pit, 'depth')
            continue;
        end

        influenceRadius_nm = max(0, pit.topRadius) + 25;
        dx = points_nm(:, 1) - pit.center(1);
        dy = points_nm(:, 2) - pit.center(2);
        nearXY = (dx .* dx + dy .* dy) <= influenceRadius_nm ^ 2;
        pitTop_nm = nominalSurfaceZ_nm + 5;
        pitBottom_nm = nominalSurfaceZ_nm - pit.depth - 5;
        nearZ = z_max >= pitBottom_nm & z_min <= pitTop_nm;
        sensitiveMask = sensitiveMask | (nearXY & nearZ);
    end
end

function sensitiveMask = classifyVpitSensitivePoints(samplePos_nm, Vpits, nominalSurfaceZ_nm)
    if isempty(Vpits)
        sensitiveMask = false(size(samplePos_nm, 1), 1);
        return;
    end

    sensitiveMask = false(size(samplePos_nm, 1), 1);
    z_values = samplePos_nm(:, 3);
    for pitIdx = 1:numel(Vpits)
        pit = Vpits{pitIdx};
        if isempty(pit) || ~isfield(pit, 'center') || ~isfield(pit, 'topRadius') || ~isfield(pit, 'depth')
            continue;
        end

        influenceRadius_nm = max(0, pit.topRadius) + 25;
        dx = samplePos_nm(:, 1) - pit.center(1);
        dy = samplePos_nm(:, 2) - pit.center(2);
        nearXY = (dx .* dx + dy .* dy) <= influenceRadius_nm ^ 2;
        pitTop_nm = nominalSurfaceZ_nm + 5;
        pitBottom_nm = nominalSurfaceZ_nm - pit.depth - 5;
        nearZ = z_values >= pitBottom_nm & z_values <= pitTop_nm;
        sensitiveMask = sensitiveMask | (nearXY & nearZ);
    end
end

function tf = hasActiveRays(active_rays)
    if isa(active_rays, 'gpuArray')
        tf = gather(any(active_rays));
    else
        tf = any(active_rays);
    end
end

function tf = canUseGpuBatching(numItems)
    if nargin < 1
        numItems = 0;
    end
    tf = false;
    if numItems < 128
        return;
    end

    try
        currentDevice = gpuDevice;
        tf = ~isempty(currentDevice) && currentDevice.DeviceSupported;
    catch
        tf = false;
    end
end

function tf = canUseGpuElectronBatching(numElectrons)
    if nargin < 1
        numElectrons = 0;
    end

    tf = numElectrons >= 64 && canUseGpuBatching(numElectrons);
end

function values = randomLike(template, sz)
    if nargin < 2
        sz = size(template);
    end

    values = rand(sz, 'like', template);
end

function value = castNumericState(value, template)
    if isa(template, 'gpuArray') && ~isa(value, 'gpuArray')
        value = gpuArray(double(value));
    elseif ~isa(template, 'gpuArray') && isa(value, 'gpuArray')
        value = gather(value);
    end
end

function value = castLogicalState(value, template)
    value = logical(value);
    if isa(template, 'gpuArray') && ~isa(value, 'gpuArray')
        value = gpuArray(value);
    elseif ~isa(template, 'gpuArray') && isa(value, 'gpuArray')
        value = gather(value);
    end
end

function n = calculate_GaN_refractive_index(lambda, Eg_GaN, Eg_InN)
    lambda_nm = lambda * 1e9;
    y = 0;
    h_nu_eV = 1240 ./ lambda_nm;
    Eg_eV = Eg_InN *y + Eg_GaN *(1 - y) - 1.4*y*(1 - y);
    ratio = h_nu_eV ./ Eg_eV;
   
    a = 9.82661 - 8.21608*y - 31.5902*y^2;
    b = 2.73591 + 0.84249*y - 6.29321*y^2;
    
    term = a * (ratio).^(-2) .* (2 - sqrt(1 + ratio) - sqrt(1 - ratio));
    term(ratio >= 1) = a * (ratio(ratio >= 1)).^(-2) .* (2 - sqrt(1 + ratio(ratio >= 1)));
    
    n_squared = term + b;
    n = sqrt(n_squared);
end

function n_InGaN = calculate_InGaN_refractive_index(lambda, x, Eg_GaN, Eg_InN)
    lambda_nm = lambda * 1e9;
    E_g_InGaN = Eg_GaN*(1 - x) + Eg_InN*x - 1.4*x*(1 - x);
    E_photon = 1240 ./ lambda_nm;
    delta_E = E_g_InGaN - Eg_GaN;
    E_shifted = E_photon - delta_E;
   
    valid_E = E_shifted ~= 0;
    shifted_lambda_nm = zeros(size(lambda_nm));
    shifted_lambda_nm(valid_E) = 1240 ./ E_shifted(valid_E);
    shifted_lambda_nm(~valid_E) = Inf;
    n_InGaN = calculate_GaN_refractive_index(shifted_lambda_nm * 1e-9,Eg_GaN, Eg_InN);
end

function alpha = calculate_InGaN_absorption(lambda, bandgap)
    lambda_nm = lambda * 1e9;
    E = 1240 ./ lambda_nm;
    E_B = 0.993 + 0.719 * bandgap;
    delta_E = 0.06;
    alpha_i = 5e3;
    alpha_0 = 5e6;
    
    exponent = (E_B - E) / delta_E;
    alpha = alpha_i + (alpha_0 - alpha_i) ./ (1 + exp(exponent));
end

function alpha = interpolated_GaN_absorption(lambda)
    persistent absorption_interpolant wavelengths_m_min wavelengths_m_max
    if isempty(absorption_interpolant)
        filename = 'GaN_Abs.xlsx';
        data = readmatrix(filename, 'Range', 'A2:B132');
        wavelengths_nm = data(:, 1);
        absorption = data(:, 2);
        [wavelengths_nm_unique, idx] = unique(wavelengths_nm, 'stable');
        absorption_unique = absorption(idx);

        [wavelengths_nm_sorted, sort_idx] = sort(wavelengths_nm_unique);
        absorption_sorted = absorption_unique(sort_idx);

        wavelengths_m = wavelengths_nm_sorted * 1e-9;
        absorption_interpolant = griddedInterpolant(wavelengths_m, absorption_sorted, 'spline');
        wavelengths_m_min = min(wavelengths_m);
        wavelengths_m_max = max(wavelengths_m);
    end
    
    if isa(lambda, 'gpuArray')
        lambda = gather(lambda);
    end
    lambda = double(lambda);
    lambda = max(min(lambda, wavelengths_m_max), wavelengths_m_min);    
    alpha = absorption_interpolant(lambda) * 1e2;
end

function [beam_positions, pixelsX, pixelsY] = generatePixelGrid(scanX, scanY, totalPixels)
    totalArea = scanX * scanY;
    pixelArea = totalArea / totalPixels;
    
    pixelSize = sqrt(pixelArea);
    pixelsX = ceil(scanX / pixelSize);
    pixelsY = ceil(scanY / pixelSize);

    x_edges = linspace(-scanX/2, scanX/2, pixelsX+1);
    y_edges = linspace(-scanY/2, scanY/2, pixelsY+1);

    x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
    y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
    
    [X, Y] = meshgrid(x_centers, y_centers);
    all_points = [X(:), Y(:)] * 1e-9;

    beam_positions = all_points;
end

function Vpits = generateVpits(scanX, scanY, density, maxDepth)
    numVpits = round(scanX * scanY * density * 1e-6);
    Vpits = cell(numVpits, 1);
    
    for i = 1:numVpits
        centerX = rand() * scanX - scanX/2;
        centerY = rand() * scanY - scanY/2;
        topRadius = 75 + rand() * 25;
        depth = maxDepth + rand() * 50;
        Vpits{i} = struct(...
            'center', [centerX, centerY],...
            'topRadius', topRadius,...
            'depth', depth); 
    end
end

function [inVpit, surfaceZ] = checkVpit(x, y, z, Vpits, nominalSurfaceZ)
    inVpit = false;
    minSurfaceZ = nominalSurfaceZ;

    for i = 1:numel(Vpits)
        v = Vpits{i};
        dx = x - v.center(1);
        dy = y - v.center(2);

        d = nominalSurfaceZ - z;
        if d < 0 || d > v.depth
            continue;
        end
        currentApothem = v.topRadius * (1 - d / v.depth);
        orientationDeg = 0;
        if isfield(v, 'orientationDeg') && ~isempty(v.orientationDeg) && isfinite(v.orientationDeg)
            orientationDeg = mod(v.orientationDeg, 60);
        end
        facetNormalsDeg = orientationDeg + (0:60:300);
        normalDistances = dx .* cosd(facetNormalsDeg) + dy .* sind(facetNormalsDeg);
        maxHexDistance = max(normalDistances);
        inHex = maxHexDistance <= currentApothem + eps;

        if inHex
            inVpit = true;
            depth_ratio = 1 - (maxHexDistance / max(v.topRadius, eps));
            depth_ratio = max(0, min(1, depth_ratio));
            surfaceZ_at_point = nominalSurfaceZ - depth_ratio * v.depth;
            minSurfaceZ = min(minSurfaceZ, surfaceZ_at_point);
        end
    end

    if inVpit
        surfaceZ = minSurfaceZ;
    else
        surfaceZ = nominalSurfaceZ;
    end
end

function result = mainMonteCarloSim(params, layer_params, Eg_GaN, Eg_InN, CL_sources, depth_planes_bins,idx, Vpits, nominalSurfaceZ)
    numRays = params.numRays;
    beamCurrent = params.beamCurrent;
    InteTime = params.InteTime;
    stepSize = 1e-9;
    maxRecursionDepth = 10;
    [layers, layerBoundaries] = Layers(layer_params,params);
    
    for i = 1:size(layers,1)
        In_comp = layers{i,2};
        bandgap = Eg_InN * In_comp + Eg_GaN*(1 - In_comp) - 1.43 * In_comp * (1 - In_comp);
        
        if contains(layers{i,1}, 'Substrate')
            layers{i,5} = @(lambda) sqrt(1 + 1.023798*lambda.^2./(lambda.^2 - (0.0614482e-6)^2) + ...
            1.058264*lambda.^2./(lambda.^2 - (0.110700e-6)^2) + ...
            5.280792*lambda.^2./(lambda.^2 - (17.92656e-6)^2));
        elseif contains(layers{i,1}, 'ITO')
            layers{i,5} = @(lambda) 1.63632 + 0.09713./(lambda*1e9).^2 + (-0.00328)./(lambda*1e9).^4;
        else
            layers{i,5} = @(lambda) calculate_InGaN_refractive_index(lambda, In_comp, Eg_GaN, Eg_InN);
        end
        
        if contains(layers{i,1}, 'ITO') || contains(layers{i,1}, 'p-EBL') || contains(layers{i,1}, 'Prestrained')
            layers{i,6} = @(lambda) 4*pi*1e8.*(0.28147./(lambda*1e9) - 0.34242./((lambda*1e9).^2) + 0.14236./((lambda*1e9).^3) - 0.01916./((lambda*1e9).^4));
        elseif contains(layers{i,1}, 'p-GaN') || contains(layers{i,1}, 'n-GaN')
            layers{i,6} = @(lambda) interpolated_GaN_absorption(lambda);
        elseif contains(layers{i,1}, 'Substrate')
            layers{i,6} = @(lambda) 0;
        else
            layers{i,6} = @(lambda) calculate_InGaN_absorption(lambda, bandgap);
        end
    end
    
    if isempty(CL_sources)
        warning('Pixel %d: No CL sources detected!', idx);
        return;
    end

    pixelParams = resolvePixelParameterContext(params, idx);
    numRays = resolveAdaptivePixelRayBudget(numRays, pixelParams);

    layerIndices = discretize(CL_sources(:,3), layerBoundaries(1:end-1));
    validLayers = unique(layerIndices);
    validLayers1 = validLayers';
    layerCounts = histcounts(layerIndices', [validLayers1, max(validLayers1)+1])';
    totalRays = sum(layerCounts);

    if totalRays == 0
        warning('Pixel %d: No valid rays in layers!', idx);
        return;
    end
    
    ratios = layerCounts / totalRays;
    samplesPerLayer = round(ratios * params.numRays);
    diff = params.numRays - sum(samplesPerLayer);

    if diff ~= 0
        [~, maxIdx] = max(layerCounts);
        samplesPerLayer(maxIdx) = samplesPerLayer(maxIdx) + diff;
    end
    
    new_CL_sources = zeros(max(params.numRays, 1), size(CL_sources, 2));
    newSourceCount = 0;
    for i = 1:length(validLayers)
        layerIdx = validLayers(i);
        mask = (layerIndices == layerIdx);
        layerPoints = CL_sources(mask,:);
        
        if size(layerPoints,1) >= samplesPerLayer(i)
            selectedIdx = randperm(size(layerPoints,1), samplesPerLayer(i));
            selectedPoints = layerPoints(selectedIdx,:);
        else
            remaining = samplesPerLayer(i) - size(layerPoints,1);
            supplement = datasample(CL_sources, remaining, 'Replace',true);
            selectedPoints = [layerPoints; supplement];
        end
        nSelected = size(selectedPoints, 1);
        if nSelected < 1
            continue;
        end
        if newSourceCount + nSelected > size(new_CL_sources, 1)
            new_CL_sources = growMatrixRows(new_CL_sources, newSourceCount + nSelected);
        end
        new_CL_sources(newSourceCount + 1:newSourceCount + nSelected, :) = selectedPoints;
        newSourceCount = newSourceCount + nSelected;
    end
    new_CL_sources = new_CL_sources(1:newSourceCount, :);
    new_CL_sources = ensureVpitShortwaveSourcesRepresented( ...
        new_CL_sources, CL_sources, numRays, pixelParams, Vpits, nominalSurfaceZ);

    selectedSourceData = new_CL_sources(1:min(numRays, end), :);
    sourcePositions = selectedSourceData(:, 1:3);
    sourceWeights = ones(size(selectedSourceData, 1), 1);
    sourceRadiativeYield = ones(size(selectedSourceData, 1), 1);
    if size(selectedSourceData, 2) >= 8
        sourceWeights = max(0.05, selectedSourceData(:, 8));
    end
    if size(selectedSourceData, 2) >= 9
        sourceRadiativeYield = max(0, min(1, selectedSourceData(:, 9)));
    end
    sourceWeights = applyPixelCalibrationRefinement(sourceWeights, selectedSourceData, pixelParams, idx);
    currentLayer_0 = findLayerIndicesFast(selectedSourceData(:,3), layerBoundaries, size(layers, 1));
    QE_weights = zeros(size(selectedSourceData, 1), 1);
    for rayIdx = 1:numel(currentLayer_0)
        QE_weights(rayIdx) = resolveLayerEmissionQuantumYield( ...
            currentLayer_0(rayIdx), layers, selectedSourceData(rayIdx, :), pixelParams);
    end
    sourceWeights = applyPixelOneStepBandMatching(sourceWeights, selectedSourceData, currentLayer_0, QE_weights, pixelParams);
    [wavelengths,currentLayer_0,numRays1, ~] = GetWavelength(selectedSourceData, layers, layerBoundaries, Eg_GaN, Eg_InN, pixelParams);
    sourceWeights = applySampledBandMatching(sourceWeights, wavelengths, QE_weights, pixelParams);
    sourceWeights = applyPixelIntensityField(sourceWeights, selectedSourceData, pixelParams);
    
    elementaryCharge = 1.602e-19;
    QE_weights = zeros(numRays1,1);

    for rayIdx = 1:numRays1
        currentLayer = currentLayer_0(rayIdx);
        QE_weights(rayIdx) = resolveLayerEmissionQuantumYield( ...
            currentLayer, layers, selectedSourceData(rayIdx, :), pixelParams);
    end

    QE_weights = max(QE_weights, 0);
    sourceSamplingWeight = size(CL_sources,1) / max(numRays, 1);
    excitationScale = (beamCurrent * InteTime / (elementaryCharge * 100)) * (1 / 100);
    photonCountWeight = 3 * QE_weights .* sourceSamplingWeight .* excitationScale .* sourceWeights .* sourceRadiativeYield;
    intensities = photonCountWeight;
    
    theta_1 = acos(sqrt(rand(numRays1, 1)));
    phi = 2*pi*rand(numRays1,1);
    directions = [sin(theta_1).*cos(phi), sin(theta_1).*sin(phi), abs(cos(theta_1))];
    polarization = rand(numRays1,1) > 0.5;
    wavelengths1 = wavelengths;
    intensities1 = intensities; 
    
    rayStruct = struct(...
    'dir', [], ...          
    'pos', [], ...         
    'intensity', [], ...    
    'pol', [], ...         
    'depth', [], ...        
    'layer', [], ...       
    'traj', [], ...         
    'lambda', [], ...     
    'isPL', [], ... 
    'isPLed', [], ...
    'initialIntensity', [] ...
    );
   
    rays = repmat(rayStruct, numRays1, 1);
 
    for i = 1:numRays1
        rays(i).dir = directions(i, :);
        rays(i).pos = sourcePositions(i, :);
        rays(i).intensity = intensities(i);
        rays(i).pol = polarization(i);
        rays(i).depth = 0;
        rays(i).layer = currentLayer_0(i);
        rays(i).traj = sourcePositions(i, :);
        rays(i).lambda = wavelengths(i);
        rays(i).isPL = false;
        rays(i).isPLed = false;
        rays(i).initialIntensity = intensities(i);
    end
 
    ray_dir = vertcat(rays.dir);
    ray_pos = vertcat(rays.pos);
    ray_intensity = [rays.intensity]';
    ray_pol = [rays.pol]';
    ray_depth = [rays.depth]';
    ray_layer = [rays.layer]';
    ray_lambda = [rays.lambda]';
    ray_isPL = [rays.isPL]';
    ray_isPLed = [rays.isPLed]';
    ray_initialIntensity = [rays.initialIntensity]';

    useGpuState = canUseGpuBatching(numRays1);
    if useGpuState
        ray_dir = gpuArray(ray_dir);
        ray_pos = gpuArray(ray_pos);
        ray_intensity = gpuArray(ray_intensity);
        ray_pol = gpuArray(ray_pol);
        ray_depth = gpuArray(ray_depth);
        ray_layer = gpuArray(ray_layer);
        ray_lambda = gpuArray(ray_lambda);
        ray_isPL = gpuArray(ray_isPL);
        ray_isPLed = gpuArray(ray_isPLed);
        ray_initialIntensity = gpuArray(ray_initialIntensity);
    end
   
    if useGpuState
        active_rays = gpuArray(true(numRays1, 1));
    else
        active_rays = true(numRays1, 1);
    end
    collectionCapacity = max(128, numRays1 * maxRecursionDepth);
    collectedWavelengths_CL = zeros(collectionCapacity, 1);
    collectedIntensities_CL = zeros(collectionCapacity, 1);
    collectedCLCount = 0;
    collectedWavelengths_PL = zeros(collectionCapacity, 1);
    collectedIntensities_PL = zeros(collectionCapacity, 1);
    collectedPLCount = 0;
    PL_positions = zeros(collectionCapacity, 3);
    PL_wavelengths = zeros(collectionCapacity, 1);
    PL_intensities = zeros(collectionCapacity, 1);
    plSourceCount = 0;
    boundaryEpsilon = 1e-12;
    
    while hasActiveRays(active_rays)
        active_idx = gather(find(active_rays));
        current_pos_gpu = ray_pos(active_idx, :);
        current_dir_gpu = ray_dir(active_idx, :);
        new_pos_gpu = current_pos_gpu + stepSize * current_dir_gpu;

        current_pos = gather(current_pos_gpu);
        new_pos = gather(new_pos_gpu);

        [cross_boundary, ~, new_layer, boundary_data] = checkBoundaryCrossing(...
            current_pos, new_pos, layerBoundaries, layers, Vpits, nominalSurfaceZ, pixelParams);

        refractedRayBatch = [];
        crossingIdxList = find(cross_boundary);
        if ~isempty(crossingIdxList)
            crossingActiveIdx = active_idx(crossingIdxList);
            crossing_dir = gather(current_dir_gpu(crossingIdxList, :));
            crossing_lambda = gather(ray_lambda(crossingActiveIdx));
            crossing_pol = gather(ray_pol(crossingActiveIdx));
            [reflected_dir, refracted_dir, reflection_fraction, refraction_fraction, skip_refraction] = ...
                handleBoundaryCrossing(crossing_dir, crossing_lambda, crossing_pol, ...
                boundary_data(crossingIdxList), layers, pixelParams);
        else
            reflected_dir = zeros(0, 3);
            refracted_dir = zeros(0, 3);
            reflection_fraction = zeros(0, 1);
            refraction_fraction = zeros(0, 1);
            skip_refraction = false(0, 1);
            crossingActiveIdx = zeros(0, 1);
            crossing_lambda = zeros(0, 1);
        end

        for crossingListIdx = 1:numel(crossingIdxList)
            localRayIdx = crossingIdxList(crossingListIdx);
            globalRayIdx = crossingActiveIdx(crossingListIdx);
            crossing = boundary_data(localRayIdx);
            boundaryPoint = crossing.point_m;

            if skip_refraction(crossingListIdx)
                ray_pos(globalRayIdx, :) = castNumericState(boundaryPoint + crossing_dir(crossingListIdx, :) * boundaryEpsilon, ray_pos);
                ray_layer(globalRayIdx) = new_layer(localRayIdx);
                continue;
            end

            reflectedIntensity = gather(ray_intensity(globalRayIdx)) * max(reflection_fraction(crossingListIdx), 0);
            transmittedIntensity = gather(ray_intensity(globalRayIdx)) * max(refraction_fraction(crossingListIdx), 0);

            if reflectedIntensity > 0
                ray_dir(globalRayIdx, :) = castNumericState(reflected_dir(crossingListIdx, :), ray_dir);
                ray_pos(globalRayIdx, :) = castNumericState(boundaryPoint + reflected_dir(crossingListIdx, :) * boundaryEpsilon, ray_pos);
                ray_intensity(globalRayIdx) = reflectedIntensity;
            else
                active_rays(globalRayIdx) = false;
            end

            if transmittedIntensity > 0
                refractedRayBatch = addRefractedRays( ...
                    boundaryPoint + refracted_dir(crossingListIdx, :) * boundaryEpsilon, ...
                    refracted_dir(crossingListIdx, :), ...
                    transmittedIntensity, ...
                    gather(ray_pol(globalRayIdx)), ...
                    gather(ray_depth(globalRayIdx)), ...
                    new_layer(localRayIdx), ...
                    crossing_lambda(crossingListIdx), ...
                    gather(ray_isPL(globalRayIdx)), ...
                    gather(ray_isPLed(globalRayIdx)), ...
                    gather(ray_initialIntensity(globalRayIdx)), ...
                    refractedRayBatch);
            end
        end

        if ~isempty(refractedRayBatch)
            [ray_dir, ray_pos, ray_intensity, ray_pol, ray_depth, ray_layer, ray_lambda, ...
                ray_isPL, ray_isPLed, ray_initialIntensity, active_rays] = appendNewRays( ...
                ray_dir, ray_pos, ray_intensity, ray_pol, ray_depth, ray_layer, ray_lambda, ...
                ray_isPL, ray_isPLed, ray_initialIntensity, active_rays, refractedRayBatch);
        end

        non_boundary_idx = active_idx(~cross_boundary);
        ray_pos(non_boundary_idx, :) = new_pos_gpu(~cross_boundary, :);

        current_active_idx = active_idx(gather(active_rays(active_idx)));
        current_active_pos = gather(ray_pos(current_active_idx, :));
        current_active_dir = gather(ray_dir(current_active_idx, :));
        current_active_intensity = gather(ray_intensity(current_active_idx));
        current_active_lambda = gather(ray_lambda(current_active_idx));
        current_active_depth = gather(ray_depth(current_active_idx));
        current_active_initialIntensity = gather(ray_initialIntensity(current_active_idx));
        current_active_layer = gather(ray_layer(current_active_idx));

        [absorbed, PL_intensity, PL_pos, PL_wavelength, PL_data] = handleAbsorptionAndPL(...
            current_active_pos, current_active_dir, current_active_intensity, ...
            current_active_lambda, current_active_layer, layers, layerBoundaries, stepSize, ...
            Eg_GaN, Eg_InN, layer_params, pixelParams, Vpits, nominalSurfaceZ);
        
        if ~isempty(PL_data.positions)
            nPLData = size(PL_data.positions, 1);
            [PL_positions, PL_intensities, PL_wavelengths] = ensureLightBufferCapacity( ...
                PL_positions, PL_intensities, PL_wavelengths, plSourceCount + nPLData, max(nPLData, 128));
            plRange = plSourceCount + 1:plSourceCount + nPLData;
            PL_positions(plRange, :) = PL_data.positions;
            PL_wavelengths(plRange) = PL_data.wavelengths;
            PL_intensities(plRange) = PL_data.intensities;
            plSourceCount = plSourceCount + nPLData;
        end
       
        ray_intensity(current_active_idx) = ray_intensity(current_active_idx) - castNumericState(absorbed, ray_intensity);
      
        pl_idx = current_active_idx(PL_intensity > 0);
        if ~isempty(pl_idx)
            localPLMask = PL_intensity > 0;
            plRays = addPLRays(PL_pos(localPLMask, :), PL_intensity(localPLMask), PL_wavelength(localPLMask), ...
                current_active_depth(localPLMask), current_active_layer(localPLMask), current_active_initialIntensity(localPLMask));
            [ray_dir, ray_pos, ray_intensity, ray_pol, ray_depth, ray_layer, ray_lambda, ...
                ray_isPL, ray_isPLed, ray_initialIntensity, active_rays] = appendNewRays( ...
                ray_dir, ray_pos, ray_intensity, ray_pol, ray_depth, ray_layer, ray_lambda, ...
                ray_isPL, ray_isPLed, ray_initialIntensity, active_rays, plRays);
        end
       
        terminated = checkTerminationConditions(current_active_pos, ...
            gather(ray_intensity(current_active_idx)), current_active_depth, layerBoundaries, ...
            current_active_initialIntensity, maxRecursionDepth);
        
        if any(terminated)
            term_idx = current_active_idx(terminated);
            term_lambda = gather(ray_lambda(term_idx));
            term_intensity = gather(ray_intensity(term_idx));
            term_isPL = gather(ray_isPL(term_idx));
            [collectedWavelengths_CL, collectedIntensities_CL, collectedCLCount] = appendSpectrumSamplesToBuffer( ...
                collectedWavelengths_CL, collectedIntensities_CL, collectedCLCount, ...
                term_lambda(~term_isPL), term_intensity(~term_isPL));
            [collectedWavelengths_PL, collectedIntensities_PL, collectedPLCount] = appendSpectrumSamplesToBuffer( ...
                collectedWavelengths_PL, collectedIntensities_PL, collectedPLCount, ...
                term_lambda(term_isPL), term_intensity(term_isPL));
            
            active_rays(term_idx) = false;
        end

        surviving_idx = current_active_idx(~terminated);
        ray_depth(surviving_idx) = ray_depth(surviving_idx) + 1;
    end
    
    PL_positions = gather(PL_positions(1:plSourceCount, :));
    PL_wavelengths = gather(PL_wavelengths(1:plSourceCount));
    PL_intensities = gather(PL_intensities(1:plSourceCount));
    
    if isempty(PL_positions)
        PL_positions = zeros(0, 3);
        PL_wavelengths = zeros(0, 1);
        PL_intensities = zeros(0, 1);
    end

    initial_CL_wavelengths = wavelengths1;
    initial_CL_intensities = intensities1;
    collectedWavelengths_CL = collectedWavelengths_CL(1:collectedCLCount);
    collectedIntensities_CL = collectedIntensities_CL(1:collectedCLCount);
    collectedWavelengths_PL = collectedWavelengths_PL(1:collectedPLCount);
    collectedIntensities_PL = collectedIntensities_PL(1:collectedPLCount);
    collectedWavelengths = [collectedWavelengths_CL; collectedWavelengths_PL];
    collectedIntensities = [collectedIntensities_CL; collectedIntensities_PL];
    
    surface_z = layerBoundaries(end-1);
    if isempty(sourcePositions)
        CL_depths = [];
        warning('No CL sources detected. Depth-resolved data will be skipped.');
    else
        CL_depths = (surface_z - sourcePositions(:,3)) * 1e9;
    end
    
    if isempty(CL_depths)
        warndlg('No CL sources detected! Check simulation parameters.',...
        'Data Warning');
        return;
    end
    
    num_depth_bins = size(depth_planes_bins,1);
    
    CL_wavelengths = wavelengths1;
    CL_intensities = intensities1;
    depth_CL = struct('Wavelengths', cell(num_depth_bins,1), 'Intensities', cell(num_depth_bins,1));
    depth_PL = struct('Wavelengths', cell(num_depth_bins,1), 'Intensities', cell(num_depth_bins,1));
    
    for i = 1:length(CL_depths)
        bin_idx = find(CL_depths(i) >= depth_planes_bins(:,1) & CL_depths(i) < depth_planes_bins(:,2), 1);
        if ~isempty(bin_idx)
            depth_CL(bin_idx).Wavelengths = [depth_CL(bin_idx).Wavelengths; CL_wavelengths(i)];
            depth_CL(bin_idx).Intensities = [depth_CL(bin_idx).Intensities; CL_intensities(i)];
        end
    end
    
    if ~isempty(PL_positions)
        PL_depths = (surface_z - PL_positions(:,3)) * 1e9;
    else
        PL_depths = [];
    end
    
    if ~isempty(PL_depths)
        for i = 1:length(PL_depths)
            bin_idx = find(PL_depths(i) >= depth_planes_bins(:,1) & PL_depths(i) < depth_planes_bins(:,2), 1);
            if ~isempty(bin_idx)
                depth_PL(bin_idx).Wavelengths = [depth_PL(bin_idx).Wavelengths; PL_wavelengths(i)];
                depth_PL(bin_idx).Intensities = [depth_PL(bin_idx).Intensities; PL_intensities(i)];
            end
        end
    end
    
    layer_CL = struct('Wavelengths', cell(size(layers,1),1), 'Intensities', cell(size(layers,1),1));
    layer_PL = struct('Wavelengths', cell(size(layers,1),1), 'Intensities', cell(size(layers,1),1));
    
    for i = 1:length(CL_depths)
        layer_idx = find_layer(sourcePositions(i,3), layerBoundaries);
        if layer_idx > 1
            layer_CL(layer_idx).Wavelengths = [layer_CL(layer_idx).Wavelengths; initial_CL_wavelengths(i)];
            layer_CL(layer_idx).Intensities = [layer_CL(layer_idx).Intensities; initial_CL_intensities(i)];
        end
    end
    
    if ~isempty(PL_positions)
        for i = 1:size(PL_positions,1)
            layer_idx = find_layer(PL_positions(i,3), layerBoundaries);
            if layer_idx > 1
                layer_PL(layer_idx).Wavelengths = [layer_PL(layer_idx).Wavelengths; PL_wavelengths(i)];
                layer_PL(layer_idx).Intensities = [layer_PL(layer_idx).Intensities; PL_intensities(i)];
            end
        end
    end
    
    [xi_modelDriven, f_modelDriven] = buildModelDrivenSpectrumEstimate( ...
        selectedSourceData, intensities1, currentLayer_0, layers, Eg_GaN, Eg_InN, pixelParams);
    [xi_sampled, f_sampled] = buildSmoothedSpectrumEstimate(collectedWavelengths, collectedIntensities, pixelParams);
    [xi_total, f_total] = blendModelDrivenAndSampledSpectra( ...
        xi_modelDriven, f_modelDriven, xi_sampled, f_sampled, collectedIntensities, pixelParams);
    result = struct();
    result.totalIntensity = sum(collectedIntensities);
    result.totalSpectrum = [xi_total(:), f_total(:)];
    result.depthIntensity = zeros(1, size(depth_planes_bins,1));
    result.depthSpectra = cell(1, size(depth_planes_bins,1));
    result.CLSources = sourcePositions; 
    result.CLIntensities = intensities;
    result.PLSources = PL_positions; 
    result.PLIntensities = PL_intensities; 
    result.CLWavelengths = CL_wavelengths;
    result.PLWavelengths = PL_wavelengths;
    if useGpuState
        result.executionBackend = 'hybrid_gpu_batch';
    else
        result.executionBackend = 'cpu';
    end

    for bin_idx = 1:size(depth_planes_bins,1)
        total_wl = [depth_CL(bin_idx).Wavelengths; depth_PL(bin_idx).Wavelengths];
        total_int = [depth_CL(bin_idx).Intensities; depth_PL(bin_idx).Intensities];
        result.depthIntensity(bin_idx) = sum(total_int);
        
        if ~isempty(total_wl)
            [xi_depth, f_depth] = buildSmoothedSpectrumEstimate(total_wl, total_int, pixelParams);
            f_depth = applyDetectorSpectralResponse(xi_depth, f_depth, pixelParams);
            result.depthSpectra{bin_idx} = [xi_depth(:), f_depth(:)];
        end
    end
    
    result.layerIntensity = zeros(1, size(layers,1));
    result.layerSpectra = cell(1, size(layers,1));
    for layer_idx = 1:size(layers,1)
        total_wl = [layer_CL(layer_idx).Wavelengths; layer_PL(layer_idx).Wavelengths];
        total_int = [layer_CL(layer_idx).Intensities; layer_PL(layer_idx).Intensities];
        result.layerIntensity(layer_idx) = sum(total_int);
        
        if ~isempty(total_wl)
            [xi_layer, f_layer] = buildSmoothedSpectrumEstimate(total_wl, total_int, pixelParams);
            f_layer = applyDetectorSpectralResponse(xi_layer, f_layer, pixelParams);
            result.layerSpectra{layer_idx} = [xi_layer(:), f_layer(:)];
        end
    end
end

function [xi_nm, spectrumValues] = buildModelDrivenSpectrumEstimate(sourceData, sourceWeights, currentLayerIdx, layers, Eg_GaN, Eg_InN, params)
    xi_nm = linspace(250, 700, 1024)';
    spectrumValues = zeros(size(xi_nm));
    if isempty(sourceData) || isempty(sourceWeights) || isempty(currentLayerIdx)
        return;
    end

    weights = max(double(sourceWeights(:)), 0);
    currentLayerIdx = double(currentLayerIdx(:));
    sourceTypeCodes = inferEmissionSourceTypeBatch(sourceData);
    validMask = isfinite(weights) & weights > 0 & isfinite(currentLayerIdx) & currentLayerIdx >= 1;
    if ~any(validMask)
        return;
    end

    sourceData = sourceData(validMask, :);
    weights = weights(validMask);
    currentLayerIdx = currentLayerIdx(validMask);
    sourceTypeCodes = sourceTypeCodes(validMask);

    layerBandgaps = zeros(size(currentLayerIdx));
    for idx = 1:numel(currentLayerIdx)
        layerBandgaps(idx) = get_bandgap(currentLayerIdx(idx), layers, Eg_GaN, Eg_InN);
    end

    localBandgaps = layerBandgaps;
    if size(sourceData, 2) >= 4
        candidateLocalBandgaps = sourceData(:, 4);
        localValidMask = isfinite(candidateLocalBandgaps);
        localBandgaps(localValidMask) = candidateLocalBandgaps(localValidMask);
    end

    hotShift_eV = zeros(size(layerBandgaps));
    if size(sourceData, 2) >= 5
        hotShift_eV = sourceData(:, 5);
        hotShift_eV(~isfinite(hotShift_eV)) = 0;
    end
    extraShift_eV = clamp(localBandgaps - layerBandgaps + hotShift_eV, -0.18, 0.18);
    extraShift_eV = clampEmissionShiftBySourceType(extraShift_eV, sourceTypeCodes, params);
    shiftBinIdx = round(extraShift_eV / 0.004);

    groupKeys = [sourceTypeCodes(:), currentLayerIdx(:), shiftBinIdx(:)];
    [~, ~, groupIds] = unique(groupKeys, 'rows', 'stable');
    wavelengthAxis_m = xi_nm * 1e-9;
    sampleSpacing_nm = median(diff(xi_nm));

    for groupIdx = 1:max(groupIds)
        groupMask = groupIds == groupIdx;
        if ~any(groupMask)
            continue;
        end

        representativeIdx = find(groupMask, 1, 'first');
        representativeRow = sourceData(representativeIdx, :);
        if numel(representativeRow) >= 4
            representativeRow(4) = NaN;
        end
        if numel(representativeRow) >= 5
            representativeRow(5) = 0;
        end

        currentLayer = currentLayerIdx(representativeIdx);
        baseSpectrumModel = resolveEmissionSpectrumModel(representativeRow, currentLayer, layers, params);

        groupShift = extraShift_eV(groupMask);
        meanShift_eV = mean(groupShift);
        shiftSpread_eV = 0;
        finiteShiftMask = isfinite(groupShift);
        if nnz(finiteShiftMask) >= 2
            shiftSpread_eV = std(groupShift(finiteShiftMask));
        end

        shiftedSpectrumModel = shiftSpectrumModelByEnergy(baseSpectrumModel, meanShift_eV);
        extraFwhm_nm = convertEnergySpreadToExtraFwhmNm(shiftedSpectrumModel, shiftSpread_eV);
        groupProfile = evaluateSpectrumModelOnAxis(wavelengthAxis_m, shiftedSpectrumModel, extraFwhm_nm);
        spectrumValues = spectrumValues + sum(weights(groupMask)) * groupProfile;
    end

    if any(spectrumValues > 0)
        instrumentFwhm_nm = estimateInstrumentResponseFwhmNm(params);
        residualBroadeningFwhm_nm = 0.65 * estimateAdaptiveSpectrumBroadeningFwhmNm(weights, params);
        totalFwhm_nm = sqrt(instrumentFwhm_nm ^ 2 + residualBroadeningFwhm_nm ^ 2);
        spectrumValues = applyInstrumentResponseConvolution(spectrumValues, sampleSpacing_nm, totalFwhm_nm);
    end
    spectrumValues = applyCalibrationWindowSuppression(xi_nm, spectrumValues, params);
end

function [xi_nm, spectrumValues] = blendModelDrivenAndSampledSpectra(modelXi_nm, modelSpectrum, sampledXi_nm, sampledSpectrum, sampledWeights, params)
    if isempty(modelXi_nm) || isempty(modelSpectrum)
        xi_nm = sampledXi_nm;
        spectrumValues = sampledSpectrum;
        spectrumValues = applyDetectorSpectralResponse(xi_nm, spectrumValues, params);
        spectrumValues = applyForwardOpticalMechanismEnvelope(xi_nm, spectrumValues, params);
        return;
    end

    xi_nm = modelXi_nm(:);
    modelSpectrum = max(double(modelSpectrum(:)), 0);
    if isempty(sampledXi_nm) || isempty(sampledSpectrum)
        spectrumValues = modelSpectrum;
        spectrumValues = applyDetectorSpectralResponse(xi_nm, spectrumValues, params);
        spectrumValues = applyForwardOpticalMechanismEnvelope(xi_nm, spectrumValues, params);
        return;
    end

    sampledSpectrumInterp = interp1(sampledXi_nm(:), double(sampledSpectrum(:)), xi_nm, 'linear', 0);
    modelShape = normalizeSpectrumArea(modelSpectrum, xi_nm);
    sampledShape = normalizeSpectrumArea(sampledSpectrumInterp, xi_nm);
    modelShape = normalizeSpectrumArea(applyCalibrationWindowSuppression(xi_nm, modelShape, params), xi_nm);
    sampledShape = normalizeSpectrumArea(applyCalibrationWindowSuppression(xi_nm, sampledShape, params), xi_nm);
    modelWeight = resolveModelDrivenSpectrumBlendWeight(sampledWeights, params);
    blendedShape = modelWeight * modelShape + (1 - modelWeight) * sampledShape;
    finalBroadeningFwhm_nm = estimateFinalSpectrumBroadeningFwhmNm(params);
    if finalBroadeningFwhm_nm > 0
        blendedShape = applyInstrumentResponseConvolution(blendedShape, median(diff(xi_nm)), finalBroadeningFwhm_nm);
        blendedShape = normalizeSpectrumArea(blendedShape, xi_nm);
    end

    targetIntegral = sum(double(sampledWeights(isfinite(sampledWeights) & sampledWeights > 0)));
    if targetIntegral <= eps
        targetIntegral = trapz(xi_nm, modelSpectrum);
    end
    spectrumValues = blendedShape * (targetIntegral / max(trapz(xi_nm, blendedShape), eps));
    spectrumValues = applyDetectorSpectralResponse(xi_nm, spectrumValues, params);
    spectrumValues = applyForwardOpticalMechanismEnvelope(xi_nm, spectrumValues, params);
end

function spectrumShape = normalizeSpectrumArea(spectrumValues, xi_nm)
    spectrumShape = max(double(spectrumValues(:)), 0);
    areaValue = trapz(xi_nm(:), spectrumShape);
    if areaValue <= eps
        return;
    end
    spectrumShape = spectrumShape / areaValue;
end

function modelWeight = resolveModelDrivenSpectrumBlendWeight(sampledWeights, params)
    effectiveSampleCount = estimateEffectiveSpectrumSampleCount(sampledWeights);
    modelWeight = 0.94 - 0.10 * min(effectiveSampleCount / 220, 1);
    if isstruct(params) && isfield(params, 'pixelLocalContext') && ~isempty(params.pixelLocalContext) && ...
            getFieldOrDefault(params.pixelLocalContext, 'isVpitLikePixel', false)
        modelWeight = modelWeight + 0.04;
    end
    modelWeight = clamp(modelWeight, 0.82, 0.98);
end

function spectrumValues = evaluateSpectrumModelOnAxis(wavelengthAxis_m, spectrumModel, extraFwhm_nm)
    components = normalizeSpectrumComponents(spectrumModel);
    if nargin < 3 || ~isfinite(extraFwhm_nm)
        extraFwhm_nm = 0;
    end

    broadenedFwhm_m = sqrt(components(:, 3) .^ 2 + (extraFwhm_nm * 1e-9) .^ 2);
    etaColumn = ensureSpectrumEtaColumn(components);
    wavelengthAxis_m = wavelengthAxis_m(:)';
    spectrumValues = zeros(numel(wavelengthAxis_m), 1);

    for componentIdx = 1:size(components, 1)
        center_m = components(componentIdx, 2);
        fwhm_m = max(broadenedFwhm_m(componentIdx), 1e-9);
        sigma_m = fwhm_m / (2 * sqrt(2 * log(2)));
        gamma_m = max(fwhm_m / 2, 1e-12);
        gaussianPart = exp(-0.5 * ((wavelengthAxis_m - center_m) / sigma_m) .^ 2) / (sigma_m * sqrt(2 * pi));
        lorentzPart = (gamma_m / pi) ./ ((wavelengthAxis_m - center_m) .^ 2 + gamma_m ^ 2);
        componentProfile = etaColumn(componentIdx) * gaussianPart + (1 - etaColumn(componentIdx)) * lorentzPart;
        spectrumValues = spectrumValues + components(componentIdx, 1) * componentProfile(:);
    end
end

function extraFwhm_nm = convertEnergySpreadToExtraFwhmNm(spectrumModel, shiftSpread_eV)
    if ~isfinite(shiftSpread_eV) || shiftSpread_eV <= 1e-6
        extraFwhm_nm = 0;
        return;
    end

    components = normalizeSpectrumComponents(spectrumModel);
    centerLambda_nm = sum(components(:, 1) .* components(:, 2)) * 1e9;
    sigmaLambda_nm = abs((centerLambda_nm ^ 2 / 1240) * shiftSpread_eV);
    extraFwhm_nm = clamp(2.355 * sigmaLambda_nm, 0, 14);
end

function extraShift_eV = clampEmissionShiftBySourceType(extraShift_eV, sourceTypeCodes, params)
    if isempty(extraShift_eV)
        return;
    end

    extraShift_eV = double(extraShift_eV(:));
    sourceTypeCodes = double(sourceTypeCodes(:));
    lowerBounds = [-0.030; -0.030; -0.040; -0.060; -0.050; -0.025];
    upperBounds = [ 0.080;  0.070;  0.075;  0.070;  0.085;  0.050];
    validTypeMask = ismember(sourceTypeCodes, 1:6);
    boundedTypes = sourceTypeCodes;
    boundedTypes(~validTypeMask) = 4;

    lowerBoundVec = lowerBounds(boundedTypes);
    upperBoundVec = upperBounds(boundedTypes);
    extraShift_eV = min(max(extraShift_eV, lowerBoundVec), upperBoundVec);

    if isstruct(params) && isfield(params, 'pixelLocalContext') && ~isempty(params.pixelLocalContext)
        centerFillGain = getFieldOrDefault(params.pixelLocalContext, 'centerFillGain', 1);
        if centerFillGain > 1.02
            centerBoostMask = boundedTypes == 4 | boundedTypes == 5 | boundedTypes == 3;
            extraShift_eV(centerBoostMask) = max(extraShift_eV(centerBoostMask), ...
                lowerBoundVec(centerBoostMask) - 0.010 * (centerFillGain - 1));
        end
    end
end

function spectrumValues = applyCalibrationWindowSuppression(xi_nm, spectrumValues, params)
    spectrumValues = max(double(spectrumValues(:)), 0);
    if isempty(spectrumValues) || ~isstruct(params) || ~isfield(params, 'calibrationProfile') || isempty(params.calibrationProfile)
        return;
    end

    profile = params.calibrationProfile;
    componentsNm = getFieldOrDefault(profile, 'componentsNm', []);
    if isempty(componentsNm) || size(componentsNm, 2) < 3
        return;
    end

    componentCenters_nm = double(componentsNm(:, 2));
    componentWidths_nm = double(componentsNm(:, 3));
    [familyCenters_nm, familyWidths_nm] = getCalibrationProtectedEmissionWindows(profile, params);
    componentCenters_nm = [componentCenters_nm(:); familyCenters_nm(:)];
    componentWidths_nm = [componentWidths_nm(:); familyWidths_nm(:)];
    validMask = isfinite(componentCenters_nm) & isfinite(componentWidths_nm) & componentWidths_nm > 0;
    componentCenters_nm = componentCenters_nm(validMask);
    componentWidths_nm = componentWidths_nm(validMask);
    if isempty(componentCenters_nm)
        return;
    end

    lowerLimit_nm = max(330, min(componentCenters_nm - 1.8 * componentWidths_nm));
    upperLimit_nm = min(700, max(componentCenters_nm + 1.8 * componentWidths_nm));
    lowerSoftWidth_nm = 18;
    upperSoftWidth_nm = 18;
    lowerSuppression = ones(size(spectrumValues));
    upperSuppression = ones(size(spectrumValues));
    lowerMask = xi_nm < lowerLimit_nm;
    upperMask = xi_nm > upperLimit_nm;
    lowerSuppression(lowerMask) = exp(-((lowerLimit_nm - xi_nm(lowerMask)) / lowerSoftWidth_nm) .^ 2);
    upperSuppression(upperMask) = exp(-((xi_nm(upperMask) - upperLimit_nm) / upperSoftWidth_nm) .^ 2);
    spectrumValues = spectrumValues .* lowerSuppression .* upperSuppression;
end

function [familyCenters_nm, familyWidths_nm] = getCalibrationProtectedEmissionWindows(profile, params)
    familyCenters_nm = [
        getFieldOrDefault(profile, 'nGaNPeakNm', getFieldOrDefault(profile, 'gaNPeakNm', 376));
        getFieldOrDefault(profile, 'pTypePeakNm', 388);
        getFieldOrDefault(profile, 'shortwavePeakNm', ...
        getFieldOrDefault(profile, 'prestrainedPeakNm', 390));
        getFieldOrDefault(profile, 'barrierPeakNm', getFieldOrDefault(profile, 'bluePeakNm', 480));
        getFieldOrDefault(profile, 'semipolarPeakNm', 510);
        getFieldOrDefault(profile, 'wellPeakNm', getFieldOrDefault(profile, 'redPeakNm', 540))];
    familyWidths_nm = [18; 18; 22; 28; 30; 32];

    if isstruct(params) && isfield(params, 'pixelLocalContext') && ~isempty(params.pixelLocalContext) && ...
            getFieldOrDefault(params.pixelLocalContext, 'isVpitLikePixel', false)
        familyWidths_nm(1:3) = familyWidths_nm(1:3) + 8;
    end

    validMask = isfinite(familyCenters_nm) & familyCenters_nm >= 330 & familyCenters_nm <= 700;
    familyCenters_nm = familyCenters_nm(validMask);
    familyWidths_nm = familyWidths_nm(validMask);
end

function [xi_nm, spectrumValues] = buildSmoothedSpectrumEstimate(wavelengths_m, intensities, params)
    xi_nm = linspace(250, 700, 1024)';
    spectrumValues = zeros(size(xi_nm));
    if isempty(wavelengths_m) || isempty(intensities)
        return;
    end

    wavelengths_nm = double(wavelengths_m(:)) * 1e9;
    weights = max(double(intensities(:)), 0);
    validMask = isfinite(wavelengths_nm) & isfinite(weights) & weights > 0;
    wavelengths_nm = wavelengths_nm(validMask);
    weights = weights(validMask);
    if isempty(wavelengths_nm)
        return;
    end

    sampleSpacing_nm = median(diff(xi_nm));
    binEdges_nm = [xi_nm(1) - 0.5 * sampleSpacing_nm; ...
        0.5 * (xi_nm(1:end-1) + xi_nm(2:end)); ...
        xi_nm(end) + 0.5 * sampleSpacing_nm];
    binIdx = discretize(wavelengths_nm, binEdges_nm);
    validBinMask = isfinite(binIdx) & binIdx >= 1 & binIdx <= numel(xi_nm);
    rawSpectrum = accumarray(binIdx(validBinMask), weights(validBinMask), [numel(xi_nm), 1], @sum, 0);
    instrumentFwhm_nm = estimateInstrumentResponseFwhmNm(params);
    adaptiveBroadeningFwhm_nm = estimateAdaptiveSpectrumBroadeningFwhmNm(weights, params);
    totalFwhm_nm = sqrt(max(instrumentFwhm_nm, 0) ^ 2 + adaptiveBroadeningFwhm_nm ^ 2);
    spectrumValues = applyInstrumentResponseConvolution(rawSpectrum, sampleSpacing_nm, totalFwhm_nm);
    spectrumValues = applyCalibrationWindowSuppression(xi_nm, spectrumValues, params);
end

function spectrumValues = applyDetectorSpectralResponse(xi_nm, spectrumValues, params)
    spectrumValues = max(double(spectrumValues(:)), 0);
    if isempty(spectrumValues) || ~isstruct(params)
        return;
    end

    if isfield(params, 'detectorResponse') && ~isempty(params.detectorResponse)
        response = interpolateDetectorResponse(xi_nm, params.detectorResponse);
    elseif isfield(params, 'calibrationProfile') && isstruct(params.calibrationProfile) && ...
            isfield(params.calibrationProfile, 'detectorResponse') && ~isempty(params.calibrationProfile.detectorResponse)
        response = interpolateDetectorResponse(xi_nm, params.calibrationProfile.detectorResponse);
    else
        structuralParams = getStructuralParameterStruct(params);
        shortwaveScale = getFieldOrDefault(structuralParams, 'detectorShortwaveResponseScale', 1);
        response = buildDefaultDetectorResponse(xi_nm, shortwaveScale);
    end

    if isempty(response)
        return;
    end
    response = max(double(response(:)), 0);
    if numel(response) ~= numel(spectrumValues) || ~any(response > 0)
        return;
    end
    response = response / max(mean(response(xi_nm >= 500 & xi_nm <= 540), 'omitnan'), eps);
    response(~isfinite(response)) = 1;
    spectrumValues = spectrumValues .* response;
end

function spectrumValues = applyForwardOpticalMechanismEnvelope(xi_nm, spectrumValues, params)
    spectrumValues = max(double(spectrumValues(:)), 0);
    if isempty(spectrumValues) || ~isstruct(params)
        return;
    end

    photonTransportEnabled = isModelToggleEnabled(params, 'enablePhotonRefractionReflection', true) || ...
        isModelToggleEnabled(params, 'enableMQWInterfaceOptics', true) || ...
        isModelToggleEnabled(params, 'enableAbsorption', true);
    secondaryEnabled = isModelToggleEnabled(params, 'enableSecondaryExcitation', true);
    structuralParams = getStructuralParameterStruct(params);
    absorptionScale = clamp(getFieldOrDefault(structuralParams, 'absorptionScale', 1), 0, 5);
    secondaryScale = clamp(getFieldOrDefault(structuralParams, 'secondaryExcitationScale', 1), 0, 5);

    electronEnergy_keV = getFieldOrDefault(params, 'electronEnergy', 6e3) / 1e3;
    if ~isfinite(electronEnergy_keV) || electronEnergy_keV <= 0
        electronEnergy_keV = 6;
    end
    activeReach = clamp(1 ./ (1 + exp(-(electronEnergy_keV - 6.0) ./ 0.9)), 0, 1);

    localContext = getFieldOrDefault(params, 'pixelLocalContext', struct());
    vpitInfluence = clamp(getFieldOrDefault(localContext, 'vpitInfluence', 0), 0, 1);
    vpitMqwRecovery = clamp(getFieldOrDefault(localContext, 'vpitMqwRecovery', 0), 0, 0.8);

    mqwMainWindow = smoothBandpassNm(xi_nm, [500, 526], 8);
    mqwRedTailWindow = smoothBandpassNm(xi_nm, [532, 570], 10);
    shortwaveWindow = smoothBandpassNm(xi_nm, [350, 430], 18);

    if photonTransportEnabled
        transportStrength = clamp(0.55 + 0.035 * electronEnergy_keV + 0.16 * vpitInfluence + 0.10 * vpitMqwRecovery, 0.55, 1.15);
        interfaceEnvelope = 1 + ...
            transportStrength * (0.025 + 0.025 * activeReach) .* mqwMainWindow - ...
            transportStrength * 0.012 .* mqwRedTailWindow;
        if isModelToggleEnabled(params, 'enableAbsorption', true)
            shortwaveEscape = exp(-0.006 * absorptionScale * transportStrength .* shortwaveWindow);
            interfaceEnvelope = interfaceEnvelope .* shortwaveEscape;
        end
        spectrumValues = spectrumValues .* max(interfaceEnvelope, 0.82);
    end

    if secondaryEnabled && secondaryScale > 0
        shortArea = trapz(xi_nm(:), spectrumValues .* shortwaveWindow);
        if isfinite(shortArea) && shortArea > 0
            conversion = clamp((0.012 + 0.018 * activeReach + 0.012 * vpitInfluence) * secondaryScale, 0, 0.09);
            reEmissionProfile = normalizedGaussianProfileNm(xi_nm, 514, 22);
            spectrumValues = spectrumValues .* (1 - min(0.045, 0.012 * secondaryScale) .* shortwaveWindow) + ...
                shortArea * conversion .* reEmissionProfile;
        end
    end
end

function window = smoothBandpassNm(xi_nm, bandNm, edgeNm)
    xi_nm = double(xi_nm(:));
    edgeNm = max(double(edgeNm), eps);
    lower = 1 ./ (1 + exp(-(xi_nm - bandNm(1)) ./ edgeNm));
    upper = 1 ./ (1 + exp((xi_nm - bandNm(2)) ./ edgeNm));
    window = lower .* upper;
end

function profile = normalizedGaussianProfileNm(xi_nm, centerNm, fwhmNm)
    xi_nm = double(xi_nm(:));
    sigmaNm = max(double(fwhmNm), eps) / (2 * sqrt(2 * log(2)));
    profile = exp(-0.5 * ((xi_nm - centerNm) ./ sigmaNm) .^ 2);
    areaValue = trapz(xi_nm, profile);
    if areaValue > eps
        profile = profile / areaValue;
    end
end

function response = interpolateDetectorResponse(xi_nm, responseData)
    response = [];
    if istable(responseData) && width(responseData) >= 2
        responseAxis = double(responseData{:, 1});
        responseValues = double(responseData{:, 2});
    elseif isnumeric(responseData) && size(responseData, 2) >= 2
        responseAxis = double(responseData(:, 1));
        responseValues = double(responseData(:, 2));
    elseif isstruct(responseData) && isfield(responseData, 'wavelengthNm') && isfield(responseData, 'response')
        responseAxis = double(responseData.wavelengthNm(:));
        responseValues = double(responseData.response(:));
    else
        return;
    end
    validMask = isfinite(responseAxis) & isfinite(responseValues) & responseValues >= 0;
    responseAxis = responseAxis(validMask);
    responseValues = responseValues(validMask);
    if numel(responseAxis) < 2
        return;
    end
    [responseAxis, order] = sort(responseAxis(:));
    responseValues = responseValues(order);
    response = interp1(responseAxis, responseValues, xi_nm(:), 'pchip', 'extrap');
end

function response = buildDefaultDetectorResponse(xi_nm, shortwaveScale)
    if ~isfinite(shortwaveScale)
        shortwaveScale = 1;
    end
    shortwaveScale = clamp(shortwaveScale, 0.70, 1.30);
    shortwaveBand = exp(-0.5 * ((xi_nm(:) - 390) / 42) .^ 2);
    deepUvRollOff = 1 - 0.16 * exp(-0.5 * ((xi_nm(:) - 330) / 34) .^ 2);
    response = deepUvRollOff .* (1 + (shortwaveScale - 1) * shortwaveBand);
    response = clamp(response, 0.55, 1.35);
end

function spectrumValues = applyInstrumentResponseConvolution(rawSpectrum, sampleSpacing_nm, instrumentFwhm_nm)
    if isempty(rawSpectrum)
        spectrumValues = rawSpectrum;
        return;
    end

    sigmaBins = max(instrumentFwhm_nm / max(sampleSpacing_nm, eps) / 2.355, 0.8);
    kernelRadius = max(3, ceil(4 * sigmaBins));
    kernelAxis = (-kernelRadius:kernelRadius)';
    lsfKernel = exp(-0.5 * (kernelAxis ./ sigmaBins) .^ 2);
    lsfKernel = lsfKernel / max(sum(lsfKernel), eps);
    spectrumValues = conv(rawSpectrum, lsfKernel, 'same');
end

function instrumentFwhm_nm = estimateInstrumentResponseFwhmNm(params)
    instrumentFwhm_nm = 12.5;
    if ~isstruct(params)
        return;
    end

    instrumentFwhm_nm = getFieldOrDefault(params, 'instrumentResponseFwhm_nm', instrumentFwhm_nm);
    if isfield(params, 'calibrationProfile') && isstruct(params.calibrationProfile) && ~isempty(params.calibrationProfile)
        profileFwhm_nm = getFieldOrDefault(params.calibrationProfile, 'instrumentResponseFwhm_nm', NaN);
        if isfinite(profileFwhm_nm) && profileFwhm_nm > 0
            instrumentFwhm_nm = profileFwhm_nm;
        else
            componentsNm = getFieldOrDefault(params.calibrationProfile, 'componentsNm', []);
            if ~isempty(componentsNm)
                componentWidths_nm = double(componentsNm(:, 3));
                componentWidths_nm = componentWidths_nm(isfinite(componentWidths_nm) & componentWidths_nm > 0);
                if ~isempty(componentWidths_nm)
                    instrumentFwhm_nm = 8.2 + 0.22 * median(componentWidths_nm);
                end
            end
        end
    end
    instrumentFwhm_nm = clamp(instrumentFwhm_nm, 9.5, 18.0);
end

function finalBroadeningFwhm_nm = estimateFinalSpectrumBroadeningFwhmNm(params)
    finalBroadeningFwhm_nm = 0.35 * estimateInstrumentResponseFwhmNm(params);
    if ~isstruct(params) || ~isfield(params, 'pixelLocalContext') || isempty(params.pixelLocalContext)
        finalBroadeningFwhm_nm = clamp(finalBroadeningFwhm_nm, 3.0, 7.0);
        return;
    end

    localContext = params.pixelLocalContext;
    familyEntropy = computeNormalizedFamilyEntropy(getFieldOrDefault(localContext, 'targetFamilyFractions', []));
    if getFieldOrDefault(localContext, 'isVpitLikePixel', false)
        centerFillGain = clamp(getFieldOrDefault(localContext, 'centerFillGain', 1), 0.85, 1.35);
        finalBroadeningFwhm_nm = finalBroadeningFwhm_nm + 1.6 + 1.1 * familyEntropy + ...
            1.2 * max(centerFillGain - 1, 0);
    else
        finalBroadeningFwhm_nm = finalBroadeningFwhm_nm + 0.5 * familyEntropy;
    end
    finalBroadeningFwhm_nm = clamp(finalBroadeningFwhm_nm, 3.0, 10.5);
end

function adaptiveBroadeningFwhm_nm = estimateAdaptiveSpectrumBroadeningFwhmNm(weights, params)
    effectiveSampleCount = estimateEffectiveSpectrumSampleCount(weights);
    lowCountBroadening_nm = clamp(34 / sqrt(max(effectiveSampleCount, 1)) - 1.8, 0, 8.5);
    adaptiveBroadeningFwhm_nm = lowCountBroadening_nm;

    if ~isstruct(params) || ~isfield(params, 'pixelLocalContext') || isempty(params.pixelLocalContext)
        return;
    end

    localContext = params.pixelLocalContext;
    familyEntropy = computeNormalizedFamilyEntropy(getFieldOrDefault(localContext, 'targetFamilyFractions', []));
    if getFieldOrDefault(localContext, 'isVpitLikePixel', false)
        adaptiveBroadeningFwhm_nm = adaptiveBroadeningFwhm_nm + 1.2 + 1.1 * familyEntropy;
    else
        adaptiveBroadeningFwhm_nm = adaptiveBroadeningFwhm_nm + 0.5 * familyEntropy;
    end
    adaptiveBroadeningFwhm_nm = clamp(adaptiveBroadeningFwhm_nm, 0, 10.5);
end

function effectiveSampleCount = estimateEffectiveSpectrumSampleCount(weights)
    weights = double(weights(:));
    weights = weights(isfinite(weights) & weights > 0);
    if isempty(weights)
        effectiveSampleCount = 0;
        return;
    end
    effectiveSampleCount = (sum(weights) ^ 2) / max(sum(weights .^ 2), eps);
end

function [cross_boundary, boundary_normal, new_layer, boundary_data] = checkBoundaryCrossing(old_pos, new_pos, layerBoundaries, layers, Vpits, nominalSurfaceZ, params)
    if nargin < 7
        params = struct();
    end
    num_rays = size(old_pos, 1);
    boundary_normal = zeros(num_rays, 3);
    boundary_data = repmat(struct( ...
        'crossed', false, ...
        'point_m', [0, 0, 0], ...
        'normal', [0, 0, 1], ...
        'boundaryType', "", ...
        'oldState', struct(), ...
        'newState', struct(), ...
        'newLayer', 0, ...
        't', 1), num_rays, 1);

    if isModelToggleEnabled(params, 'enableVpitTransportCorrection', true)
        sensitiveMask = classifyVpitSensitiveSegments(old_pos, new_pos, Vpits, nominalSurfaceZ);
    else
        sensitiveMask = false(num_rays, 1);
    end
    bulkMask = ~sensitiveMask;

    if any(bulkMask)
        [bulkCross, bulkNormal, bulkNewLayer, bulkBoundaryData] = ...
            checkPlanarBoundaryCrossingFast(old_pos(bulkMask, :), new_pos(bulkMask, :), layerBoundaries, layers);
        cross_boundary(bulkMask) = bulkCross;
        boundary_normal(bulkMask, :) = bulkNormal;
        new_layer(bulkMask) = bulkNewLayer;
        boundary_data(bulkMask) = bulkBoundaryData;
    end

    sensitiveIdx = find(sensitiveMask);
    for localIdx = 1:numel(sensitiveIdx)
        i = sensitiveIdx(localIdx);
        crossing = precies.vpit('firstBoundaryCrossing', old_pos(i, :), new_pos(i, :), layers, layerBoundaries, Vpits, nominalSurfaceZ);
        boundary_data(i) = crossing;
        cross_boundary(i) = crossing.crossed;
        boundary_normal(i, :) = crossing.normal;
        new_layer(i) = crossing.newLayer;
    end
end

function [cross_boundary, boundary_normal, new_layer, boundary_data] = checkPlanarBoundaryCrossingFast(old_pos, new_pos, layerBoundaries, layers)
    num_rays = size(old_pos, 1);
    boundary_normal = zeros(num_rays, 3);
    boundary_data = repmat(struct( ...
        'crossed', false, ...
        'point_m', [0, 0, 0], ...
        'normal', [0, 0, 1], ...
        'boundaryType', "layer_plane", ...
        'oldState', struct(), ...
        'newState', struct(), ...
        'newLayer', 0, ...
        't', 1), num_rays, 1);

    old_layer = findLayerIndicesFast(old_pos(:, 3), layerBoundaries, size(layers, 1));
    new_layer = findLayerIndicesFast(new_pos(:, 3), layerBoundaries, size(layers, 1));
    cross_boundary = old_layer ~= new_layer;

    crossingIdx = find(cross_boundary);
    for idx = reshape(crossingIdx, 1, [])
        boundaryZ = inferPlanarBoundaryZ(old_pos(idx, 3), new_pos(idx, 3), layerBoundaries);
        t = 1;
        if abs(new_pos(idx, 3) - old_pos(idx, 3)) > eps
            t = (boundaryZ - old_pos(idx, 3)) / (new_pos(idx, 3) - old_pos(idx, 3));
            t = max(0, min(1, t));
        end
        point = old_pos(idx, :) + t * (new_pos(idx, :) - old_pos(idx, :));
        if new_pos(idx, 3) > old_pos(idx, 3)
            normal = [0, 0, 1];
        else
            normal = [0, 0, -1];
        end

        boundary_normal(idx, :) = normal;
        boundary_data(idx) = struct( ...
            'crossed', true, ...
            'point_m', point, ...
            'normal', normal, ...
            'boundaryType', "layer_plane", ...
            'oldState', buildBulkStateFromLayerIndex(old_layer(idx), layers), ...
            'newState', buildBulkStateFromLayerIndex(new_layer(idx), layers), ...
            'newLayer', new_layer(idx), ...
            't', t);
    end
end

function [reflected_dir, refracted_dir, reflection_fraction, refraction_fraction, skip_refraction] = ...
    handleBoundaryCrossing(dir, lambda, pol, boundary_data, layers, params)
    if nargin < 6
        params = struct();
    end
    num_rays = size(dir, 1);
    reflected_dir = dir;
    refracted_dir = dir;
    reflection_fraction = zeros(num_rays, 1);
    refraction_fraction = zeros(num_rays, 1);
    skip_refraction = false(num_rays, 1);
    
    for i = 1:num_rays
        crossing = boundary_data(i);
        if ~crossing.crossed
            continue;
        end

        normal = crossing.normal;
        normalNorm = norm(normal);
        if normalNorm <= eps
            continue;
        end
        normal = normal / normalNorm;

        if ~isModelToggleEnabled(params, 'enablePhotonRefractionReflection', true)
            skip_refraction(i) = true;
            reflection_fraction(i) = 0;
            refraction_fraction(i) = 1;
            refracted_dir(i, :) = dir(i, :);
            continue;
        end

        if shouldSkipMQWRefraction(crossing, layers, params)
            skip_refraction(i) = true;
            reflection_fraction(i) = 0;
            refraction_fraction(i) = 1;
            refracted_dir(i, :) = dir(i, :);
            continue;
        end

        n1_func = getStateRefractiveIndexFunc(crossing.oldState);
        n2_func = getStateRefractiveIndexFunc(crossing.newState);

        cos_theta1 = abs(max(min(dot(dir(i, :), normal), 1), -1));
        R = fresnelReflectance(n1_func, n2_func, lambda(i), cos_theta1, pol(i));

        reflection_fraction(i) = R;
        refraction_fraction(i) = max(0, 1 - R);
        reflected_dir(i, :) = dir(i, :) - 2 * dot(dir(i, :), normal) * normal;
        reflected_dir(i, :) = reflected_dir(i, :) / max(norm(reflected_dir(i, :)), eps);

        if refraction_fraction(i) > 0
            [candidateRefractedDir, hasTotalInternalReflection] = calculateRefractedDirection( ...
                dir(i, :), normal, n1_func(lambda(i)), n2_func(lambda(i)));
            if hasTotalInternalReflection
                refraction_fraction(i) = 0;
                reflection_fraction(i) = 1;
            else
                refracted_dir(i, :) = candidateRefractedDir;
            end
        end
    end
end

function mqw_transition = isMQWTransition(old_layer, new_layer, layers)
   
    mqw_transition = false(size(old_layer));
    
    for i = 1:length(old_layer)
        if old_layer(i) >= 1 && old_layer(i) <= size(layers, 1) && ...
           new_layer(i) >= 1 && new_layer(i) <= size(layers, 1)
            
            old_name = layers{old_layer(i), 1};
            new_name = layers{new_layer(i), 1};
     
            if contains(old_name, 'MQW') && contains(new_name, 'MQW')
                mqw_transition(i) = true;
            end
        end
    end
end

function tf = shouldSkipMQWRefraction(crossing, layers, params)
    if nargin < 3
        params = struct();
    end
    if isModelToggleEnabled(params, 'enableMQWInterfaceOptics', true)
        tf = false;
        return;
    end
    tf = crossing.boundaryType == "layer_plane" && ...
        strcmp(crossing.oldState.materialKind, 'bulk') && ...
        strcmp(crossing.newState.materialKind, 'bulk') && ...
        isMQWTransition(crossing.oldState.layerIndex, crossing.newState.layerIndex, layers);
end

function nFunc = getStateRefractiveIndexFunc(state)
    if isfield(state, 'nFunc') && ~isempty(state.nFunc)
        nFunc = state.nFunc;
    else
        nFunc = @(~) 1.0;
    end
end

function [refractedDir, hasTotalInternalReflection] = calculateRefractedDirection(dir, normal, n1, n2)
    eta = n1 / max(n2, eps);
    cosTheta1 = abs(max(min(dot(dir, normal), 1), -1));
    k = 1 - eta^2 * (1 - cosTheta1^2);

    if k <= 0
        refractedDir = dir;
        hasTotalInternalReflection = true;
        return;
    end

    refractedDir = eta * dir + (eta * cosTheta1 - sqrt(k)) * normal;
    refractedDir = refractedDir / max(norm(refractedDir), eps);
    hasTotalInternalReflection = false;
end

function [absorbed, PL_intensity, PL_pos, PL_wavelength, PL_data] = ...
    handleAbsorptionAndPL(pos, dir, intensity, lambda, layer, layers, layerBoundaries, stepSize, ...
    Eg_GaN, Eg_InN, ~, params, Vpits, nominalSurfaceZ)
    absorbed = zeros(size(pos, 1), 1);
    PL_intensity = zeros(size(pos, 1), 1);
    PL_pos = zeros(size(pos));
    PL_wavelength = zeros(size(pos, 1), 1);
    PL_data = struct('positions', [], 'wavelengths', [], 'intensities', []);

    if ~isModelToggleEnabled(params, 'enableAbsorption', true)
        return;
    end

    absorptionScale = getStructuralOrModelParameter(params, 'absorptionScale', 1);
    absorptionScale = clamp(absorptionScale, 0, 5);
    secondaryScale = getStructuralOrModelParameter(params, 'secondaryExcitationScale', 1);
    if ~isModelToggleEnabled(params, 'enableSecondaryExcitation', true)
        secondaryScale = 0;
    end
    secondaryScale = clamp(secondaryScale, 0, 5);

    samplePos_nm = (pos + 0.5 * stepSize * dir) * 1e9;
    if isModelToggleEnabled(params, 'enableVpitTransportCorrection', true)
        sensitiveMask = classifyVpitSensitivePoints(samplePos_nm, Vpits, nominalSurfaceZ);
    else
        sensitiveMask = false(size(pos, 1), 1);
    end
    bulkMask = ~sensitiveMask;

    bulkIdx = find(bulkMask);
    bulkLayers = layer(bulkIdx);
    uniqueLayers = unique(bulkLayers(:)');
    for layerIdx = uniqueLayers
        if layerIdx < 1 || layerIdx > size(layers, 1) || isempty(layers{layerIdx, 6})
            continue;
        end

        layerMask = bulkLayers == layerIdx;
        layerBulkIdx = bulkIdx(layerMask);
        alpha = absorptionScale * max(0, layers{layerIdx, 6}(lambda(layerBulkIdx)));
        absorption_prob = 1 - exp(-alpha * stepSize);
        absorptionHits = rand(numel(layerBulkIdx), 1) < absorption_prob(:);
        if ~any(absorptionHits)
            continue;
        end

        absorbedIdx = layerBulkIdx(absorptionHits);
        absorbedProb = absorption_prob(absorptionHits);
        absorbed(absorbedIdx) = intensity(absorbedIdx) .* absorbedProb;

        bandgap = get_bandgap(layerIdx, layers, Eg_GaN, Eg_InN);
        photon_energy = 1240e-9 ./ lambda(absorbedIdx);
        plMask = photon_energy >= bandgap;
        if secondaryScale <= 0 || ~any(plMask)
            continue;
        end

        state = buildBulkStateFromLayerIndex(layerIdx, layers);
        state.bandgap_eV = bandgap;
        plIdx = absorbedIdx(plMask);
        PL_intensity(plIdx) = absorbed(plIdx) * resolveSecondaryEmissionYield(state, layers) * secondaryScale;
        PL_pos(plIdx, :) = pos(plIdx, :);
        for wavelengthIdx = reshape(plIdx, 1, [])
            PL_wavelength(wavelengthIdx) = sampleStateEmissionWavelength(state, layers, params, Eg_GaN, Eg_InN);
        end
        PL_data.positions = [PL_data.positions; pos(plIdx, :)];
        PL_data.wavelengths = [PL_data.wavelengths; PL_wavelength(plIdx)];
        PL_data.intensities = [PL_data.intensities; PL_intensity(plIdx)];
    end

    sensitiveIdx = find(sensitiveMask);
    for localIdx = 1:numel(sensitiveIdx)
        i = sensitiveIdx(localIdx);
        state = precies.vpit('resolvePointState', samplePos_nm(i, :), layers, layerBoundaries, Vpits, nominalSurfaceZ);

        if strcmp(state.materialKind, 'cavity') || strcmp(state.materialKind, 'air') || isempty(state.alphaFunc)
            continue;
        end

        alpha = absorptionScale * max(0, state.alphaFunc(lambda(i)));
        absorption_prob = 1 - exp(-alpha * stepSize);

        if rand() < absorption_prob
            absorbed(i) = intensity(i) * absorption_prob;

            if isfinite(state.bandgap_eV)
                bandgap = state.bandgap_eV;
            else
                bandgap = get_bandgap(layer(i), layers, Eg_GaN, Eg_InN);
            end

            photon_energy = 1240e-9 / lambda(i);
            if secondaryScale > 0 && photon_energy >= bandgap
                PL_intensity(i) = absorbed(i) * resolveSecondaryEmissionYield(state, layers) * secondaryScale;
                PL_pos(i, :) = pos(i, :);
                PL_wavelength(i) = sampleStateEmissionWavelength(state, layers, params, Eg_GaN, Eg_InN);
                PL_data.positions = [PL_data.positions; pos(i, :)];
                PL_data.wavelengths = [PL_data.wavelengths; PL_wavelength(i)];
                PL_data.intensities = [PL_data.intensities; PL_intensity(i)];
            end
        end
    end
end

function yieldValue = resolveSecondaryEmissionYield(state, layers)
    yieldValue = 0.20;
    if isstruct(state)
        if isfield(state, 'materialKind') && strcmp(state.materialKind, 'semipolar_shell')
            yieldValue = 0.28;
        end
        layerName = '';
        if isfield(state, 'layerName') && ~isempty(state.layerName)
            layerName = char(string(state.layerName));
        elseif isfield(state, 'layerIndex') && state.layerIndex >= 1 && state.layerIndex <= size(layers, 1)
            layerName = char(string(layers{state.layerIndex, 1}));
        end
        if contains(layerName, 'MQW-Well')
            yieldValue = max(yieldValue, 0.30);
        elseif contains(layerName, 'MQW-Barrier')
            yieldValue = max(yieldValue, 0.24);
        elseif contains(layerName, 'GaN')
            yieldValue = max(yieldValue, 0.18);
        end
    end
    yieldValue = clamp(yieldValue, 0.08, 0.42);
end

function lambda = sampleStateEmissionWavelength(state, layers, params, Eg_GaN, Eg_InN)
    sourceRow = [0, 0, 0, state.bandgap_eV, 0, classifyCarrierEmissionSource(state), 0];
    spectrumModel = resolveEmissionSpectrumModel(sourceRow, state.layerIndex, layers, params);
    lambda = sampleEmissionWavelength(spectrumModel);

    if isfinite(state.bandgap_eV)
        referenceBandgap = get_bandgap(state.layerIndex, layers, Eg_GaN, Eg_InN);
        photonEnergy = 1240e-9 / lambda + max(0, state.bandgap_eV - referenceBandgap);
        lambda = 1240e-9 / max(photonEnergy, 0.05);
    end
end

function terminated = checkTerminationConditions(pos, intensity, depth, layerBoundaries, initialIntensity, maxDepth)
    terminated = pos(:, 3) < layerBoundaries(1) | ...
        pos(:, 3) >= layerBoundaries(end-1) | ...
        intensity < 0.0001 .* initialIntensity | ...
        depth > maxDepth;
end

function new_rays = addRefractedRays(pos, dir, intensity, pol, depth, layer, lambda, isPL, isPLed, initialIntensity, existingBatch)
    if nargin < 11 || isempty(existingBatch)
        existingBatch = struct( ...
            'dir', zeros(0, 3), ...
            'pos', zeros(0, 3), ...
            'intensity', zeros(0, 1), ...
            'pol', false(0, 1), ...
            'depth', zeros(0, 1), ...
            'layer', zeros(0, 1), ...
            'lambda', zeros(0, 1), ...
            'isPL', false(0, 1), ...
            'isPLed', false(0, 1), ...
            'initialIntensity', zeros(0, 1));
    end

    new_rays = existingBatch;
    new_rays.dir = [new_rays.dir; dir];
    new_rays.pos = [new_rays.pos; pos];
    new_rays.intensity = [new_rays.intensity; intensity(:)];
    new_rays.pol = [new_rays.pol; logical(pol(:))];
    new_rays.depth = [new_rays.depth; depth(:) + 1];
    new_rays.layer = [new_rays.layer; layer(:)];
    new_rays.lambda = [new_rays.lambda; lambda(:)];
    new_rays.isPL = [new_rays.isPL; logical(isPL(:))];
    new_rays.isPLed = [new_rays.isPLed; logical(isPLed(:))];
    new_rays.initialIntensity = [new_rays.initialIntensity; initialIntensity(:)];
end

function new_rays = addPLRays(pos, intensity, wavelength, depth, layer, initialIntensity)
    num_rays = size(pos, 1);
    theta = acos(sqrt(rand(num_rays, 1)));
    phi = 2 * pi * rand(num_rays, 1);
    dir = [sin(theta) .* cos(phi), sin(theta) .* sin(phi), abs(cos(theta))];
    pol = rand(num_rays, 1) > 0.5;
    isPL = true(num_rays, 1);
    isPLed = true(num_rays, 1);
    
    new_rays = struct(...
        'dir', dir, ...
        'pos', pos, ...
        'intensity', intensity, ...
        'pol', pol, ...
        'depth', depth + 1, ...
        'layer', layer, ...
        'lambda', wavelength, ...
        'isPL', isPL, ...
        'isPLed', isPLed, ...
        'initialIntensity', initialIntensity ...
    );
end

function [ray_dir, ray_pos, ray_intensity, ray_pol, ray_depth, ray_layer, ray_lambda, ...
    ray_isPL, ray_isPLed, ray_initialIntensity, active_rays] = appendNewRays( ...
    ray_dir, ray_pos, ray_intensity, ray_pol, ray_depth, ray_layer, ray_lambda, ...
    ray_isPL, ray_isPLed, ray_initialIntensity, active_rays, new_rays)

    if isempty(new_rays) || ~isstruct(new_rays) || isempty(new_rays.pos)
        return;
    end

    ray_dir = [ray_dir; castNumericState(new_rays.dir, ray_dir)];
    ray_pos = [ray_pos; castNumericState(new_rays.pos, ray_pos)];
    ray_intensity = [ray_intensity; castNumericState(new_rays.intensity(:), ray_intensity)];
    ray_pol = [ray_pol; logical(new_rays.pol(:))];
    ray_depth = [ray_depth; castNumericState(new_rays.depth(:), ray_depth)];
    ray_layer = [ray_layer; new_rays.layer(:)];
    ray_lambda = [ray_lambda; castNumericState(new_rays.lambda(:), ray_lambda)];
    ray_isPL = [ray_isPL; castLogicalState(new_rays.isPL(:), ray_isPL)];
    ray_isPLed = [ray_isPLed; castLogicalState(new_rays.isPLed(:), ray_isPLed)];
    ray_initialIntensity = [ray_initialIntensity; castNumericState(new_rays.initialIntensity(:), ray_initialIntensity)];
    active_rays = [active_rays; castLogicalState(true(size(new_rays.pos, 1), 1), active_rays)];
end
