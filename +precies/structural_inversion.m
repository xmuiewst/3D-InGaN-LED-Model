function varargout = structural_inversion(action, varargin)

    switch lower(string(action))
        case "run"
            varargout{1} = runStructuralInversion(varargin{:});
        case "ablation"
            varargout{1} = runAblation(varargin{:});
        case "forwardablation"
            varargout{1} = runForwardAblation(varargin{:});
        case "forwardfit"
            varargout{1} = runForwardFit(varargin{:});
        case "leaveonevoltageout"
            varargout{1} = runLeaveOneVoltageOut(varargin{:});
        case "sensitivity"
            varargout{1} = runSensitivity(varargin{:});
        case "uncertainty"
            varargout{1} = runUncertainty(varargin{:});
        case "depthresponse"
            varargout{1} = estimateDepthResponse(varargin{:});
        case "defaultoptions"
            varargout{1} = defaultOptions();
        case "selftest"
            varargout{1} = selfTest();
        otherwise
            error('precies:structural_inversion:InvalidAction', ...
                'Unsupported action: %s', action);
    end
end

function result = runStructuralInversion(experimentalData, baselineData, options)
    if nargin < 3
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);
    [expSeries, baseSeries] = normalizePairedSeries(experimentalData, baselineData, opts);
    bounds = buildBounds(opts);
    initialVector = initialParameterVector(opts, bounds);

    initialScore = objectiveForVector(initialVector, bounds, expSeries, baseSeries, opts, struct());
    if opts.enableOptimization
        optimOptions = optimset( ...
            'Display', opts.optimizerDisplay, ...
            'MaxIter', opts.maxIterations, ...
            'MaxFunEvals', opts.maxFunctionEvaluations, ...
            'TolX', opts.toleranceX, ...
            'TolFun', opts.toleranceFun);
        objective = @(u) objectiveForVector(unconstrainedToParams(u, bounds), ...
            bounds, expSeries, baseSeries, opts, struct());
        [uBest, bestScore, exitFlag, optimOutput] = fminsearch( ...
            objective, paramsToUnconstrained(initialVector, bounds), optimOptions);
        bestVector = unconstrainedToParams(uBest, bounds);
        iterations = getFieldOrDefault(optimOutput, 'iterations', 0);
    else
        bestVector = initialVector;
        bestScore = initialScore;
        exitFlag = 1;
        iterations = 0;
    end

    bestParams = vectorToParameterStruct(bestVector, opts.parameterNames);
    fittedSeries = synthesizeSeries(baseSeries, bestParams, opts, struct());
    initialSeries = synthesizeSeries(baseSeries, vectorToParameterStruct(initialVector, opts.parameterNames), opts, struct());
    initialMetrics = evaluateSeriesMetrics(expSeries, initialSeries, opts);
    fittedMetrics = evaluateSeriesMetrics(expSeries, fittedSeries, opts);

    result = struct();
    result.fixedInputPriors = opts.fixedInputPriors;
    result.independentValidationReferences = opts.independentValidationReferences;
    result.parameterNames = opts.parameterNames;
    result.bounds = bounds;
    result.initialParameters = vectorToParameterStruct(initialVector, opts.parameterNames);
    result.bestParameters = bestParams;
    result.initialObjective = initialScore;
    result.objective = bestScore;
    result.exitFlag = exitFlag;
    result.iterations = iterations;
    result.initialMetrics = initialMetrics;
    result.metrics = fittedMetrics;
    result.fittedSeries = fittedSeries;
    result.parameterTable = buildParameterTable(result.bestParameters, bounds, opts);
    result.independentValidation = buildIndependentValidationTable(result.bestParameters, opts.independentValidationReferences);
    result.depthResponse = estimateDepthResponse(result.bestParameters, opts);
    result.sideViewCLValidation = buildSideViewCLValidation(result.depthResponse, opts.independentValidationReferences);
    result.sensitivity = runSensitivity(expSeries, baseSeries, result, opts);
    result.uncertainty = runUncertainty(expSeries, baseSeries, result, opts);
    result.options = opts;

    if opts.writeOutputs && ~isempty(opts.outputDirectory)
        writeOutputs(result, opts);
    end
end

function result = runAblation(experimentalData, baselineData, structuralResult, options)
    if nargin < 4
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);
    [expSeries, baseSeries] = normalizePairedSeries(experimentalData, baselineData, opts);
    bestParams = resolveBestParameters(structuralResult, opts);

    scenarioNames = {'M0', 'M1', 'M2', 'M3'};
    scenarios = { ...
        struct('disablePEInteractionVolume', true, 'disableRefractionReflection', true, ...
        'disableAbsorption', true, 'disableSecondaryExcitation', true), ...
        struct('disableRefractionReflection', true, 'disableAbsorption', true, ...
        'disableSecondaryExcitation', true), ...
        struct('disableSecondaryExcitation', true), ...
        struct()};

    metricsCell = cell(numel(scenarios), 1);
    objectives = zeros(numel(scenarios), 1);
    for idx = 1:numel(scenarios)
        fittedSeries = synthesizeSeries(baseSeries, bestParams, opts, scenarios{idx});
        metricsCell{idx} = evaluateSeriesMetrics(expSeries, fittedSeries, opts);
        objectives(idx) = objectiveFromMetrics(metricsCell{idx}, opts);
    end

    result = struct();
    result.scenarioNames = scenarioNames(:);
    result.metrics = vertcat(metricsCell{:});
    result.objective = objectives;
    result.summaryTable = buildAblationTable(result);
    result.note = ['Use precomputed forward simulations with the same scenario labels when available. ' ...
        'The built-in fallback is a structural-response surrogate for fast screening.'];
end

function result = runForwardAblation(experimentalData, simulationOptions, structuralResult, options)
    if nargin < 4
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);
    bestParams = resolveBestParameters(structuralResult, opts);
    scenarios = buildForwardAblationScenarios();

    scenarioNames = string({scenarios.name})';
    metricsCell = cell(numel(scenarios), 1);
    objectives = zeros(numel(scenarios), 1);
    simulationDataCell = cell(numel(scenarios), 1);

    for idx = 1:numel(scenarios)
        scenarioOptions = applyForwardAblationScenario(simulationOptions, bestParams, scenarios(idx));
        if opts.progressEveryScenario > 0
            fprintf('forward ablation: %s (%d/%d)\n', scenarios(idx).name, idx, numel(scenarios));
        end
        simData = precies.simulation('simulateOptions', scenarioOptions);
        metricExpData = cropDataForMetrics(experimentalData, opts.roi);
        metricSimData = cropDataForMetrics(simData, opts.roi);
        metricOpts = struct( ...
            'ganWindow', opts.ganWindow, ...
            'mqwWindow', opts.mqwWindow, ...
            'wavelengthAxis', opts.wavelengthAxis, ...
            'voltage', inferForwardAblationVoltage(scenarioOptions), ...
            'enableOptimization', false);
        metricsCell{idx} = precies.pixelwise_inversion('metrics', metricExpData, metricSimData, metricOpts);
        objectives(idx) = objectiveFromMetrics(metricsCell{idx}, opts);
        if opts.storeForwardAblationData
            simulationDataCell{idx} = simData;
        end
    end

    metrics = vertcat(metricsCell{:});
    result = struct();
    result.scenarioNames = scenarioNames;
    result.metrics = metrics;
    result.objective = objectives;
    result.summaryTable = buildForwardAblationTable(scenarioNames, metrics, objectives);
    result.simulationData = simulationDataCell;
    result.bestParameters = bestParams;
    result.roi = opts.roi;
    result.note = ['Forward ablation actually reruns the forward Monte Carlo/ray-tracing model for each model toggle. ' ...
        'Use sufficient rays/electrons for publication-quality numbers.'];
end

function result = runForwardFit(experimentalData, simulationOptionsSeries, structuralSeed, options)
    if nargin < 3
        structuralSeed = struct();
    end
    if nargin < 4
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);
    expSeries = normalizeSeries(experimentalData, opts);
    simOptionsCell = normalizeForwardSimulationOptionsSeries(simulationOptionsSeries);
    nPairs = min(numel(expSeries), numel(simOptionsCell));
    if nPairs < 1
        result = emptyForwardFitResult(opts, 'No simulation options were supplied.');
        return;
    end

    seedParams = resolveBestParameters(structuralSeed, opts);
    bounds = buildBounds(opts);
    seedVector = parameterStructToVector(seedParams, opts.parameterNames);
    seedVector = min(bounds(:, 2), max(bounds(:, 1), seedVector));
    [candidateVectors, candidateLabel] = buildForwardFitCandidates(seedVector, bounds, opts);

    nCandidates = size(candidateVectors, 1);
    candidateObjective = nan(nCandidates, 1);
    candidateNcc = nan(nCandidates, 1);
    candidateNrmse = nan(nCandidates, 1);
    candidateResidualRms = nan(nCandidates, 1);
    candidateMetrics = cell(nCandidates, 1);
    candidateSimData = cell(nCandidates, nPairs);

    for candidateIdx = 1:nCandidates
        params = vectorToParameterStruct(candidateVectors(candidateIdx, :)', opts.parameterNames);
        if opts.progressEveryScenario > 0
            fprintf('forward structural fit: %s (%d/%d)\n', candidateLabel(candidateIdx), candidateIdx, nCandidates);
        end
        [metrics, objective, simDataCell] = evaluateForwardFitCandidate( ...
            expSeries(1:nPairs), simOptionsCell(1:nPairs), params, opts);
        candidateMetrics{candidateIdx} = metrics;
        candidateObjective(candidateIdx) = objective;
        candidateNcc(candidateIdx) = mean(metrics.summaryTable.Average_NCC, 'omitnan');
        candidateNrmse(candidateIdx) = mean(metrics.summaryTable.nRMSE, 'omitnan');
        candidateResidualRms(candidateIdx) = mean(metrics.summaryTable.MeanResidualRMS, 'omitnan');
        if opts.storeForwardFitData
            candidateSimData(candidateIdx, :) = simDataCell(:)';
        end
    end

    [bestObjective, bestIdx] = min(candidateObjective);
    if isempty(bestIdx) || ~isfinite(bestObjective)
        bestIdx = 1;
    end
    bestVector = candidateVectors(bestIdx, :)';
    bestParams = vectorToParameterStruct(bestVector, opts.parameterNames);

    result = struct();
    result.method = 'bounded full-forward candidate search';
    result.note = ['This fit reruns simulateOptions for a finite set of global structural-parameter candidates. ' ...
        'It is stronger than the surrogate inversion, but still not a dense posterior or continuous full MC optimizer.'];
    result.seedParameters = vectorToParameterStruct(seedVector, opts.parameterNames);
    result.bestParameters = bestParams;
    result.objective = candidateObjective(bestIdx);
    result.metrics = candidateMetrics{bestIdx};
    result.candidateTable = buildForwardFitCandidateTable( ...
        candidateLabel, candidateVectors, candidateObjective, candidateNcc, candidateNrmse, candidateResidualRms, opts);
    result.simulationData = candidateSimData;
    result.parameterTable = buildParameterTable(bestParams, bounds, opts);
    result.depthResponse = estimateDepthResponse(bestParams, opts);
end

function result = runLeaveOneVoltageOut(experimentalData, baselineData, options)
    if nargin < 3
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);
    [expSeries, baseSeries] = normalizePairedSeries(experimentalData, baselineData, opts);
    n = numel(expSeries);
    heldVoltage = zeros(n, 1);
    trainObjective = zeros(n, 1);
    heldObjective = zeros(n, 1);
    avgNCC = zeros(n, 1);
    avgNRMSE = zeros(n, 1);
    trained = repmat(struct(), n, 1);

    looOpts = opts;
    looOpts.maxIterations = min(opts.maxIterations, opts.leaveOneOutMaxIterations);
    looOpts.maxFunctionEvaluations = min(opts.maxFunctionEvaluations, opts.leaveOneOutMaxFunctionEvaluations);
    looOpts.uncertaintySamples = 0;
    looOpts.sensitivityStepFraction = 0;

    for holdIdx = 1:n
        trainMask = true(1, n);
        trainMask(holdIdx) = false;
        trainResult = runStructuralInversion({expSeries(trainMask).data}, {baseSeries(trainMask).data}, looOpts);
        heldFit = synthesizeSeries(baseSeries(holdIdx), trainResult.bestParameters, opts, struct());
        heldMetrics = evaluateSeriesMetrics(expSeries(holdIdx), heldFit, opts);
        trained(holdIdx).bestParameters = trainResult.bestParameters;
        heldVoltage(holdIdx) = expSeries(holdIdx).voltage;
        trainObjective(holdIdx) = trainResult.objective;
        heldObjective(holdIdx) = objectiveFromMetrics(heldMetrics, opts);
        avgNCC(holdIdx) = heldMetrics.summaryTable.Average_NCC(1);
        avgNRMSE(holdIdx) = heldMetrics.summaryTable.nRMSE(1);
    end

    result = struct();
    result.heldVoltage = heldVoltage;
    result.trainObjective = trainObjective;
    result.heldObjective = heldObjective;
    result.averageNCC = avgNCC;
    result.averageNRMSE = avgNRMSE;
    result.trained = trained;
    result.summaryTable = table(heldVoltage, trainObjective, heldObjective, avgNCC, avgNRMSE, ...
        'VariableNames', {'HeldVoltage_kV', 'TrainObjective', 'HeldObjective', 'Average_NCC', 'nRMSE'});
end

function result = runSensitivity(experimentalData, baselineData, structuralResult, options)
    if nargin < 4
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);
    if opts.sensitivityStepFraction <= 0
        result = emptySensitivityResult(opts);
        return;
    end
    [expSeries, baseSeries] = normalizePairedSeries(experimentalData, baselineData, opts);
    bestParams = resolveBestParameters(structuralResult, opts);
    bounds = buildBounds(opts);
    bestVector = parameterStructToVector(bestParams, opts.parameterNames);
    baseObjective = objectiveForVector(bestVector, bounds, expSeries, baseSeries, opts, struct());

    n = numel(opts.parameterNames);
    lowObjective = zeros(n, 1);
    highObjective = zeros(n, 1);
    sensitivityScore = zeros(n, 1);
    for idx = 1:n
        span = bounds(idx, 2) - bounds(idx, 1);
        delta = opts.sensitivityStepFraction * span;
        lowVector = bestVector;
        highVector = bestVector;
        lowVector(idx) = max(bounds(idx, 1), bestVector(idx) - delta);
        highVector(idx) = min(bounds(idx, 2), bestVector(idx) + delta);
        lowObjective(idx) = objectiveForVector(lowVector, bounds, expSeries, baseSeries, opts, struct());
        highObjective(idx) = objectiveForVector(highVector, bounds, expSeries, baseSeries, opts, struct());
        sensitivityScore(idx) = max(abs([lowObjective(idx), highObjective(idx)] - baseObjective)) / max(delta, eps);
    end

    parameter = string(opts.parameterNames(:));
    bestValue = bestVector(:);
    result = struct();
    result.baseObjective = baseObjective;
    result.table = table(parameter, bestValue, lowObjective, highObjective, sensitivityScore, ...
        'VariableNames', {'Parameter', 'BestValue', 'LowerPerturbObjective', 'UpperPerturbObjective', 'SensitivityScore'});
end

function result = runUncertainty(experimentalData, baselineData, structuralResult, options)
    if nargin < 4
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);
    if opts.uncertaintySamples <= 0
        result = emptyUncertaintyResult(opts);
        return;
    end
    [expSeries, baseSeries] = normalizePairedSeries(experimentalData, baselineData, opts);
    bestParams = resolveBestParameters(structuralResult, opts);
    bounds = buildBounds(opts);
    bestVector = parameterStructToVector(bestParams, opts.parameterNames);
    nParams = numel(bestVector);
    span = bounds(:, 2) - bounds(:, 1);
    samples = zeros(opts.uncertaintySamples, nParams);
    scores = zeros(opts.uncertaintySamples, 1);
    samples(1, :) = bestVector(:)';
    scores(1) = objectiveForVector(bestVector, bounds, expSeries, baseSeries, opts, struct());

    for sampleIdx = 2:opts.uncertaintySamples
        candidate = bestVector(:) + opts.uncertaintyStepFraction * span(:) .* randn(nParams, 1);
        candidate = min(bounds(:, 2), max(bounds(:, 1), candidate));
        samples(sampleIdx, :) = candidate(:)';
        scores(sampleIdx) = objectiveForVector(candidate, bounds, expSeries, baseSeries, opts, struct());
    end

    scoreThreshold = min(scores) + opts.uncertaintyObjectiveWindow * max(std(scores, 0, 'omitnan'), eps);
    accepted = scores <= scoreThreshold;
    if nnz(accepted) < 3
        [~, order] = sort(scores, 'ascend');
        accepted(order(1:min(3, numel(order)))) = true;
    end

    acceptedSamples = samples(accepted, :);
    ciLow = prctile(acceptedSamples, 16, 1)';
    ciHigh = prctile(acceptedSamples, 84, 1)';
    parameter = string(opts.parameterNames(:));
    bestValue = bestVector(:);
    result = struct();
    result.samples = samples;
    result.objective = scores;
    result.acceptedMask = accepted;
    result.parameterCI = table(parameter, bestValue, ciLow, ciHigh, ...
        'VariableNames', {'Parameter', 'BestValue', 'CI16', 'CI84'});
end

function depthResponse = estimateDepthResponse(parameters, options)
    if nargin < 2
        options = struct();
    end
    opts = mergeOptions(defaultOptions(), options);
    params = resolveBestParameters(parameters, opts);
    voltages = opts.voltagesKeV(:);
    if isempty(voltages) || any(~isfinite(voltages))
        voltages = (4:10)';
    end
    depthGridNm = (opts.depthResponseRangeNm(1):opts.depthResponseStepNm:opts.depthResponseRangeNm(2))';
    responses = zeros(numel(depthGridNm), numel(voltages));
    fwhmNm = zeros(numel(voltages), 1);
    centroidNm = zeros(numel(voltages), 1);

    diffusionNm = getFieldOrDefault(params, 'carrierDiffusionLengthNm', 60);
    activeThicknessNm = getFieldOrDefault(params, 'mqwActiveThicknessNm', 150) + ...
        getFieldOrDefault(params, 'mqwVerticalOffsetNm', 0);
    for idx = 1:numel(voltages)
        penetrationNm = estimatePenetrationDepthNm(voltages(idx), opts);
        sigmaNm = max(8, 0.28 * penetrationNm + 0.45 * diffusionNm);
        primary = exp(-0.5 * ((depthGridNm - penetrationNm) / sigmaNm) .^ 2);
        activeGate = 1 ./ (1 + exp((depthGridNm - activeThicknessNm) / max(6, diffusionNm * 0.35)));
        response = primary .* (0.35 + 0.65 * activeGate);
        response = response / max(sum(response), eps);
        responses(:, idx) = response;
        centroidNm(idx) = sum(depthGridNm .* response) / max(sum(response), eps);
        fwhmNm(idx) = estimateFwhm(depthGridNm, response);
    end

    depthResponse = struct();
    depthResponse.depthGridNm = depthGridNm;
    depthResponse.voltagesKeV = voltages;
    depthResponse.response = responses;
    depthResponse.summaryTable = table(voltages, centroidNm, fwhmNm, ...
        'VariableNames', {'Voltage_kV', 'ResponseCentroid_nm', 'DepthResponseFWHM_nm'});
end

function opts = defaultOptions()
    opts = struct();
    opts.voltagesKeV = [];
    opts.wavelengthAxis = [];
    opts.defaultWavelengthRange = [300, 800];
    opts.ganWindow = [350, 405];
    opts.mqwWindow = [430, 580];
    opts.enableOptimization = true;
    opts.maxIterations = 120;
    opts.maxFunctionEvaluations = 520;
    opts.leaveOneOutMaxIterations = 70;
    opts.leaveOneOutMaxFunctionEvaluations = 260;
    opts.toleranceX = 1e-4;
    opts.toleranceFun = 1e-5;
    opts.optimizerDisplay = 'off';
    opts.objectiveNccWeight = 0.12;
    opts.objectivePeakWeight = 0.08;
    opts.objectiveWorstPeakWeight = 0.09;
    opts.objectiveWorstNRMSEWeight = 0.06;
    opts.objectiveRatioWeight = 0.018;
    opts.sensitivityStepFraction = 0.08;
    opts.uncertaintySamples = 48;
    opts.uncertaintyStepFraction = 0.10;
    opts.uncertaintyObjectiveWindow = 0.85;
    opts.pTypeWindow = [374, 398];
    opts.vpitShortwaveWindow = [374, 398];
    opts.depthResponseRangeNm = [0, 450];
    opts.depthResponseStepNm = 2;
    opts.absorptionPathScaleNm = 120;
    opts.roi = [];
    opts.progressEveryScenario = 0;
    opts.storeForwardAblationData = false;
    opts.forwardFitParameterNames = { ...
        'mqwActiveThicknessNm', ...
        'absorptionScale', ...
        'secondaryExcitationScale', ...
        'ganMqwIqeRatio', ...
        'vpitShortwaveChannelScale'};
    opts.forwardFitStepFraction = 0.12;
    opts.forwardFitMaxCandidates = 11;
    opts.storeForwardFitData = false;
    opts.writeOutputs = false;
    opts.outputDirectory = fullfile(pwd, 'outputs', 'structural_inversion');
    opts.parameterNames = { ...
        'vpitDepthNm', ...
        'vpitEffectiveRadiusNm', ...
        'vpitQuenchStrength', ...
        'mqwActiveThicknessNm', ...
        'carrierDiffusionLengthNm', ...
        'absorptionScale', ...
        'secondaryExcitationScale', ...
        'ganMqwIqeRatio', ...
        'vpitShortwaveChannelScale', ...
        'pTypeMqwIqeRatio'};
    opts.parameterBounds = struct( ...
        'vpitDepthNm', [25, 140], ...
        'vpitEffectiveRadiusNm', [25, 180], ...
        'vpitQuenchStrength', [0, 1.8], ...
        'mqwActiveThicknessNm', [70, 180], ...
        'mqwVerticalOffsetNm', [-35, 35], ...
        'carrierDiffusionLengthNm', [10, 140], ...
        'absorptionScale', [0.15, 3.0], ...
        'secondaryExcitationScale', [0, 3.0], ...
        'ganMqwIqeRatio', [0.15, 4.0], ...
        'vpitShortwaveChannelScale', [0, 4.0], ...
        'vpitInCompositionGradientScale', [0, 3.0], ...
        'vpitQcseScreeningScale', [0, 3.0], ...
        'vpitStrainRelaxationScale', [0, 3.0], ...
        'vpitSidewallNonradiativeScale', [0, 3.0], ...
        'detectorShortwaveResponseScale', [0.80, 1.25], ...
        'pTypeMqwIqeRatio', [0.05, 5.0]);
    opts.initialGuess = struct( ...
        'vpitDepthNm', 55, ...
        'vpitEffectiveRadiusNm', 80, ...
        'vpitQuenchStrength', 0.45, ...
        'mqwActiveThicknessNm', 122, ...
        'mqwVerticalOffsetNm', 0, ...
        'carrierDiffusionLengthNm', 60, ...
        'absorptionScale', 1.0, ...
        'secondaryExcitationScale', 1.0, ...
        'ganMqwIqeRatio', 1.0, ...
        'vpitShortwaveChannelScale', 1.40, ...
        'vpitInCompositionGradientScale', 1.20, ...
        'vpitQcseScreeningScale', 1.10, ...
        'vpitStrainRelaxationScale', 0.90, ...
        'vpitSidewallNonradiativeScale', 0.55, ...
        'detectorShortwaveResponseScale', 1.0, ...
        'pTypeMqwIqeRatio', 1.0);
    opts.fixedInputPriors = struct( ...
        'layerSequence', 'Provided by epitaxy/TEM and treated as fixed model geometry.', ...
        'materialFamily', 'GaN/InGaN/sapphire material assignment.', ...
        'refractiveIndexRange', 'Literature/calibration constrained optical constants.', ...
        'nominalMQWArchitecture', 'Used as a prior, not claimed as blind recovery.');
    opts.independentValidationReferences = struct( ...
        'TEM_vpitDepthNm', NaN, ...
        'TEM_mqwThicknessNm', NaN, ...
        'sideViewCL', 'Optional independent spectral-depth reference.', ...
        'sideViewCLProfile', [], ...
        'heldOutVoltage', 'Leave-one-voltage-out prediction metrics.');
end

function [expSeries, baseSeries] = normalizePairedSeries(experimentalData, baselineData, opts)
    expSeries = normalizeSeries(experimentalData, opts);
    baseSeries = normalizeSeries(baselineData, opts);
    n = min(numel(expSeries), numel(baseSeries));
    expSeries = expSeries(1:n);
    baseSeries = baseSeries(1:n);
    if n < 1
        error('precies:structural_inversion:EmptySeries', 'At least one paired data set is required.');
    end
    for idx = 1:n
        if ~isfinite(expSeries(idx).voltage) && isfinite(baseSeries(idx).voltage)
            expSeries(idx).voltage = baseSeries(idx).voltage;
        elseif ~isfinite(baseSeries(idx).voltage) && isfinite(expSeries(idx).voltage)
            baseSeries(idx).voltage = expSeries(idx).voltage;
        end
    end
end

function series = normalizeSeries(data, opts)
    if iscell(data)
        n = numel(data);
        series = repmat(struct('data', [], 'voltage', NaN), n, 1);
        for idx = 1:n
            series(idx).data = data{idx};
            series(idx).voltage = inferVoltage(data{idx}, idx, opts);
        end
        return;
    end
    if isstruct(data) && numel(data) > 1
        n = numel(data);
        series = repmat(struct('data', [], 'voltage', NaN), n, 1);
        for idx = 1:n
            if isfield(data, 'data')
                series(idx).data = data(idx).data;
                series(idx).voltage = inferVoltage(data(idx), idx, opts);
            else
                series(idx).data = data(idx);
                series(idx).voltage = inferVoltage(data(idx), idx, opts);
            end
        end
        return;
    end
    if isstruct(data) && isfield(data, 'data') && isscalar(data)
        series = struct('data', data.data, 'voltage', inferVoltage(data, 1, opts));
    else
        series = struct('data', data, 'voltage', inferVoltage(data, 1, opts));
    end
end

function voltage = inferVoltage(data, idx, opts)
    voltage = NaN;
    if isstruct(data)
        candidates = {'voltage', 'voltageKeV', 'voltage_kV', 'Voltage_kV'};
        for fieldIdx = 1:numel(candidates)
            fieldName = candidates{fieldIdx};
            if isfield(data, fieldName) && ~isempty(data.(fieldName))
                voltage = double(data.(fieldName));
                return;
            end
        end
    end
    if ~isempty(opts.voltagesKeV) && idx <= numel(opts.voltagesKeV)
        voltage = double(opts.voltagesKeV(idx));
    end
end

function bounds = buildBounds(opts)
    n = numel(opts.parameterNames);
    bounds = zeros(n, 2);
    for idx = 1:n
        fieldName = opts.parameterNames{idx};
        if ~isfield(opts.parameterBounds, fieldName)
            error('precies:structural_inversion:MissingBounds', 'Missing bounds for %s.', fieldName);
        end
        bounds(idx, :) = double(opts.parameterBounds.(fieldName));
    end
end

function vector = initialParameterVector(opts, bounds)
    vector = zeros(numel(opts.parameterNames), 1);
    for idx = 1:numel(opts.parameterNames)
        fieldName = opts.parameterNames{idx};
        if isfield(opts.initialGuess, fieldName)
            vector(idx) = double(opts.initialGuess.(fieldName));
        else
            vector(idx) = mean(bounds(idx, :));
        end
        vector(idx) = min(bounds(idx, 2), max(bounds(idx, 1), vector(idx)));
    end
end

function score = objectiveForVector(vector, bounds, expSeries, baseSeries, opts, scenario)
    vector = min(bounds(:, 2), max(bounds(:, 1), vector(:)));
    params = vectorToParameterStruct(vector, opts.parameterNames);
    fittedSeries = synthesizeSeries(baseSeries, params, opts, scenario);
    metrics = evaluateSeriesMetrics(expSeries, fittedSeries, opts);
    score = objectiveFromMetrics(metrics, opts);
end

function score = objectiveFromMetrics(metrics, opts)
    T = metrics.summaryTable;
    nrmse = mean(T.nRMSE, 'omitnan');
    worstNRMSE = max(T.nRMSE, [], 'omitnan');
    nccLoss = mean(1 - T.Average_NCC, 'omitnan');
    dMQW = getMetricTableColumn(T, {'AbsDeltaLambdaMQW_nm', 'Abs_DeltaLambda_MQW_nm'}, zeros(height(T), 1));
    dGaN = getMetricTableColumn(T, {'AbsDeltaLambdaGaN_nm', 'Abs_DeltaLambda_GaN_nm'}, zeros(height(T), 1));
    ratioErr = getMetricTableColumn(T, {'AbsRatioError_percent', 'Abs_Error_IGaN_over_IMQW_percent'}, zeros(height(T), 1));
    peakLoss = mean(abs(dMQW) + abs(dGaN), 'omitnan') / 30;
    worstPeakLoss = max(abs(dMQW), [], 'omitnan') / 20;
    ratioLoss = mean(log1p(abs(ratioErr) / 40), 'omitnan');
    score = nrmse + opts.objectiveNccWeight * nccLoss + ...
        opts.objectivePeakWeight * peakLoss + ...
        getFieldOrDefault(opts, 'objectiveWorstPeakWeight', 0) * worstPeakLoss + ...
        getFieldOrDefault(opts, 'objectiveWorstNRMSEWeight', 0) * worstNRMSE + ...
        opts.objectiveRatioWeight * ratioLoss;
    if ~isfinite(score)
        score = 1e6;
    end
end

function fittedSeries = synthesizeSeries(baseSeries, params, opts, scenario)
    fittedSeries = baseSeries;
    for idx = 1:numel(baseSeries)
        voltage = baseSeries(idx).voltage;
        fittedSeries(idx).data = applyStructuralModel(baseSeries(idx).data, params, voltage, opts, scenario);
    end
end

function correctedData = applyStructuralModel(baseData, params, voltage, opts, scenario)
    [wavelengthNm, cube] = dataToCube(baseData, opts);
    vpitLikelihood = detectVpitLikelihood(wavelengthNm, cube, opts);
    penetrationNm = estimatePenetrationDepthNm(voltage, opts);
    diffusionNm = getFieldOrDefault(params, 'carrierDiffusionLengthNm', 60);
    activeThicknessNm = getFieldOrDefault(params, 'mqwActiveThicknessNm', 150) + ...
        getFieldOrDefault(params, 'mqwVerticalOffsetNm', 0);
    vpitDepthNm = getFieldOrDefault(params, 'vpitDepthNm', 55);
    vpitRadiusNm = getFieldOrDefault(params, 'vpitEffectiveRadiusNm', 80);
    vpitQuench = getFieldOrDefault(params, 'vpitQuenchStrength', 0.45);
    absorptionScale = getFieldOrDefault(params, 'absorptionScale', 1);
    secondaryScale = getFieldOrDefault(params, 'secondaryExcitationScale', 1);
    ganMqwRatio = getFieldOrDefault(params, 'ganMqwIqeRatio', 1);
    vpitShortwaveScale = getFieldOrDefault(params, 'vpitShortwaveChannelScale', 1);
    pTypeMqwRatio = getFieldOrDefault(params, 'pTypeMqwIqeRatio', ganMqwRatio);
    inGradientScale = getFieldOrDefault(params, 'vpitInCompositionGradientScale', ...
        getFieldOrDefault(params, 'vpitInDepletionScale', 1));
    qcseScale = getFieldOrDefault(params, 'vpitQcseScreeningScale', 1);
    strainRelaxScale = getFieldOrDefault(params, 'vpitStrainRelaxationScale', 1);
    sidewallNonradiativeScale = getFieldOrDefault(params, 'vpitSidewallNonradiativeScale', ...
        getFieldOrDefault(params, 'dopingNonradiativeScale', 1));
    detectorShortwaveScale = getFieldOrDefault(params, 'detectorShortwaveResponseScale', 1);
    if ~isfinite(ganMqwRatio)
        ganMqwRatio = 1;
    end
    if ~isfinite(vpitShortwaveScale)
        vpitShortwaveScale = 1;
    end
    if ~isfinite(pTypeMqwRatio)
        pTypeMqwRatio = ganMqwRatio;
    end
    if ~isfinite(inGradientScale)
        inGradientScale = 1;
    end
    if ~isfinite(qcseScale)
        qcseScale = 1;
    end
    if ~isfinite(strainRelaxScale)
        strainRelaxScale = 1;
    end
    if ~isfinite(sidewallNonradiativeScale)
        sidewallNonradiativeScale = 1;
    end
    if ~isfinite(detectorShortwaveScale)
        detectorShortwaveScale = 1;
    end
    if ~isfinite(voltage)
        voltage = 6;
    end

    if getFieldOrDefault(scenario, 'disablePEInteractionVolume', false)
        diffusionNm = max(4, 0.22 * diffusionNm);
    end
    if getFieldOrDefault(scenario, 'disableAbsorption', false)
        absorptionScale = 0;
    end
    if getFieldOrDefault(scenario, 'disableSecondaryExcitation', false)
        secondaryScale = 0;
    end
    if getFieldOrDefault(scenario, 'disableVpitTransport', false)
        vpitLikelihood(:) = 0;
        vpitQuench = 0;
    end

    ganWindow = softBandpass(wavelengthNm, opts.ganWindow, 7);
    pTypeWindow = softBandpass(wavelengthNm, opts.pTypeWindow, 6);
    highBandgapBlueShiftNm = clamp( ...
        5.2 * inGradientScale + 3.2 * qcseScale + 3.8 * strainRelaxScale, 0, 18);
    vpitShortwaveWindow = softBandpass(wavelengthNm, opts.vpitShortwaveWindow - highBandgapBlueShiftNm, 6);
    mqwWindow = softBandpass(wavelengthNm, opts.mqwWindow, 12);
    mqwMainWindow = softBandpass(wavelengthNm, [500, 526], 8);
    mqwRedTailWindow = softBandpass(wavelengthNm, [532, 570], 10);
    bgWindow = max(0, 1 - max(ganWindow, mqwWindow));

    activeOnset = sigmoid((penetrationNm - activeThicknessNm) / max(6, 0.45 * diffusionNm));
    vpitOnset = sigmoid((penetrationNm - vpitDepthNm) / max(5, 0.38 * diffusionNm));
    lowVoltageWeight = sigmoid((6.5 - voltage) / 1.35);
    highVoltageMqwRecovery = sigmoid((voltage - 7.2) / 0.75);
    radiusInfluence = clamp(vpitRadiusNm / 80, 0.45, 2.2);
    vpitTerm = clamp(vpitLikelihood * radiusInfluence, 0, 1);

    highBandgapDrive = clamp( ...
        0.20 + 0.86 * max(vpitShortwaveScale - 1, 0), 0, 2.4) .* ...
        (1 + 0.11 * inGradientScale + 0.07 * qcseScale + 0.08 * strainRelaxScale) .* ...
        max(0.45, 1 - 0.11 * sidewallNonradiativeScale);
    vpitLongwaveCompetitionRaw = ...
        lowVoltageWeight .* vpitTerm .* ...
        (0.12 + 0.22 * max(vpitShortwaveScale - 1, 0) + ...
        0.05 * inGradientScale + 0.04 * qcseScale + 0.04 * strainRelaxScale + ...
        0.08 * vpitQuench - 0.03 * sidewallNonradiativeScale);
    vpitLongwaveCompetitionMap = clamp(vpitLongwaveCompetitionRaw .* (1 - 0.55 * highVoltageMqwRecovery), ...
        0, 0.50);
    vpitMqwRecoveryMap = clamp( ...
        highVoltageMqwRecovery .* vpitTerm .* (0.28 + 0.34 * vpitOnset + 0.10 * activeOnset) .* ...
        max(0.45, 1 - 0.12 * vpitQuench), ...
        0, 0.62);

    ganGainMap = ganMqwRatio * (0.72 + 1.10 * activeOnset + ...
        0.42 * vpitTerm * vpitOnset + ...
        0.46 * lowVoltageWeight .* vpitTerm .* (0.55 + 0.45 * vpitOnset) .* highBandgapDrive);
    vpitQuenchMap = clamp( ...
        0.42 * vpitQuench * vpitTerm .* (1 - 0.35 * vpitOnset) .* (1 - 0.40 * highVoltageMqwRecovery), ...
        0, 0.88);
    mqwGainMap = (1.15 - 0.34 * activeOnset) .* ...
        (1 - vpitQuenchMap) .* ...
        (1 - vpitLongwaveCompetitionMap) .* ...
        (1 + 0.75 * vpitMqwRecoveryMap);
    bgGainMap = 0.90 + 0.10 * sigmoid((penetrationNm - 0.5 * activeThicknessNm) / max(8, diffusionNm));
    pTypeRelativeGain = pTypeMqwRatio / max(ganMqwRatio, 0.05);
    pTypeBoostMap = clamp( ...
        1 + 0.42 * (pTypeRelativeGain - 1) * lowVoltageWeight .* (0.28 + 0.72 * vpitTerm), ...
        0.20, 4.0);
    vpitShortwaveBoostMap = clamp( ...
        1 + lowVoltageWeight .* vpitTerm .* ...
        (0.55 + 0.45 * vpitOnset) .* (1 - 0.16 * vpitQuench) .* highBandgapDrive, ...
        0.18, 5.2);

    if getFieldOrDefault(scenario, 'disableRefractionReflection', false)
        ganGainMap = ones(size(vpitLikelihood)) * (0.90 + 0.20 * activeOnset);
        mqwGainMap = ones(size(vpitLikelihood)) * (0.95 - 0.12 * activeOnset);
        bgGainMap = 1.0;
    end

    pathNm = opts.absorptionPathScaleNm + 0.35 * max(activeThicknessNm, 0);
    shortAbsorptionProfile = ((max(wavelengthNm, 260) / 405) .^ -2.2) .* softBandpass(wavelengthNm, [300, 450], 18);
    absorptionTransfer = exp(-0.0014 * absorptionScale * pathNm * shortAbsorptionProfile);
    interfaceTransfer = 1 + ...
        (0.025 + 0.035 * activeOnset) * mqwMainWindow(:)' - ...
        0.014 * mqwRedTailWindow(:)';
    interfaceTransfer = max(interfaceTransfer, 0.82);
    if getFieldOrDefault(scenario, 'disableRefractionReflection', false)
        interfaceTransfer = ones(size(interfaceTransfer));
    end
    secondaryTransfer = 1 + (0.055 + 0.055 * activeOnset + 0.025 * vpitOnset) * secondaryScale * mqwMainWindow(:)' .* ...
        (0.45 + 0.55 * mean(vpitLikelihood(:), 'omitnan'));
    secondaryShortwaveDepletion = 1 - clamp(0.015 * secondaryScale, 0, 0.06) * softBandpass(wavelengthNm, [350, 430], 18);
    detectorTransfer = 1 + (clamp(detectorShortwaveScale, 0.80, 1.25) - 1) * ...
        softBandpass(wavelengthNm, [350, 430], 28);
    detectorTransfer = detectorTransfer(:)';

    correctedCube = zeros(size(cube));
    spectralTransferBase = bgWindow(:)' * bgGainMap + 0;
    for row = 1:size(cube, 1)
        for col = 1:size(cube, 2)
            ganGain = ganGainMap(row, col);
            mqwGain = mqwGainMap(row, col);
            shortExcess = ganGain * ( ...
                pTypeWindow(:)' * (pTypeBoostMap(row, col) - 1) + ...
                vpitShortwaveWindow(:)' * (vpitShortwaveBoostMap(row, col) - 1));
            transfer = bgWindow(:)' * bgGainMap + ganWindow(:)' * ganGain + mqwWindow(:)' * mqwGain;
            transfer = transfer + shortExcess;
            transfer = max(transfer, 0.02);
            transfer = transfer .* absorptionTransfer(:)' .* interfaceTransfer .* secondaryTransfer .* secondaryShortwaveDepletion(:)' .* detectorTransfer;
            if all(transfer <= eps)
                transfer = spectralTransferBase + 1;
            end
            correctedCube(row, col, :) = reshape(squeeze(cube(row, col, :))' .* transfer, 1, 1, []);
        end
    end
    correctedData = cubeToData(wavelengthNm, correctedCube);
end

function metrics = evaluateSeriesMetrics(expSeries, fittedSeries, opts)
    n = numel(expSeries);
    perVoltageCell = cell(n, 1);
    for idx = 1:n
        metricOpts = struct( ...
            'ganWindow', opts.ganWindow, ...
            'mqwWindow', opts.mqwWindow, ...
            'wavelengthAxis', opts.wavelengthAxis, ...
            'voltage', expSeries(idx).voltage, ...
            'enableOptimization', false);
        perVoltageCell{idx} = precies.pixelwise_inversion('metrics', expSeries(idx).data, fittedSeries(idx).data, metricOpts);
    end
    perVoltage = vertcat(perVoltageCell{:});
    metrics = struct();
    metrics.perVoltage = perVoltage;
    metrics.summaryTable = buildSeriesSummaryTable(perVoltage);
end

function T = buildSeriesSummaryTable(metricsArray)
    n = numel(metricsArray);
    voltage = zeros(n, 1);
    avgNCC = zeros(n, 1);
    stdNCC = zeros(n, 1);
    avgNRMSE = zeros(n, 1);
    stdNRMSE = zeros(n, 1);
    dMQW = zeros(n, 1);
    dGaN = zeros(n, 1);
    dFwhm = zeros(n, 1);
    ratioErr = zeros(n, 1);
    residualRms = zeros(n, 1);
    signedResidualArea = zeros(n, 1);
    for idx = 1:n
        summary = metricsArray(idx).summary;
        voltage(idx) = metricsArray(idx).voltage;
        avgNCC(idx) = summary.averageNCC;
        stdNCC(idx) = summary.stdNCC;
        avgNRMSE(idx) = summary.averageNRMSE;
        stdNRMSE(idx) = summary.stdNRMSE;
        dMQW(idx) = summary.meanAbsDeltaLambdaMQW;
        dGaN(idx) = summary.meanAbsDeltaLambdaGaN;
        dFwhm(idx) = summary.meanAbsDeltaFwhmMQW;
        ratioErr(idx) = summary.meanAbsRatioError;
        residualSummary = summarizeResidualMaps(metricsArray(idx));
        residualRms(idx) = residualSummary.meanResidualRms;
        signedResidualArea(idx) = residualSummary.meanSignedResidualArea;
    end
    T = table(voltage, avgNCC, stdNCC, avgNRMSE, stdNRMSE, dMQW, dGaN, dFwhm, ratioErr, ...
        residualRms, signedResidualArea, ...
        'VariableNames', {'Voltage_kV', 'Average_NCC', 'SD_NCC', 'nRMSE', 'SD_nRMSE', ...
        'AbsDeltaLambdaMQW_nm', 'AbsDeltaLambdaGaN_nm', 'AbsDeltaFWHMMQW_nm', ...
        'AbsRatioError_percent', 'MeanResidualRMS', 'MeanSignedResidualArea'});
end

function [wavelengthNm, cube] = dataToCube(data, opts)
    if isstruct(data) && isfield(data, 'wavelengthIntensityMaps') && ~isempty(data.wavelengthIntensityMaps)
        cube = double(data.wavelengthIntensityMaps);
        if isfield(data, 'wavelengthAxis') && ~isempty(data.wavelengthAxis)
            wavelengthNm = double(data.wavelengthAxis(:));
        else
            wavelengthNm = defaultAxis(size(cube, 3), opts);
        end
        return;
    end

    if isstruct(data) && isfield(data, 'totalSpectra') && ~isempty(data.totalSpectra)
        [wavelengthNm, cube] = spectraCellToCube(data.totalSpectra, opts);
        return;
    end

    if iscell(data)
        [wavelengthNm, cube] = spectraCellToCube(data, opts);
        return;
    end

    if isnumeric(data) && ndims(data) == 3
        cube = double(data);
        wavelengthNm = defaultAxis(size(cube, 3), opts);
        return;
    end

    error('precies:structural_inversion:InvalidData', ...
        'Data must be totalSpectra, wavelengthIntensityMaps, a spectra cell array, or a 3-D numeric cube.');
end

function [wavelengthNm, cube] = spectraCellToCube(spectraCell, opts)
    [rows, cols] = size(spectraCell);
    wavelengthNm = [];
    for idx = 1:numel(spectraCell)
        spectrum = spectraCell{idx};
        if ~isempty(spectrum) && size(spectrum, 2) >= 2
            wavelengthNm = double(spectrum(:, 1));
            break;
        end
    end
    if isempty(wavelengthNm)
        wavelengthNm = defaultAxis(512, opts);
    end
    if ~isempty(opts.wavelengthAxis)
        wavelengthNm = double(opts.wavelengthAxis(:));
    end
    cube = zeros(rows, cols, numel(wavelengthNm));
    for row = 1:rows
        for col = 1:cols
            spectrum = spectraCell{row, col};
            if isempty(spectrum) || size(spectrum, 2) < 2
                continue;
            end
            cube(row, col, :) = interp1(double(spectrum(:, 1)), double(spectrum(:, 2)), ...
                wavelengthNm, 'linear', 0);
        end
    end
end

function wavelengthNm = defaultAxis(n, opts)
    if ~isempty(opts.wavelengthAxis)
        wavelengthNm = double(opts.wavelengthAxis(:));
    else
        wavelengthNm = linspace(opts.defaultWavelengthRange(1), opts.defaultWavelengthRange(2), n)';
    end
end

function data = cubeToData(wavelengthNm, cube)
    [rows, cols, ~] = size(cube);
    totalSpectra = cell(rows, cols);
    for row = 1:rows
        for col = 1:cols
            totalSpectra{row, col} = [wavelengthNm(:), squeeze(cube(row, col, :))];
        end
    end
    data = struct( ...
        'totalSpectra', {totalSpectra}, ...
        'wavelengthAxis', wavelengthNm(:), ...
        'wavelengthIntensityMaps', cube, ...
        'spectralRows', rows, ...
        'spectralCols', cols);
end

function likelihood = detectVpitLikelihood(wavelengthNm, cube, opts)
    shortMask = wavelengthNm >= opts.ganWindow(1) & wavelengthNm <= opts.ganWindow(2);
    longMask = wavelengthNm >= opts.mqwWindow(1) & wavelengthNm <= opts.mqwWindow(2);
    shortMap = sum(cube(:, :, shortMask), 3);
    longMap = sum(cube(:, :, longMask), 3);
    ratioMap = shortMap ./ max(longMap, eps);
    center = median(ratioMap(:), 'omitnan');
    spread = iqr(ratioMap(:));
    if ~isfinite(spread) || spread <= eps
        spread = max(std(ratioMap(:), 0, 'omitnan'), eps);
    end
    likelihood = sigmoid((ratioMap - center) / max(spread, eps));
    likelihood(~isfinite(likelihood)) = 0;
end

function depthNm = estimatePenetrationDepthNm(voltageKeV, opts)
    if ~isfinite(voltageKeV)
        voltageKeV = 6;
    end
    depthNm = 12.5 * voltageKeV ^ 1.55;
    if isfield(opts, 'penetrationDepthTable') && ~isempty(opts.penetrationDepthTable)
        tbl = opts.penetrationDepthTable;
        if istable(tbl) && all(ismember({'Voltage_keV', 'MeanDepth_nm'}, tbl.Properties.VariableNames))
            depthNm = interp1(tbl.Voltage_keV, tbl.MeanDepth_nm, voltageKeV, 'linear', 'extrap');
        elseif isnumeric(tbl) && size(tbl, 2) >= 2
            depthNm = interp1(tbl(:, 1), tbl(:, 2), voltageKeV, 'linear', 'extrap');
        end
    end
    depthNm = max(1, depthNm);
end

function y = softBandpass(x, window, edgeNm)
    x = x(:);
    left = sigmoid((x - window(1)) / edgeNm);
    right = sigmoid((window(2) - x) / edgeNm);
    y = left .* right;
end

function y = sigmoid(x)
    y = 1 ./ (1 + exp(-x));
end

function fwhm = estimateFwhm(x, y)
    y = y(:);
    x = x(:);
    if isempty(y) || max(y) <= 0
        fwhm = NaN;
        return;
    end
    mask = y >= 0.5 * max(y);
    if ~any(mask)
        fwhm = NaN;
    else
        fwhm = max(x(mask)) - min(x(mask));
    end
end

function params = vectorToParameterStruct(vector, names)
    params = struct();
    for idx = 1:numel(names)
        params.(names{idx}) = vector(idx);
    end
end

function vector = parameterStructToVector(params, names)
    vector = zeros(numel(names), 1);
    for idx = 1:numel(names)
        vector(idx) = getFieldOrDefault(params, names{idx}, NaN);
    end
end

function params = resolveBestParameters(input, opts)
    if isstruct(input) && isfield(input, 'bestParameters')
        params = input.bestParameters;
    elseif isstruct(input)
        params = input;
    else
        params = opts.initialGuess;
    end
    params = mergeOptions(opts.initialGuess, params);
end

function u = paramsToUnconstrained(params, bounds)
    params = min(bounds(:, 2) - eps, max(bounds(:, 1) + eps, params(:)));
    scaled = (params - bounds(:, 1)) ./ max(bounds(:, 2) - bounds(:, 1), eps);
    scaled = min(1 - 1e-8, max(1e-8, scaled));
    u = log(scaled ./ (1 - scaled));
end

function params = unconstrainedToParams(u, bounds)
    scaled = 1 ./ (1 + exp(-u(:)));
    params = bounds(:, 1) + scaled .* (bounds(:, 2) - bounds(:, 1));
end

function T = buildParameterTable(params, bounds, opts)
    parameter = string(opts.parameterNames(:));
    value = parameterStructToVector(params, opts.parameterNames);
    lowerBound = bounds(:, 1);
    upperBound = bounds(:, 2);
    role = repmat("inverted structural/output parameter", numel(parameter), 1);
    T = table(parameter, value, lowerBound, upperBound, role, ...
        'VariableNames', {'Parameter', 'Value', 'LowerBound', 'UpperBound', 'Role'});
end

function T = buildIndependentValidationTable(params, refs)
    quantity = strings(2, 1);
    validationSource = strings(2, 1);
    invertedValue = nan(2, 1);
    referenceValue = nan(2, 1);
    unit = strings(2, 1);
    count = 0;

    [count, quantity, validationSource, invertedValue, referenceValue, unit] = appendValidationRow( ...
        count, quantity, validationSource, invertedValue, referenceValue, unit, ...
        "V-pit depth", "vpitDepthNm", "TEM_vpitDepthNm", "nm", params, refs);
    [count, quantity, validationSource, invertedValue, referenceValue, unit] = appendValidationRow( ...
        count, quantity, validationSource, invertedValue, referenceValue, unit, ...
        "MQW active thickness", "mqwActiveThicknessNm", "TEM_mqwThicknessNm", "nm", params, refs);

    if count < 1
        T = table(string.empty(0, 1), string.empty(0, 1), nan(0, 1), nan(0, 1), ...
            nan(0, 1), nan(0, 1), string.empty(0, 1), ...
            'VariableNames', {'Quantity', 'ValidationSource', 'InvertedValue', ...
            'ReferenceValue', 'Difference', 'RelativeDifference_percent', 'Unit'});
        return;
    end
    quantity = quantity(1:count);
    validationSource = validationSource(1:count);
    invertedValue = invertedValue(1:count);
    referenceValue = referenceValue(1:count);
    unit = unit(1:count);
    difference = invertedValue - referenceValue;
    relativeDifference = difference ./ max(abs(referenceValue), eps) * 100;
    T = table(quantity, validationSource, invertedValue, referenceValue, difference, relativeDifference, unit, ...
        'VariableNames', {'Quantity', 'ValidationSource', 'InvertedValue', ...
        'ReferenceValue', 'Difference', 'RelativeDifference_percent', 'Unit'});
end

function [count, quantityList, sourceList, invertedList, referenceList, unitList] = appendValidationRow( ...
        count, quantityList, sourceList, invertedList, referenceList, unitList, ...
        quantity, paramName, refName, unit, params, refs)
    if ~isstruct(params) || ~isfield(params, paramName) || ~isfinite(params.(paramName))
        return;
    end
    if ~isstruct(refs) || ~isfield(refs, refName) || ~isfinite(refs.(refName))
        return;
    end
    count = count + 1;
    quantityList(count) = quantity;
    sourceList(count) = refName;
    invertedList(count) = double(params.(paramName));
    referenceList(count) = double(refs.(refName));
    unitList(count) = unit;
end

function result = buildSideViewCLValidation(depthResponse, refs)
    result = emptySideViewCLValidation();
    [depthNm, intensity] = parseSideViewCLProfile(refs);
    if isempty(depthNm) || isempty(intensity) || ~isstruct(depthResponse) || ...
            ~isfield(depthResponse, 'depthGridNm') || ~isfield(depthResponse, 'response')
        return;
    end

    modelDepth = double(depthResponse.depthGridNm(:));
    modelProfile = mean(double(depthResponse.response), 2, 'omitnan');
    if all(modelProfile <= 0) || numel(modelDepth) < 2
        return;
    end

    [depthNm, order] = sort(depthNm(:));
    intensity = intensity(order);
    valid = isfinite(depthNm) & isfinite(intensity);
    depthNm = depthNm(valid);
    intensity = intensity(valid);
    if numel(depthNm) < 2
        return;
    end

    expProfile = interp1(depthNm, intensity, modelDepth, 'linear', 0);
    expProfile = normalizeProfile(expProfile);
    modelProfile = normalizeProfile(modelProfile);

    residual = expProfile - modelProfile;
    nccValue = calculateProfileNcc(expProfile, modelProfile);
    nrmseValue = sqrt(mean(residual .^ 2, 'omitnan')) / max(max(expProfile) - min(expProfile), eps);
    expCentroid = sum(modelDepth .* expProfile, 'omitnan') / max(sum(expProfile, 'omitnan'), eps);
    modelCentroid = sum(modelDepth .* modelProfile, 'omitnan') / max(sum(modelProfile, 'omitnan'), eps);
    [~, expPeakIdx] = max(expProfile);
    [~, modelPeakIdx] = max(modelProfile);
    expFwhm = estimateFwhm(modelDepth, expProfile);
    modelFwhm = estimateFwhm(modelDepth, modelProfile);

    result = struct();
    result.depthNm = modelDepth;
    result.experimentalProfile = expProfile;
    result.modelProfile = modelProfile;
    result.residual = residual;
    result.summaryTable = table(nccValue, nrmseValue, modelDepth(modelPeakIdx) - modelDepth(expPeakIdx), ...
        modelCentroid - expCentroid, modelFwhm - expFwhm, ...
        'VariableNames', {'NCC', 'nRMSE', 'PeakDepthError_nm', ...
        'CentroidError_nm', 'FWHMError_nm'});
    result.profileTable = table(modelDepth, expProfile, modelProfile, residual, ...
        'VariableNames', {'Depth_nm', 'ExperimentalProfile', 'ModelProfile', 'Residual'});
end

function result = emptySideViewCLValidation()
    result = struct();
    result.depthNm = [];
    result.experimentalProfile = [];
    result.modelProfile = [];
    result.residual = [];
    result.summaryTable = table();
    result.profileTable = table();
end

function [depthNm, intensity] = parseSideViewCLProfile(refs)
    depthNm = [];
    intensity = [];
    if ~isstruct(refs)
        return;
    end
    if isfield(refs, 'sideViewCLDepthNm') && isfield(refs, 'sideViewCLIntensity')
        depthNm = double(refs.sideViewCLDepthNm(:));
        intensity = double(refs.sideViewCLIntensity(:));
        return;
    end
    if isfield(refs, 'sideViewCLProfile') && ~isempty(refs.sideViewCLProfile)
        profile = refs.sideViewCLProfile;
        if isnumeric(profile) && size(profile, 2) >= 2
            depthNm = double(profile(:, 1));
            intensity = double(profile(:, 2));
            return;
        end
        if isstruct(profile) && isfield(profile, 'depthNm') && isfield(profile, 'intensity')
            depthNm = double(profile.depthNm(:));
            intensity = double(profile.intensity(:));
            return;
        end
    end
    if isfield(refs, 'sideViewCL') && isstruct(refs.sideViewCL)
        profile = refs.sideViewCL;
        if isfield(profile, 'profile') && isnumeric(profile.profile) && size(profile.profile, 2) >= 2
            depthNm = double(profile.profile(:, 1));
            intensity = double(profile.profile(:, 2));
        elseif isfield(profile, 'depthNm') && isfield(profile, 'intensity')
            depthNm = double(profile.depthNm(:));
            intensity = double(profile.intensity(:));
        end
    end
end

function y = normalizeProfile(y)
    y = double(y(:));
    y(~isfinite(y)) = 0;
    y = y - min(y);
    if max(y) > 0
        y = y / max(y);
    end
end

function nccValue = calculateProfileNcc(referenceProfile, modelProfile)
    ref = referenceProfile(:) - mean(referenceProfile(:), 'omitnan');
    mdl = modelProfile(:) - mean(modelProfile(:), 'omitnan');
    nccValue = sum(ref .* mdl, 'omitnan') / ...
        (sqrt(sum(ref .^ 2, 'omitnan') * sum(mdl .^ 2, 'omitnan')) + eps);
    if ~isfinite(nccValue)
        nccValue = NaN;
    end
end

function T = buildAblationTable(result)
    n = numel(result.scenarioNames);
    scenario = string(result.scenarioNames(:));
    objective = result.objective(:);
    avgNCC = zeros(n, 1);
    avgNRMSE = zeros(n, 1);
    dMQW = zeros(n, 1);
    dGaN = zeros(n, 1);
    residualRms = zeros(n, 1);
    signedResidualArea = zeros(n, 1);
    for idx = 1:n
        Tm = result.metrics(idx).summaryTable;
        avgNCC(idx) = mean(Tm.Average_NCC, 'omitnan');
        avgNRMSE(idx) = mean(Tm.nRMSE, 'omitnan');
        dMQW(idx) = mean(Tm.AbsDeltaLambdaMQW_nm, 'omitnan');
        dGaN(idx) = mean(Tm.AbsDeltaLambdaGaN_nm, 'omitnan');
        residualSummary = summarizeResidualMaps(result.metrics(idx));
        residualRms(idx) = residualSummary.meanResidualRms;
        signedResidualArea(idx) = residualSummary.meanSignedResidualArea;
    end
    T = table(scenario, objective, avgNCC, avgNRMSE, dMQW, dGaN, residualRms, signedResidualArea, ...
        'VariableNames', {'Scenario', 'Objective', 'Average_NCC', 'nRMSE', ...
        'AbsDeltaLambdaMQW_nm', 'AbsDeltaLambdaGaN_nm', 'MeanResidualRMS', 'MeanSignedResidualArea'});
end

function scenarios = buildForwardAblationScenarios()
    scenarios = repmat(struct('name', '', 'toggles', struct()), 4, 1);
    scenarios(1).name = 'M0';
    scenarios(1).toggles = struct( ...
        'enablePEInteractionVolume', false, ...
        'enablePhotonRefractionReflection', false, ...
        'enableMQWInterfaceOptics', false, ...
        'enableAbsorption', false, ...
        'enableSecondaryExcitation', false);
    scenarios(2).name = 'M1';
    scenarios(2).toggles = struct( ...
        'enablePhotonRefractionReflection', false, ...
        'enableMQWInterfaceOptics', false, ...
        'enableAbsorption', false, ...
        'enableSecondaryExcitation', false);
    scenarios(3).name = 'M2';
    scenarios(3).toggles = struct('enableSecondaryExcitation', false);
    scenarios(4).name = 'M3';
    scenarios(4).toggles = struct();
end

function scenarioOptions = applyForwardAblationScenario(baseOptions, bestParams, scenario)
    if nargin < 1 || ~isstruct(baseOptions)
        baseOptions = struct();
    end
    scenarioOptions = baseOptions;
    currentStructural = getFieldOrDefault(scenarioOptions, 'structuralParameters', struct());
    scenarioOptions.structuralParameters = mergeOptions(currentStructural, bestParams);
    currentToggles = getFieldOrDefault(scenarioOptions, 'modelToggles', struct());
    scenarioOptions.modelToggles = mergeOptions(currentToggles, scenario.toggles);
    if isfield(scenarioOptions, 'params') && isstruct(scenarioOptions.params)
        scenarioOptions.params.structuralParameters = scenarioOptions.structuralParameters;
        scenarioOptions.params.modelToggles = scenarioOptions.modelToggles;
    end
end

function simOptionsCell = normalizeForwardSimulationOptionsSeries(simulationOptionsSeries)
    if isempty(simulationOptionsSeries)
        simOptionsCell = {};
    elseif iscell(simulationOptionsSeries)
        simOptionsCell = simulationOptionsSeries(:);
    elseif isstruct(simulationOptionsSeries) && numel(simulationOptionsSeries) > 1
        simOptionsCell = num2cell(simulationOptionsSeries(:));
    elseif isstruct(simulationOptionsSeries)
        simOptionsCell = {simulationOptionsSeries};
    else
        error('precies:structural_inversion:InvalidSimulationOptions', ...
            'Simulation options must be a struct or cell array of structs.');
    end
end

function [candidateVectors, candidateLabel] = buildForwardFitCandidates(seedVector, bounds, opts)
    selectedNames = opts.forwardFitParameterNames;
    if isempty(selectedNames)
        selectedNames = opts.parameterNames;
    end
    candidates = seedVector(:)';
    labels = "seed";
    for idx = 1:numel(selectedNames)
        paramIdx = find(strcmp(opts.parameterNames, selectedNames{idx}), 1);
        if isempty(paramIdx)
            continue;
        end
        span = bounds(paramIdx, 2) - bounds(paramIdx, 1);
        delta = opts.forwardFitStepFraction * span;
        low = seedVector;
        low(paramIdx) = max(bounds(paramIdx, 1), low(paramIdx) - delta);
        high = seedVector;
        high(paramIdx) = min(bounds(paramIdx, 2), high(paramIdx) + delta);
        candidates = [candidates; low(:)'; high(:)'];
        labels = [labels; string(selectedNames{idx}) + " low"; string(selectedNames{idx}) + " high"];
    end
    [candidateVectors, uniqueIdx] = unique(round(candidates, 10), 'rows', 'stable');
    labels = labels(uniqueIdx);
    maxCandidates = max(1, round(opts.forwardFitMaxCandidates));
    if size(candidateVectors, 1) > maxCandidates
        candidateVectors = candidateVectors(1:maxCandidates, :);
        labels = labels(1:maxCandidates);
    end
    candidateLabel = labels(:);
end

function [metrics, objective, simDataCell] = evaluateForwardFitCandidate(expSeries, simOptionsCell, params, opts)
    nPairs = min(numel(expSeries), numel(simOptionsCell));
    perVoltageCell = cell(nPairs, 1);
    simDataCell = cell(nPairs, 1);
    scenario = struct('name', 'Forward structural fit', 'toggles', struct());
    for idx = 1:nPairs
        scenarioOptions = applyForwardAblationScenario(simOptionsCell{idx}, params, scenario);
        simData = precies.simulation('simulateOptions', scenarioOptions);
        simDataCell{idx} = simData;
        metricExpData = cropDataForMetrics(expSeries(idx).data, opts.roi);
        metricSimData = cropDataForMetrics(simData, opts.roi);
        metricOpts = struct( ...
            'ganWindow', opts.ganWindow, ...
            'mqwWindow', opts.mqwWindow, ...
            'wavelengthAxis', opts.wavelengthAxis, ...
            'voltage', expSeries(idx).voltage, ...
            'enableOptimization', false);
        perVoltageCell{idx} = precies.pixelwise_inversion('metrics', metricExpData, metricSimData, metricOpts);
    end
    perVoltage = vertcat(perVoltageCell{:});
    metrics = struct();
    metrics.perVoltage = perVoltage;
    metrics.summaryTable = buildSeriesSummaryTable(perVoltage);
    objective = objectiveFromMetrics(metrics, opts);
end

function result = emptyForwardFitResult(opts, note)
    result = struct();
    result.method = 'bounded full-forward candidate search';
    result.note = note;
    result.bestParameters = opts.initialGuess;
    result.objective = NaN;
    result.metrics = [];
    result.candidateTable = table();
    result.simulationData = {};
end

function T = buildForwardFitCandidateTable(labels, vectors, objective, avgNCC, avgNRMSE, residualRms, opts)
    candidate = (1:size(vectors, 1))';
    label = string(labels(:));
    T = table(candidate, label, objective(:), avgNCC(:), avgNRMSE(:), residualRms(:), ...
        'VariableNames', {'Candidate', 'Label', 'Objective', 'Average_NCC', 'nRMSE', 'MeanResidualRMS'});
    for idx = 1:numel(opts.parameterNames)
        T.(opts.parameterNames{idx}) = vectors(:, idx);
    end
end

function dataOut = cropDataForMetrics(dataIn, roi)
    dataOut = dataIn;
    if isempty(roi) || numel(roi) ~= 4
        return;
    end
    roi = round(double(roi(:)'));
    rowRange = roi(1):roi(2);
    colRange = roi(3):roi(4);

    if isstruct(dataOut)
        if isfield(dataOut, 'totalSpectra') && ~isempty(dataOut.totalSpectra)
            [rowRangeLocal, colRangeLocal] = clipRowColRange(rowRange, colRange, size(dataOut.totalSpectra, 1), size(dataOut.totalSpectra, 2));
            dataOut.totalSpectra = dataOut.totalSpectra(rowRangeLocal, colRangeLocal);
            dataOut.spectralRows = numel(rowRangeLocal);
            dataOut.spectralCols = numel(colRangeLocal);
        end
        if isfield(dataOut, 'wavelengthIntensityMaps') && ~isempty(dataOut.wavelengthIntensityMaps)
            [rowRangeLocal, colRangeLocal] = clipRowColRange(rowRange, colRange, size(dataOut.wavelengthIntensityMaps, 1), size(dataOut.wavelengthIntensityMaps, 2));
            dataOut.wavelengthIntensityMaps = dataOut.wavelengthIntensityMaps(rowRangeLocal, colRangeLocal, :);
            dataOut.spectralRows = numel(rowRangeLocal);
            dataOut.spectralCols = numel(colRangeLocal);
        end
    elseif iscell(dataOut)
        [rowRangeLocal, colRangeLocal] = clipRowColRange(rowRange, colRange, size(dataOut, 1), size(dataOut, 2));
        dataOut = dataOut(rowRangeLocal, colRangeLocal);
    elseif isnumeric(dataOut) && ndims(dataOut) == 3
        [rowRangeLocal, colRangeLocal] = clipRowColRange(rowRange, colRange, size(dataOut, 1), size(dataOut, 2));
        dataOut = dataOut(rowRangeLocal, colRangeLocal, :);
    end
end

function [rowRangeLocal, colRangeLocal] = clipRowColRange(rowRange, colRange, rows, cols)
    rowRangeLocal = rowRange(rowRange >= 1 & rowRange <= rows);
    colRangeLocal = colRange(colRange >= 1 & colRange <= cols);
    if isempty(rowRangeLocal)
        rowRangeLocal = 1:rows;
    end
    if isempty(colRangeLocal)
        colRangeLocal = 1:cols;
    end
end

function voltage = inferForwardAblationVoltage(options)
    voltage = NaN;
    if isfield(options, 'voltageKeV') && ~isempty(options.voltageKeV)
        voltage = double(options.voltageKeV);
    elseif isfield(options, 'voltage_kV') && ~isempty(options.voltage_kV)
        voltage = double(options.voltage_kV);
    elseif isfield(options, 'params') && isstruct(options.params) && ...
            isfield(options.params, 'electronEnergy') && ~isempty(options.params.electronEnergy)
        voltage = double(options.params.electronEnergy) / 1e3;
    end
end

function T = buildForwardAblationTable(scenarioNames, metrics, objectives)
    n = numel(scenarioNames);
    objective = objectives(:);
    avgNCC = zeros(n, 1);
    avgNRMSE = zeros(n, 1);
    dMQW = zeros(n, 1);
    dGaN = zeros(n, 1);
    dFwhm = zeros(n, 1);
    ratioErr = zeros(n, 1);
    residualRms = zeros(n, 1);
    signedResidualArea = zeros(n, 1);
    for idx = 1:n
        summary = metrics(idx).summary;
        avgNCC(idx) = summary.averageNCC;
        avgNRMSE(idx) = summary.averageNRMSE;
        dMQW(idx) = summary.meanAbsDeltaLambdaMQW;
        dGaN(idx) = summary.meanAbsDeltaLambdaGaN;
        dFwhm(idx) = summary.meanAbsDeltaFwhmMQW;
        ratioErr(idx) = summary.meanAbsRatioError;
        residualSummary = summarizeResidualMaps(metrics(idx));
        residualRms(idx) = residualSummary.meanResidualRms;
        signedResidualArea(idx) = residualSummary.meanSignedResidualArea;
    end
    T = table(scenarioNames(:), objective, avgNCC, avgNRMSE, dMQW, dGaN, dFwhm, ratioErr, ...
        residualRms, signedResidualArea, ...
        'VariableNames', {'Scenario', 'Objective', 'Average_NCC', 'nRMSE', ...
        'AbsDeltaLambdaMQW_nm', 'AbsDeltaLambdaGaN_nm', 'AbsDeltaFWHMMQW_nm', ...
        'AbsRatioError_percent', 'MeanResidualRMS', 'MeanSignedResidualArea'});
end

function residualSummary = summarizeResidualMaps(metric)
    residualSummary = struct('meanResidualRms', NaN, 'meanSignedResidualArea', NaN);
    if isstruct(metric) && isfield(metric, 'residualRmsMap') && ~isempty(metric.residualRmsMap)
        residualSummary.meanResidualRms = mean(metric.residualRmsMap(:), 'omitnan');
    end
    if isstruct(metric) && isfield(metric, 'signedResidualAreaMap') && ~isempty(metric.signedResidualAreaMap)
        residualSummary.meanSignedResidualArea = mean(metric.signedResidualAreaMap(:), 'omitnan');
    end
end

function values = getMetricTableColumn(T, candidateNames, defaultValues)
    values = defaultValues;
    if ~istable(T)
        return;
    end
    for idx = 1:numel(candidateNames)
        name = candidateNames{idx};
        if ismember(name, T.Properties.VariableNames)
            values = T.(name);
            return;
        end
    end
end

function result = emptySensitivityResult(opts)
    parameter = string(opts.parameterNames(:));
    nanVec = nan(numel(parameter), 1);
    result = struct('baseObjective', NaN, 'table', table(parameter, nanVec, nanVec, nanVec, nanVec, ...
        'VariableNames', {'Parameter', 'BestValue', 'LowerPerturbObjective', ...
        'UpperPerturbObjective', 'SensitivityScore'}));
end

function result = emptyUncertaintyResult(opts)
    parameter = string(opts.parameterNames(:));
    nanVec = nan(numel(parameter), 1);
    result = struct('samples', [], 'objective', [], 'acceptedMask', [], ...
        'parameterCI', table(parameter, nanVec, nanVec, nanVec, ...
        'VariableNames', {'Parameter', 'BestValue', 'CI16', 'CI84'}));
end

function options = mergeOptions(defaults, override)
    options = defaults;
    if nargin < 2 || ~isstruct(override)
        return;
    end
    fields = fieldnames(override);
    for idx = 1:numel(fields)
        name = fields{idx};
        if isstruct(override.(name)) && isfield(options, name) && isstruct(options.(name))
            options.(name) = mergeOptions(options.(name), override.(name));
        else
            options.(name) = override.(name);
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

function value = clamp(value, lowerBound, upperBound)
    value = max(lowerBound, min(upperBound, value));
end

function writeOutputs(result, opts)
    if ~exist(opts.outputDirectory, 'dir')
        mkdir(opts.outputDirectory);
    end
    writetable(result.parameterTable, fullfile(opts.outputDirectory, 'structural_parameters.csv'));
    if isfield(result, 'independentValidation') && ~isempty(result.independentValidation)
        writetable(result.independentValidation, fullfile(opts.outputDirectory, 'independent_validation.csv'));
    end
    if isfield(result, 'sideViewCLValidation') && isstruct(result.sideViewCLValidation) && ...
            isfield(result.sideViewCLValidation, 'summaryTable') && height(result.sideViewCLValidation.summaryTable) > 0
        writetable(result.sideViewCLValidation.summaryTable, fullfile(opts.outputDirectory, 'sideview_cl_validation.csv'));
        writetable(result.sideViewCLValidation.profileTable, fullfile(opts.outputDirectory, 'sideview_cl_profile_fit.csv'));
    end
    writetable(result.metrics.summaryTable, fullfile(opts.outputDirectory, 'structural_fit_metrics.csv'));
    writetable(result.sensitivity.table, fullfile(opts.outputDirectory, 'structural_sensitivity.csv'));
    writetable(result.uncertainty.parameterCI, fullfile(opts.outputDirectory, 'structural_uncertainty_ci.csv'));
    writetable(result.depthResponse.summaryTable, fullfile(opts.outputDirectory, 'depth_response.csv'));
    save(fullfile(opts.outputDirectory, 'structural_inversion_result.mat'), 'result');
end

function result = selfTest()
    opts = defaultOptions();
    opts.voltagesKeV = [5, 7, 10];
    opts.wavelengthAxis = linspace(330, 650, 240)';
    opts.maxIterations = 45;
    opts.maxFunctionEvaluations = 170;
    opts.uncertaintySamples = 8;
    opts.sensitivityStepFraction = 0.05;
    opts.independentValidationReferences.sideViewCLProfile = [ ...
        opts.depthResponseRangeNm(1):25:opts.depthResponseRangeNm(2); ...
        exp(-0.5 * (((opts.depthResponseRangeNm(1):25:opts.depthResponseRangeNm(2)) - 120) / 55) .^ 2)]';
    [experimentalSeries, baselineSeries] = buildSyntheticSeries(opts);
    result = runStructuralInversion(experimentalSeries, baselineSeries, opts);
    initialNRMSE = mean(result.initialMetrics.summaryTable.nRMSE, 'omitnan');
    fittedNRMSE = mean(result.metrics.summaryTable.nRMSE, 'omitnan');
    assert(fittedNRMSE <= initialNRMSE, ...
        'Structural inversion did not improve the global nRMSE in self-test.');
    ablationResult = runAblation(experimentalSeries, baselineSeries, result, opts);
    assert(height(ablationResult.summaryTable) >= 4, 'Ablation summary was not generated.');
    assert(isfield(result, 'sideViewCLValidation') && height(result.sideViewCLValidation.summaryTable) == 1, ...
        'Side-view CL validation summary was not generated.');
    looOpts = opts;
    looOpts.maxIterations = 12;
    looOpts.maxFunctionEvaluations = 60;
    looOpts.uncertaintySamples = 0;
    looResult = runLeaveOneVoltageOut(experimentalSeries, baselineSeries, looOpts);
    assert(height(looResult.summaryTable) == numel(opts.voltagesKeV), ...
        'Leave-one-voltage-out summary has the wrong size.');
    fprintf('precies.structural_inversion selftest: nRMSE %.4f -> %.4f\n', initialNRMSE, fittedNRMSE);
end

function [experimentalSeries, baselineSeries] = buildSyntheticSeries(opts)
    rows = 3;
    cols = 3;
    wl = opts.wavelengthAxis(:);
    [xGrid, yGrid] = meshgrid(linspace(-1, 1, cols), linspace(-1, 1, rows));
    vpitMap = exp(-2.5 * (xGrid .^ 2 + yGrid .^ 2));
    baselineSeries = cell(numel(opts.voltagesKeV), 1);
    experimentalSeries = cell(numel(opts.voltagesKeV), 1);
    trueParams = opts.initialGuess;
    trueParams.vpitDepthNm = 68;
    trueParams.vpitEffectiveRadiusNm = 95;
    trueParams.vpitQuenchStrength = 0.72;
    trueParams.mqwActiveThicknessNm = 158;
    trueParams.carrierDiffusionLengthNm = 48;
    trueParams.absorptionScale = 1.35;
    trueParams.secondaryExcitationScale = 1.45;
    trueParams.ganMqwIqeRatio = 1.28;
    trueParams.vpitShortwaveChannelScale = 1.65;
    trueParams.pTypeMqwIqeRatio = 1.85;

    for voltageIdx = 1:numel(opts.voltagesKeV)
        voltage = opts.voltagesKeV(voltageIdx);
        ganAmp = 0.10 + 0.08 * voltage;
        mqwAmp = 1.0 - 0.025 * voltage;
        cube = zeros(rows, cols, numel(wl));
        for row = 1:rows
            for col = 1:cols
                gan = ganAmp * (1 + 0.65 * vpitMap(row, col)) * exp(-0.5 * ((wl - 390) / 12) .^ 2);
                mqw = mqwAmp * (1 - 0.18 * vpitMap(row, col)) * exp(-0.5 * ((wl - 520) / 22) .^ 2);
                cube(row, col, :) = gan + mqw + 0.02;
            end
        end
        baseData = cubeToData(wl, cube);
        baseData.voltage = voltage;
        expData = applyStructuralModel(baseData, trueParams, voltage, opts, struct());
        expData.voltage = voltage;
        baselineSeries{voltageIdx} = baseData;
        experimentalSeries{voltageIdx} = expData;
    end
end
