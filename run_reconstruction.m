function result = run_reconstruction(options)
if nargin < 1
    options = struct();
end

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
oldDir = pwd;
cd(rootDir);
cleanupObj = onCleanup(@() cd(oldDir));
opts = precies.workflow('defaultOptions', rootDir);
opts = precies.workflow('mergeStructs', opts, options);
if ~exist(opts.outputDir, 'dir')
    mkdir(opts.outputDir);
end

fprintf('Loading example TIFF...\n');
example = precies.workflow('loadExampleExperiment', rootDir, opts);
params = opts.params;
params.scanX = example.importedTiffData.scanX_nm;
params.scanY = example.importedTiffData.scanY_nm;
params.totalPixels = example.importedTiffData.pixelsX * example.importedTiffData.pixelsY;
params.electronEnergy = opts.voltageKeV * 1e3;
params.mqw_barriers = max(0, round(params.mqw_pairs) - 1);
fprintf('Example data: %d x %d pixels, %.1f x %.1f nm, %d V-pits.\n', ...
    example.importedTiffData.pixelsY, example.importedTiffData.pixelsX, ...
    example.importedTiffData.scanY_nm, example.importedTiffData.scanX_nm, ...
    numel(example.importedTiffData.vPitRadii_nm));

[calibrationProfile, localField] = precies.workflow('calibrateExample', ...
    example.experimentalData, example.importedTiffData, params);
params.calibrationProfile = calibrationProfile;
if ~isempty(localField)
    params.localParameterField = localField;
end

simulationOptions = precies.workflow('buildSimulationOptions', example.importedTiffData, params, opts);
fprintf('Running 5 kV simulation over %d pixels...\n', params.totalPixels);
simulationData = precies.simulation('simulateOptions', simulationOptions);
fprintf('Simulation complete.\n');

mappingFiles = precies.workflow('writeSpectralTiff', simulationData, opts.outputDir, '5KV_simulation_mapping');
fprintf('Simulation mapping saved to %s\n', mappingFiles.stackTif);

metricOptions = struct( ...
    'voltage', opts.voltageKeV, ...
    'roi', [1, example.importedTiffData.pixelsY, 1, example.importedTiffData.pixelsX], ...
    'wavelengthAxis', example.experimentalData.wavelengthAxis, ...
    'enableOptimization', false, ...
    'storeResidualCube', false);
metrics = precies.pixelwise_inversion('metrics', example.experimentalData, simulationData, metricOptions);
metricsPath = fullfile(opts.outputDir, '5KV_metrics.csv');
writetable(metrics.summaryTable, metricsPath);
fprintf('Metrics saved to %s\n', metricsPath);
fprintf('Average NCC: %.4f\n', metrics.summary.averageNCC);
fprintf('nRMSE: %.4f\n', metrics.summary.averageNRMSE);
table1 = precies.workflow('buildTable1Validation', metrics.summaryTable);
table1Path = fullfile(opts.outputDir, 'Table1_validation_metrics.csv');
writetable(table1, table1Path);
fprintf('Table 1 validation metrics saved to %s\n', table1Path);

activityPackagePath = fullfile(opts.outputDir, 'activity_structure_package.mat');
activityPackage = precies.workflow('exportActivityStructure', ...
    simulationData, example.importedTiffData, params, activityPackagePath);
fprintf('Activity-Structure package saved to %s\n', activityPackagePath);

result = struct();
result.rootDir = rootDir;
result.outputDir = opts.outputDir;
result.exampleTiff = example.tiffPath;
result.params = params;
result.simulationOptions = simulationOptions;
result.simulationData = simulationData;
result.metrics = metrics;
result.mappingFiles = mappingFiles;
result.activityPackagePath = activityPackagePath;
result.activityPackage = activityPackage;

runAblation = logical(opts.runAblation);
if opts.promptForAblation && usejava('desktop')
    choice = questdlg('Run the forward ablation experiment now?', ...
        'Forward Ablation', 'Yes', 'No', 'No');
    runAblation = strcmp(choice, 'Yes');
end

if runAblation
    fprintf('Running forward ablation...\n');
    ablationOptions = struct( ...
        'roi', metricOptions.roi, ...
        'wavelengthAxis', metricOptions.wavelengthAxis, ...
        'progressEveryScenario', 1, ...
        'storeForwardAblationData', false);
    ablationResult = precies.structural_inversion('forwardAblation', ...
        example.experimentalData, simulationOptions, struct(), ablationOptions);
    table2 = precies.workflow('buildTable2Ablation', ablationResult.summaryTable);
    table2Files = precies.workflow('writeTable2Files', table2, opts.outputDir);
    result.ablationResult = ablationResult;
    result.table2 = table2;
    result.table2Files = table2Files;
    fprintf('Forward ablation table saved to %s\n', table2Files.csv);
else
    fprintf('Forward ablation skipped.\n');
end

save(fullfile(opts.outputDir, '5KV_reconstruction_result.mat'), 'result', '-v7.3');
fprintf('Run package saved to %s\n', fullfile(opts.outputDir, '5KV_reconstruction_result.mat'));
end
