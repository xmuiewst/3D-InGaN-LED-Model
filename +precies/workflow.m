function varargout = workflow(action, varargin)
switch lower(string(action))
    case "defaultoptions"
        varargout{1} = defaultOptions(varargin{:});
    case "mergestructs"
        varargout{1} = mergeStructs(varargin{:});
    case "loadexampleexperiment"
        varargout{1} = loadExampleExperiment(varargin{:});
    case "calibrateexample"
        [varargout{1:nargout}] = calibrateExample(varargin{:});
    case "buildsimulationoptions"
        varargout{1} = buildSimulationOptions(varargin{:});
    case "writespectraltiff"
        varargout{1} = writeSpectralTiff(varargin{:});
    case "buildtable1validation"
        varargout{1} = buildTable1Validation(varargin{:});
    case "exportactivitystructure"
        varargout{1} = exportActivityStructure(varargin{:});
    case "buildtable2ablation"
        varargout{1} = buildTable2Ablation(varargin{:});
    case "writetable2files"
        varargout{1} = writeTable2Files(varargin{:});
    case "loadactivitypackage"
        varargout{1} = loadActivityPackage(varargin{:});
    case "buildlegacyvolume"
        varargout{1} = buildLegacyVolume(varargin{:});
    case "buildactivitystructurevolume"
        varargout{1} = buildActivityStructureVolume(varargin{:});
    case "showlegacysliceviewer"
        showLegacySliceViewer(varargin{:});
    case "showactivitystructurevolume"
        showActivityStructureVolume(varargin{:});
    case "resolvecolormap"
        varargout{1} = resolveColormap(varargin{:});
    otherwise
        error('precies:workflow:InvalidAction', 'Unsupported workflow action: %s', action);
end
end

function opts = defaultOptions(rootDir)
if nargin < 1 || isempty(rootDir)
    rootDir = pwd;
end
opts = struct();
opts.tiffFile = fullfile(rootDir, 'data', '5KV-Example.tif');
opts.outputDir = fullfile(rootDir, 'output');
opts.voltageKeV = 5;
opts.rngSeed = 20260512;
opts.progressEveryPixels = 1;
opts.promptForAblation = true;
opts.runAblation = false;
opts.params = struct( ...
    'InteTime', 0.3, ...
    'numRays', 1000, ...
    'electronEnergy', 5e3, ...
    'beamCurrent', 100e-12, ...
    'depth_step', 100, ...
    'mqw_pairs', 10, ...
    'mqw_barriers', 9, ...
    'scanX', 500, ...
    'scanY', 500, ...
    'totalPixels', 100, ...
    'gridPoints3D', 100, ...
    'barrierThick', 8, ...
    'wellThick', 5, ...
    'elementaryCharge', 1.602e-19, ...
    'numElectrons', 100);
end

function out = mergeStructs(base, override)
out = base;
if nargin < 2 || isempty(override) || ~isstruct(override)
    return;
end
names = fieldnames(override);
for idx = 1:numel(names)
    name = names{idx};
    value = override.(name);
    if isstruct(value) && isfield(out, name) && isstruct(out.(name)) && isscalar(value) && isscalar(out.(name))
        out.(name) = mergeStructs(out.(name), value);
    else
        out.(name) = value;
    end
end
end

function example = loadExampleExperiment(rootDir, opts)
tiffPath = findExampleTiff(rootDir, opts.tiffFile);
mapData = loadTIFFForVpitsFile(tiffPath);
spectraData = mapData.spectraData;
rows = mapData.spectralRows;
cols = mapData.spectralCols;
centers = mapData.vPitCentroids_nm;
radii = mapData.vPitRadii_nm;
orientations = mapData.vPitOrientations_deg;
if isempty(centers)
    centers = zeros(0, 2);
    radii = zeros(0, 1);
    orientations = zeros(0, 1);
end
validMask = all(isfinite(centers), 2) & isfinite(radii(:)) & radii(:) > 0;
centers = centers(validMask, :);
radii = radii(validMask);
if numel(orientations) >= numel(validMask)
    orientations = orientations(validMask);
else
    orientations = zeros(numel(radii), 1);
end
importedTiffData = struct( ...
    'scanX_nm', mapData.scanX_nm, ...
    'scanY_nm', mapData.scanY_nm, ...
    'pixelsX', cols, ...
    'pixelsY', rows, ...
    'vPitCentroids_nm', centers, ...
    'vPitRadii_nm', radii(:), ...
    'vPitOrientations_deg', orientations(:), ...
    'totalSpectra', {spectraData}, ...
    'spectralRows', rows, ...
    'spectralCols', cols, ...
    'wavelengthAxis', mapData.wavelengthAxis, ...
    'sourceFile', mapData.sourceFile, ...
    'sourceName', mapData.sourceName);
experimentalData = struct( ...
    'totalSpectra', {spectraData}, ...
    'depthSpectra', {{}}, ...
    'layerSpectra', {{}}, ...
    'isImportedTiff', true, ...
    'spectralRows', rows, ...
    'spectralCols', cols, ...
    'wavelengthAxis', mapData.wavelengthAxis, ...
    'sourceFile', mapData.sourceFile, ...
    'sourceName', mapData.sourceName);
example = struct( ...
    'tiffPath', tiffPath, ...
    'mapData', mapData, ...
    'importedTiffData', importedTiffData, ...
    'experimentalData', experimentalData);
end

function tiffPath = findExampleTiff(rootDir, requestedPath)
if nargin >= 2 && ~isempty(requestedPath) && isfile(requestedPath)
    tiffPath = requestedPath;
    return;
end
candidate = fullfile(rootDir, 'data', '5KV-Example.tif');
if isfile(candidate)
    tiffPath = candidate;
    return;
end
matches = dir(fullfile(rootDir, '**', '5KV-Example.tif'));
if isempty(matches)
    error('precies:workflow:MissingExampleTiff', '5KV-Example.tif was not found under %s.', rootDir);
end
tiffPath = fullfile(matches(1).folder, matches(1).name);
end

function [profile, localField] = calibrateExample(experimentalData, importedTiffData, params)
roi = [1, importedTiffData.pixelsY, 1, importedTiffData.pixelsX];
commonAxisNm = getCommonAxisNm(experimentalData);
[meanSpectrum, validPixels] = buildMeanNormalizedSpectrum(experimentalData.totalSpectra, roi, commonAxisNm);
fitResult = fitAdaptiveGaussianMixture(commonAxisNm, meanSpectrum);
profile = buildWorkflowCalibrationProfile(experimentalData, roi, fitResult, params);
profile.validPixelCount = validPixels;
profile.wavelengthAxis = commonAxisNm(:);
profile.meanSpectrum = meanSpectrum(:);
profile.fittedSpectrum = fitResult.fittedSpectrum(:);
localField = buildExperimentalLocalParameterField(experimentalData, commonAxisNm, profile);
end

function profile = buildWorkflowCalibrationProfile(experimentalData, roi, fitResult, params)
baseConfig = precies.config('wavelength', params);
layerParams = baseConfig.layerParams;
defaultNGaNSpectrum = getLayerSpectrum(layerParams, 'n-GaN');
defaultPrestrainedSpectrum = getLayerSpectrum(layerParams, 'Prestrained');
defaultBarrierSpectrum = getLayerSpectrum(layerParams, 'MQW-Barrier');
defaultWellSpectrum = getLayerSpectrum(layerParams, 'MQW-Well');
defaultPTypeSpectrum = getLayerSpectrum(layerParams, 'p-EBL');
defaultPrestrainedIn = getLayerComposition(layerParams, 'Prestrained');
defaultWellIn = getLayerComposition(layerParams, 'MQW-Well');
defaultBarrierIn = getLayerComposition(layerParams, 'MQW-Barrier');
defaultEblIn = getLayerComposition(layerParams, 'p-EBL');
componentsNm = fitResult.componentsNm;
familyDefaults = struct( ...
    'nGaN', convertSpectrumModelToComponentsNm(defaultNGaNSpectrum), ...
    'prestrained', convertSpectrumModelToComponentsNm(defaultPrestrainedSpectrum), ...
    'barrier', convertSpectrumModelToComponentsNm(defaultBarrierSpectrum), ...
    'well', convertSpectrumModelToComponentsNm(defaultWellSpectrum), ...
    'pType', convertSpectrumModelToComponentsNm(defaultPTypeSpectrum));
mainPeakNm = getFieldOrDefault(fitResult.windowMetrics, 'mainPeakNm', fitResult.redComponent(2));
assignments = assignComponentsToEmissionFamilies(componentsNm, familyDefaults, mainPeakNm);
nGaNComponentsNm = assignments.nGaNComponentsNm;
prestrainedComponentsNm = assignments.prestrainedComponentsNm;
barrierComponentsNm = assignments.barrierComponentsNm;
wellComponentsNm = assignments.wellComponentsNm;
semipolarComponentsNm = assignments.semipolarComponentsNm;
pTypeComponentsNm = assignments.pTypeComponentsNm;
nGaNComponent = dominantComponentFromGroup(nGaNComponentsNm);
prestrainedComponent = dominantComponentFromGroup(prestrainedComponentsNm);
barrierComponent = dominantComponentFromGroup(barrierComponentsNm);
wellComponent = dominantComponentFromGroup(wellComponentsNm);
semipolarComponent = dominantComponentFromGroup(semipolarComponentsNm);
pTypeComponent = dominantComponentFromGroup(pTypeComponentsNm);
bluePeakNm = barrierComponent(2);
redPeakNm = wellComponent(2);
blueFwhmNm = barrierComponent(3);
redFwhmNm = wellComponent(3);
blueWeight = barrierComponent(1);
redWeight = wellComponent(1);
wellSpectrum = buildGaussianMixtureSpectrum(wellComponentsNm);
barrierSpectrum = buildGaussianMixtureSpectrum(barrierComponentsNm);
nGaNSpectrum = buildGaussianMixtureSpectrum(nGaNComponentsNm);
prestrainedSpectrum = buildGaussianMixtureSpectrum(prestrainedComponentsNm);
pTypeSpectrum = buildGaussianMixtureSpectrum(pTypeComponentsNm);
wellInComposition = invertInComposition(redPeakNm, defaultWellIn);
barrierInComposition = invertInComposition(bluePeakNm, defaultBarrierIn);
barrierInComposition = min(barrierInComposition, wellInComposition - 0.02);
barrierInComposition = clampValue(barrierInComposition, 0.01, max(0.02, wellInComposition - 0.02));
prestrainedInComposition = clampValue(invertInComposition(prestrainedComponent(2), defaultPrestrainedIn), 0.02, 0.18);
eblInComposition = clampValue(invertInComposition(pTypeComponent(2), defaultEblIn), 0.05, 0.30);
semipolarSpectrum = buildGaussianMixtureSpectrum(semipolarComponentsNm);
profile = struct( ...
    'profileVersion', 8, ...
    'createdAt', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'entryMode', 'wavelength', ...
    'sourceFile', getFieldOrDefault(experimentalData, 'sourceFile', ''), ...
    'sourceName', getFieldOrDefault(experimentalData, 'sourceName', 'Imported TIFF'), ...
    'roi', roi, ...
    'fitError', fitResult.fitError, ...
    'qualityThreshold', 0.12, ...
    'wellInComposition', clampValue(wellInComposition, 0.05, 0.95), ...
    'barrierInComposition', barrierInComposition, ...
    'prestrainedInComposition', prestrainedInComposition, ...
    'eblInComposition', eblInComposition, ...
    'wellThick', params.wellThick, ...
    'barrierThick', params.barrierThick, ...
    'nGaNSpectrum', nGaNSpectrum, ...
    'gaNSpectrum', nGaNSpectrum, ...
    'prestrainedSpectrum', prestrainedSpectrum, ...
    'wellSpectrum', wellSpectrum, ...
    'barrierSpectrum', barrierSpectrum, ...
    'pTypeSpectrum', pTypeSpectrum, ...
    'semipolarSpectrum', semipolarSpectrum, ...
    'matchMetric', 'normalized_rmse_max_normalized', ...
    'fitMethod', 'roi_mean_adaptive_multi_pseudo_voigt', ...
    'numPeaks', fitResult.numPeaks, ...
    'componentsNm', componentsNm, ...
    'nGaNComponentsNm', nGaNComponentsNm, ...
    'gaNComponentsNm', nGaNComponentsNm, ...
    'prestrainedComponentsNm', prestrainedComponentsNm, ...
    'barrierComponentsNm', barrierComponentsNm, ...
    'wellComponentsNm', wellComponentsNm, ...
    'pTypeComponentsNm', pTypeComponentsNm, ...
    'semipolarComponentsNm', semipolarComponentsNm, ...
    'dominantPeakNm', fitResult.dominantComponent(2), ...
    'nGaNPeakNm', nGaNComponent(2), ...
    'gaNPeakNm', nGaNComponent(2), ...
    'prestrainedPeakNm', prestrainedComponent(2), ...
    'barrierPeakNm', barrierComponent(2), ...
    'wellPeakNm', wellComponent(2), ...
    'semipolarPeakNm', semipolarComponent(2), ...
    'pTypePeakNm', pTypeComponent(2), ...
    'bluePeakNm', bluePeakNm, ...
    'redPeakNm', redPeakNm, ...
    'blueFwhmNm', blueFwhmNm, ...
    'redFwhmNm', redFwhmNm, ...
    'blueWeight', blueWeight, ...
    'redWeight', redWeight, ...
    'windowMetrics', fitResult.windowMetrics, ...
    'shortwavePeakNm', fitResult.shortwavePeakNm, ...
    'shoulderPeakNm', fitResult.shoulderPeakNm, ...
    'mainPeakNm', mainPeakNm, ...
    'shortwaveFitError', fitResult.shortwaveFitError, ...
    'shoulderFitError', fitResult.shoulderFitError, ...
    'shortwaveToMainAreaRatio', fitResult.shortwaveToMainAreaRatio, ...
    'shoulderToMainAreaRatio', fitResult.shoulderToMainAreaRatio, ...
    'mainPeakArea', fitResult.mainPeakArea);
end

function simulationOptions = buildSimulationOptions(importedTiffData, params, opts)
centers = importedTiffData.vPitCentroids_nm;
radii = importedTiffData.vPitRadii_nm;
orientations = importedTiffData.vPitOrientations_deg;
rawVpits = cell(numel(radii), 1);
for idx = 1:numel(radii)
    orientationDeg = 0;
    if numel(orientations) >= idx && isfinite(orientations(idx))
        orientationDeg = orientations(idx);
    end
    rawVpits{idx} = struct( ...
        'center', centers(idx, :), ...
        'topRadius', radii(idx), ...
        'orientationDeg', orientationDeg);
end
simulationOptions = struct();
simulationOptions.entryMode = 'wavelength';
simulationOptions.voltageKeV = opts.voltageKeV;
simulationOptions.params = params;
simulationOptions.Vpits = rawVpits;
simulationOptions.defaultVpitCount = 0;
simulationOptions.rngSeed = opts.rngSeed;
simulationOptions.progressEveryPixels = opts.progressEveryPixels;
end

function files = writeSpectralTiff(simulationData, outputDir, stem)
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
[wavelengthAxis, cube] = spectraToCube(simulationData.totalSpectra);
stackPath = fullfile(outputDir, [stem '.tif']);
integratedPath = fullfile(outputDir, [stem '_integrated.tif']);
axisPath = fullfile(outputDir, [stem '_wavelength_axis.csv']);
scalePath = fullfile(outputDir, [stem '_scale.csv']);
if isfile(stackPath)
    delete(stackPath);
end
if isfile(integratedPath)
    delete(integratedPath);
end
integrated = sum(cube, 3, 'omitnan');
spectraData = cubeToSpectraCell(wavelengthAxis, cube);
params = getFieldOrDefault(simulationData, 'simulationParams', struct());
scanX = getFieldOrDefault(params, 'scanX', size(cube, 2));
scanY = getFieldOrDefault(params, 'scanY', size(cube, 1));
writeOdemisCompatibleTiff(stackPath, spectraData, integrated, wavelengthAxis, scanX, scanY, 'Total mapping');
integratedAxis = mean(wavelengthAxis(isfinite(wavelengthAxis)));
if ~isfinite(integratedAxis)
    integratedAxis = 0;
end
integratedSpectra = integratedMapToSpectraCell(integrated, integratedAxis);
writeOdemisCompatibleTiff(integratedPath, integratedSpectra, integrated, integratedAxis, scanX, scanY, 'Integrated mapping');
scaleValue = max(cube(:), [], 'omitnan');
if ~isfinite(scaleValue) || scaleValue <= 0
    scaleValue = 1;
end
integratedScale = max(integrated(:), [], 'omitnan');
if ~isfinite(integratedScale) || integratedScale <= 0
    integratedScale = 1;
end
writetable(table(wavelengthAxis(:), 'VariableNames', {'Wavelength_nm'}), axisPath);
writetable(table(scaleValue, integratedScale, 'VariableNames', {'SpectralStackMaxRaw', 'IntegratedMapMaxRaw'}), scalePath);
files = struct('stackTif', stackPath, 'integratedTif', integratedPath, 'wavelengthAxisCsv', axisPath, 'scaleCsv', scalePath);
end

function spectraData = cubeToSpectraCell(wavelengthAxis, cube)
[rows, cols, ~] = size(cube);
spectraData = cell(rows, cols);
for rowIdx = 1:rows
    for colIdx = 1:cols
        spectraData{rowIdx, colIdx} = [wavelengthAxis(:), reshape(double(cube(rowIdx, colIdx, :)), [], 1)];
    end
end
end

function spectraData = integratedMapToSpectraCell(integratedMap, wavelengthValue)
[rows, cols] = size(integratedMap);
spectraData = cell(rows, cols);
for rowIdx = 1:rows
    for colIdx = 1:cols
        spectraData{rowIdx, colIdx} = [wavelengthValue, double(integratedMap(rowIdx, colIdx))];
    end
end
end

function writeOdemisCompatibleTiff(fullpath, spectraData, mappingData, wavelengthAxis, scanX_nm, scanY_nm, exportLabel)
warnState = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
cleanupObj = onCleanup(@() warning(warnState));
spectralRows = size(spectraData, 1);
spectralCols = size(spectraData, 2);
surveyRows = spectralRows;
surveyCols = spectralCols;
pixelSizeX_um = double(scanX_nm) / max(surveyCols, 1) / 1000;
pixelSizeY_um = double(scanY_nm) / max(surveyRows, 1) / 1000;
if ~isfinite(pixelSizeX_um) || pixelSizeX_um <= 0
    pixelSizeX_um = 0.024132355839844;
end
if ~isfinite(pixelSizeY_um) || pixelSizeY_um <= 0
    pixelSizeY_um = 0.024132355839844;
end
wavelengthAxis = double(wavelengthAxis(:));
nWavelengths = numel(wavelengthAxis);
if nWavelengths < 1
    error('precies:workflow:EmptyWavelengthAxis', 'No valid wavelength axis was found.');
end
mappingData = double(mappingData);
if size(mappingData, 1) ~= surveyRows || size(mappingData, 2) ~= surveyCols
    mappingData = imresize(mappingData, [surveyRows, surveyCols]);
end
wavelengthIntensityMaps = zeros(spectralRows, spectralCols, nWavelengths, 'uint16');
globalSpectralMax = 0;
for rowIdx = 1:spectralRows
    for colIdx = 1:spectralCols
        spectrum = spectraData{rowIdx, colIdx};
        if ~isempty(spectrum) && size(spectrum, 2) >= 2
            globalSpectralMax = max(globalSpectralMax, max(double(spectrum(:, 2))));
        end
    end
end
if globalSpectralMax <= 0
    globalSpectralMax = 1;
end
spectralScalingFactor = min(1, 65535 / globalSpectralMax);
for rowIdx = 1:spectralRows
    for colIdx = 1:spectralCols
        spectrum = spectraData{rowIdx, colIdx};
        if ~isempty(spectrum) && size(spectrum, 1) == nWavelengths && size(spectrum, 2) >= 2
            intensities = double(spectrum(:, 2)) * spectralScalingFactor;
            wavelengthIntensityMaps(rowIdx, colIdx, :) = uint16(max(0, min(65535, intensities)));
        end
    end
end
currentDate = char(datetime("now", "Format", "yyyy-MM-dd'T'HH:mm:ss"));
omeXML = buildOdemisOmeXml(currentDate, pixelSizeX_um, pixelSizeY_um, surveyRows, surveyCols, spectralRows, spectralCols, nWavelengths, wavelengthAxis);
t = Tiff(fullpath, 'w8');
cleanupTiff = onCleanup(@() closeTiffIfOpen(t));
tagstruct = struct();
tagstruct.ImageLength = surveyRows;
tagstruct.ImageWidth = surveyCols;
tagstruct.Photometric = Tiff.Photometric.RGB;
tagstruct.BitsPerSample = 8;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.SamplesPerPixel = 3;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'Odemis 3.2-38-gd05792c6a-dirty';
tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Inch;
tagstruct.ImageDescription = omeXML;
tagstruct.SubFileType = 1;
tagstruct.Compression = Tiff.Compression.LZW;
tagstruct.PageName = 'Composited image';
tagstruct.DateTime = char(datetime("now", "Format", "yyyy:MM:dd HH:mm:ss"));
t.setTag(tagstruct);
previewData = repmat(scaleToUint8(mappingData), 1, 1, 3);
t.write(previewData);
t.writeDirectory();
tagstruct.SamplesPerPixel = 1;
tagstruct.BitsPerSample = 16;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.XResolution = 434290;
tagstruct.YResolution = 434290;
tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Centimeter;
tagstruct.SubFileType = 0;
tagstruct.Make = 'pcie-6251';
tagstruct.Model = 'Unknown (driver 3.2-38-gd05792c6a-dirty (driver ni_pcimio v0.7.76, linux 5.4.0))';
tagstruct.PageName = 'Secondary electrons survey';
tagstruct.XPosition = 99.943140064710590;
tagstruct.YPosition = 98.646753643830962;
tagstruct.Software = 'Odemis 3.2-38-gd05792c6a-dirty';
tagstruct.DateTime = char(datetime("now", "Format", "yyyy:MM:dd HH:mm:ss"));
t.setTag(tagstruct);
mappingData16 = scaleToParentUint16(mappingData);
t.write(mappingData16);
t.writeDirectory();
tagstruct.PageName = 'Secondary electrons concurrent';
t.setTag(tagstruct);
t.write(mappingData16);
tagstruct.XResolution = 500000;
tagstruct.YResolution = 500000;
tagstruct.ImageLength = spectralRows;
tagstruct.ImageWidth = spectralCols;
tagstruct.PageName = 'Spectrum with Spectrometer';
tagstruct.Make = 'Andor Newton DU920P_BU2 (s/n: 26340)';
tagstruct.Model = 'PCB: 0/0, firmware: 20.24, EPROM: 0/0 (driver driver: ''0.0.0.0'', SDK: ''2.104.30000.0'')';
tagstruct.XPosition = 99.942719122124387;
tagstruct.YPosition = 98.646638092380442;
for wavelengthIdx = 1:nWavelengths
    t.writeDirectory();
    t.setTag(tagstruct);
    t.write(wavelengthIntensityMaps(:, :, wavelengthIdx));
end
t.close();
clear cleanupTiff;
fprintf('Odemis-compatible TIFF saved: %s (%s)\n', fullpath, exportLabel);
end

function closeTiffIfOpen(t)
try
    t.close();
catch
end
end

function data8 = scaleToUint8(data)
data = double(data);
data(~isfinite(data)) = 0;
maxValue = max(data(:));
if maxValue <= 0
    data8 = zeros(size(data), 'uint8');
else
    data8 = uint8(max(0, min(255, data ./ maxValue * 255)));
end
end

function data16 = scaleToParentUint16(data)
data = double(data);
data(~isfinite(data)) = 0;
maxValue = max(data(:));
if maxValue > 65535
    data = data * (65535 / maxValue);
elseif maxValue == 0
    data(:) = 0;
end
data16 = uint16(max(0, min(65535, data)));
end

function omeXML = buildOdemisOmeXml(currentDate, pixelSizeX_um, pixelSizeY_um, surveyRows, surveyCols, spectralRows, spectralCols, nWavelengths, wavelengthAxis)
extraSettingsStr = ['{"SEM E-beam full": {"accelVoltage": 4000.0, "blanker": true, ' ...
    '"dwellTime": 1.02e-05, "external": true, "horizontalFoV": ' num2str(surveyCols * pixelSizeX_um * 1e-6) ', ' ...
    '"magnification": 22337.748688007505, "pixelSize": [' num2str(pixelSizeX_um * 1e-6) ', ' num2str(pixelSizeY_um * 1e-6) '], ' ...
    '"power": 1, "probeCurrent": 8.9e-11, "resolution": [' num2str(surveyCols) ', ' num2str(surveyRows) '], "rotation": 0.0, ' ...
    '"scale": [8, 8], "translation": [0, 0]}, ' ...
    '"Optical Path Properties": {"focusDistance": 0.0005, "holeDiameter": 0.0006, ' ...
    '"magnification": 0.35, "numericalAperture": 0.2, "parabolaF": 0.0025, ' ...
    '"polePosition": [1641, 971], "refractiveIndex": 1.0, "rotation": 4.712388980384587, ' ...
    '"xMax": 0.01325}, ' ...
    '"Spectrometer": {"exposureTime": 0.2, "binning": [1, 255], "resolution": [' num2str(nWavelengths) ', 1]}}'];
omeXML = ['<?xml version="1.0" encoding="UTF-8"?>' ...
    '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2012-06" ' ...
    'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ' ...
    'xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2012-06 ' ...
    'http://www.openmicroscopy.org/Schemas/OME/2012-06/ome.xsd">' ...
    '<Instrument ID="Instrument:0">' ...
    '<Microscope Manufacturer="Delmic" Model="SPARC" />' ...
    '<Detector ID="Detector:1" Model="pcie-6251" />' ...
    '<Objective CalibratedMagnification="22337.748688007504825" ID="Objective:1" />' ...
    '<Detector ID="Detector:2" Model="pcie-6251" />' ...
    '<Objective CalibratedMagnification="22337.748688007504825" ID="Objective:2" />' ...
    '<Detector ID="Detector:3" Model="Andor Newton DU920P_BU2 (s/n: 26340)" />' ...
    '<Objective CalibratedMagnification="0.350000000000000" ID="Objective:3" Model="Optical Path Properties" />' ...
    '</Instrument>' ...
    '<Image ID="Image:0" Name="Composited image preview">' ...
    '<AcquisitionDate>' currentDate '</AcquisitionDate>' ...
    '<Transform A00="1.000000000000000" A01="0.000000000000000" A02="0.000000000000000" ' ...
    'A10="0.000000000000000" A11="1.000000000000000" A12="0.000000000000000" />' ...
    '<Pixels DimensionOrder="CXYZT" ID="Pixels:0" ' ...
    'PhysicalSizeX="' num2str(pixelSizeX_um) '" ' ...
    'PhysicalSizeY="' num2str(pixelSizeY_um) '" ' ...
    'SizeC="3" SizeT="1" SizeX="' num2str(surveyCols) '" ' ...
    'SizeY="' num2str(surveyRows) '" SizeZ="1" Type="uint8">' ...
    '<Channel ID="Channel:0:0" Name="Composited image preview" SamplesPerPixel="1" />' ...
    '<TiffData FirstC="0" FirstT="0" FirstZ="0" IFD="0" PlaneCount="1" />' ...
    '<Plane PositionX="0" PositionY="0" TheC="0" TheT="0" TheZ="0" />' ...
    '</Pixels>' ...
    '</Image>' ...
    '<Image ID="Image:1" Name="Secondary electrons survey">' ...
    '<AcquisitionDate>' currentDate '</AcquisitionDate>' ...
    '<ObjectiveSettings ID="Objective:1" />' ...
    '<ExtraSettings>' extraSettingsStr '</ExtraSettings>' ...
    '<Pixels DimensionOrder="XYZTC" ID="Pixels:1" ' ...
    'PhysicalSizeX="' num2str(pixelSizeX_um) '" ' ...
    'PhysicalSizeY="' num2str(pixelSizeY_um) '" ' ...
    'SizeC="1" SizeT="1" SizeX="' num2str(surveyCols) '" ' ...
    'SizeY="' num2str(surveyRows) '" SizeZ="1" Type="uint16">' ...
    '<Channel ID="Channel:1:0" Name="Secondary electrons survey" />' ...
    '<TiffData FirstC="0" FirstT="0" FirstZ="0" IFD="1" PlaneCount="1" />' ...
    '<Plane ExposureTime="0.000010200000000" IntegrationCount="12" PositionX="0" PositionY="0" TheC="0" TheT="0" TheZ="0" />' ...
    '</Pixels>' ...
    '</Image>' ...
    '<Image ID="Image:2" Name="Secondary electrons concurrent">' ...
    '<AcquisitionDate>' currentDate '</AcquisitionDate>' ...
    '<ObjectiveSettings ID="Objective:2" />' ...
    '<ExtraSettings>' extraSettingsStr '</ExtraSettings>' ...
    '<Pixels DimensionOrder="XYZTC" ID="Pixels:2" ' ...
    'PhysicalSizeX="' num2str(pixelSizeX_um) '" ' ...
    'PhysicalSizeY="' num2str(pixelSizeY_um) '" ' ...
    'SizeC="1" SizeT="1" SizeX="' num2str(spectralCols) '" ' ...
    'SizeY="' num2str(spectralRows) '" SizeZ="' num2str(nWavelengths) '" Type="uint16">' ...
    '<Channel ID="Channel:2:0" Name="Spectrum with Spectrometer" EmissionWavelength="' num2str(mean(wavelengthAxis)) ...
    '" SpectrumRange="' num2str(min(wavelengthAxis)) '-' num2str(max(wavelengthAxis)) ' nm" />' ...
    '<TiffData FirstC="0" FirstT="0" FirstZ="0" IFD="2" PlaneCount="1" />' ...
    '<Plane ExposureTime="0.2" IntegrationCount="1" PositionX="0" PositionY="0" TheC="0" TheT="0" TheZ="0" />' ...
    '</Pixels>' ...
    '</Image>' ...
    '</OME>'];
end

function [axisNm, cube] = spectraToCube(spectra)
[rows, cols] = size(spectra);
axisNm = [];
for rowIdx = 1:rows
    for colIdx = 1:cols
        spectrum = spectra{rowIdx, colIdx};
        if ~isempty(spectrum) && size(spectrum, 2) >= 2
            axisNm = double(spectrum(:, 1));
            break;
        end
    end
    if ~isempty(axisNm)
        break;
    end
end
if isempty(axisNm)
    error('precies:workflow:EmptySpectra', 'No valid spectra were found in the simulation output.');
end
cube = zeros(rows, cols, numel(axisNm));
for rowIdx = 1:rows
    for colIdx = 1:cols
        spectrum = spectra{rowIdx, colIdx};
        if isempty(spectrum) || size(spectrum, 2) < 2
            continue;
        end
        if numel(spectrum(:, 1)) == numel(axisNm) && max(abs(double(spectrum(:, 1)) - axisNm)) < 1e-9
            cube(rowIdx, colIdx, :) = reshape(double(spectrum(:, 2)), 1, 1, []);
        else
            cube(rowIdx, colIdx, :) = reshape(interp1(double(spectrum(:, 1)), double(spectrum(:, 2)), axisNm, 'linear', 0), 1, 1, []);
        end
    end
end
end

function table1 = buildTable1Validation(summaryTable)
voltage = getColumn(summaryTable, 'Voltage_kV', NaN);
averageNcc = getColumn(summaryTable, 'Average_NCC', NaN);
sdNcc = getColumn(summaryTable, 'SD_NCC', NaN);
nrmsePercent = 100 * getColumn(summaryTable, 'nRMSE', NaN);
sdNrmsePercent = 100 * getColumn(summaryTable, 'SD_nRMSE', NaN);
deltaLambdaMqw = getFirstColumn(summaryTable, {'AbsDeltaLambdaMQW_nm', 'Abs_DeltaLambda_MQW_nm'}, NaN);
deltaLambdaGan = getFirstColumn(summaryTable, {'AbsDeltaLambdaGaN_nm', 'Abs_DeltaLambda_GaN_nm'}, NaN);
deltaFwhmMqw = getFirstColumn(summaryTable, {'AbsDeltaFWHMMQW_nm', 'Abs_DeltaFWHM_MQW_nm'}, NaN);
table1 = table(voltage, averageNcc, sdNcc, nrmsePercent, sdNrmsePercent, ...
    deltaLambdaMqw, deltaLambdaGan, deltaFwhmMqw, ...
    'VariableNames', {'Voltage_kV', 'Average_NCC', 'SD_NCC', 'nRMSE_percent', ...
    'SD_nRMSE_percent', 'AbsDeltaLambdaMQW_nm', 'AbsDeltaLambdaGaN_nm', ...
    'AbsDeltaFWHMMQW_nm'});
end

function values = getFirstColumn(inputTable, names, defaultValue)
values = defaultValue * ones(height(inputTable), 1);
for idx = 1:numel(names)
    if ismember(names{idx}, inputTable.Properties.VariableNames)
        values = double(inputTable.(names{idx}));
        return;
    end
end
end

function values = getColumn(inputTable, name, defaultValue)
if ismember(name, inputTable.Properties.VariableNames)
    values = double(inputTable.(name));
else
    values = defaultValue * ones(height(inputTable), 1);
end
end

function activityPackage = exportActivityStructure(simulationData, importedTiffData, params, outputPath)
[X, Y, Z, V, W] = collectPoints(simulationData);
if isempty(X)
    error('precies:workflow:NoActivityPoints', 'No source point cloud was found in the simulation output.');
end
layerTable = buildLayerTable(simulationData);
vpitTable = buildVpitTable(simulationData);
facetTable = buildFacetTable(simulationData);
gridPoints = max(20, round(params.gridPoints3D));
xGrid = linspace(-params.scanX / 2, params.scanX / 2, gridPoints);
yGrid = linspace(-params.scanY / 2, params.scanY / 2, gridPoints);
[activityMaps, activityMapTable] = buildActivityMaps(X, Y, V, W, xGrid, yGrid);
metadata = table( ...
    ["ExportType"; "Source"; "Voltage_kV"; "PointCount"; "ScanX_nm"; "ScanY_nm"; "GridPoints"; "VPitCount"], ...
    ["Activity-Structure input package"; "run_reconstruction"; string(params.electronEnergy / 1e3); string(numel(V)); string(params.scanX); string(params.scanY); string(gridPoints); string(height(vpitTable))], ...
    'VariableNames', {'Key', 'Value'});
activityPackage = struct( ...
    'metadata', metadata, ...
    'points', table(X(:), Y(:), Z(:), V(:), W(:), 'VariableNames', {'X_nm', 'Y_nm', 'Z_nm', 'Intensity', 'Wavelength_nm'}), ...
    'layers', layerTable, ...
    'vPits', vpitTable, ...
    'vPitFacets', facetTable, ...
    'xGrid_nm', xGrid, ...
    'yGrid_nm', yGrid, ...
    'activityMaps', activityMaps, ...
    'activityMapTable', activityMapTable, ...
    'importedTiffData', importedTiffData, ...
    'description', 'Activity-Structure input package exported by the reviewer reproduction workflow.');
save(outputPath, 'activityPackage', '-v7.3');
end

function [X, Y, Z, V, W] = collectPoints(simulationData)
positions = simulationData.positions;
intensities = simulationData.intensities;
X = positions(:, 1) * 1e9;
Y = positions(:, 2) * 1e9;
Z = positions(:, 3) * 1e9;
V = intensities(:);
count = min([numel(X), numel(Y), numel(Z), numel(V)]);
X = X(1:count);
Y = Y(1:count);
Z = Z(1:count);
V = V(1:count);
W = zeros(count, 1);
if isfield(simulationData, 'wavelengths') && ~isempty(simulationData.wavelengths)
    Wraw = simulationData.wavelengths(:);
    if max(Wraw, [], 'omitnan') < 10
        Wraw = Wraw * 1e9;
    end
    if numel(Wraw) >= count
        W = Wraw(1:count);
    end
end
valid = isfinite(X) & isfinite(Y) & isfinite(Z) & isfinite(V) & V > 0;
X = X(valid);
Y = Y(valid);
Z = Z(valid);
V = V(valid);
W = W(valid);
end

function layerTable = buildLayerTable(simulationData)
names = string(simulationData.layerNames(:));
boundaries = simulationData.layerBoundaries(:) * 1e9;
finiteCount = min(numel(names), max(0, numel(boundaries) - 1));
names = names(1:finiteCount);
zBottom = boundaries(1:finiteCount);
zTop = boundaries(2:finiteCount + 1);
thickness = zTop - zBottom;
layerIndex = (1:finiteCount)';
channel = strings(finiteCount, 1);
weight = zeros(finiteCount, 1);
nominalWavelength = nan(finiteCount, 1);
for idx = 1:finiteCount
    [channel(idx), weight(idx), nominalWavelength(idx)] = classifyLayer(names(idx));
end
layerTable = table(layerIndex, names, zBottom, zTop, thickness, channel, weight, nominalWavelength, ...
    'VariableNames', {'LayerIndex', 'LayerName', 'ZBottom_nm', 'ZTop_nm', 'Thickness_nm', ...
    'ActivityChannel', 'ActivityWeight', 'NominalWavelength_nm'});
end

function [channel, weight, wavelength] = classifyLayer(layerName)
text = char(layerName);
if contains(text, 'MQW-Well')
    channel = "MQW_470_570";
    weight = 1.00;
    wavelength = 515;
elseif contains(text, 'MQW-Barrier')
    channel = "Shoulder_405_470";
    weight = 0.18;
    wavelength = 445;
elseif contains(text, 'Prestrained')
    channel = "Shoulder_405_470";
    weight = 0.09;
    wavelength = 420;
elseif contains(text, 'p-EBL')
    channel = "Shortwave_370_405";
    weight = 0.06;
    wavelength = 405;
elseif contains(text, 'GaN')
    channel = "Shortwave_370_405";
    weight = 0.06;
    wavelength = 390;
else
    channel = "Inactive";
    weight = 0;
    wavelength = NaN;
end
end

function vpitTable = buildVpitTable(simulationData)
vpitTable = table([], [], [], [], [], [], [], [], ...
    'VariableNames', {'VPitIndex', 'CenterX_nm', 'CenterY_nm', 'TopRadius_nm', ...
    'Depth_nm', 'Orientation_deg', 'SurfaceZ_nm', 'ApexZ_nm'});
if ~isfield(simulationData, 'Vpits') || isempty(simulationData.Vpits)
    return;
end
rows = zeros(numel(simulationData.Vpits), 8);
for idx = 1:numel(simulationData.Vpits)
    pit = simulationData.Vpits{idx};
    orientation = fieldOrDefault(pit, 'orientationDeg', NaN);
    rows(idx, :) = [idx, pit.center(1), pit.center(2), pit.topRadius, pit.depth, orientation, pit.surfaceZ_nm, pit.apex(3)];
end
vpitTable = array2table(rows, 'VariableNames', vpitTable.Properties.VariableNames);
end

function facetTable = buildFacetTable(simulationData)
facetTable = table([], [], [], [], [], [], [], [], [], [], [], [], [], [], ...
    'VariableNames', {'VPitIndex', 'FacetIndex', 'StartX_nm', 'StartY_nm', 'StartZ_nm', ...
    'EndX_nm', 'EndY_nm', 'EndZ_nm', 'ApexX_nm', 'ApexY_nm', 'ApexZ_nm', ...
    'NormalX', 'NormalY', 'NormalZ'});
if ~isfield(simulationData, 'Vpits') || isempty(simulationData.Vpits)
    return;
end
rows = zeros(0, 14);
for idx = 1:numel(simulationData.Vpits)
    pit = simulationData.Vpits{idx};
    if ~isfield(pit, 'facets')
        continue;
    end
    for facetIdx = 1:numel(pit.facets)
        facet = pit.facets(facetIdx);
        rows(end + 1, :) = [idx, facetIdx, facet.vertexStart, facet.vertexEnd, pit.apex, facet.planeNormalOutward];
    end
end
if ~isempty(rows)
    facetTable = array2table(rows, 'VariableNames', facetTable.Properties.VariableNames);
end
end

function [maps, mapTable] = buildActivityMaps(X, Y, V, W, xGrid, yGrid)
valid = isfinite(X) & isfinite(Y) & isfinite(V) & V > 0;
allMap = gridActivityMap(X, Y, V, xGrid, yGrid, valid);
hasWavelength = ~isempty(W) && numel(W) == numel(V) && any(isfinite(W) & W > 0);
if hasWavelength
    mqwMap = gridActivityMap(X, Y, V, xGrid, yGrid, valid & W >= 470 & W <= 570);
    shoulderMap = gridActivityMap(X, Y, V, xGrid, yGrid, valid & W >= 405 & W < 470);
    shortMap = gridActivityMap(X, Y, V, xGrid, yGrid, valid & W >= 370 & W < 405);
else
    mqwMap = allMap;
    shoulderMap = 0.25 * allMap;
    shortMap = 0.20 * allMap;
end
maps = struct('all', normalizeMap(allMap), 'mqw', normalizeMap(max(mqwMap, 0.65 * allMap)), ...
    'shoulder', normalizeMap(max(shoulderMap, 0.22 * allMap)), 'short', normalizeMap(max(shortMap, 0.20 * allMap)));
[Xmesh, Ymesh] = meshgrid(xGrid, yGrid);
mapTable = table(Xmesh(:), Ymesh(:), maps.all(:), maps.mqw(:), maps.shoulder(:), maps.short(:), ...
    'VariableNames', {'X_nm', 'Y_nm', 'AllActivity', 'MQW_470_570_Activity', ...
    'Shoulder_405_470_Activity', 'Shortwave_370_405_Activity'});
end

function map = gridActivityMap(X, Y, V, xGrid, yGrid, mask)
map = zeros(numel(yGrid), numel(xGrid));
if nnz(mask) < 1
    return;
end
xEdges = gridEdges(xGrid);
yEdges = gridEdges(yGrid);
ix = discretize(X(mask), xEdges);
iy = discretize(Y(mask), yEdges);
values = V(mask);
valid = isfinite(ix) & isfinite(iy) & isfinite(values);
if ~any(valid)
    return;
end
map = accumarray([iy(valid), ix(valid)], values(valid), size(map), @sum, 0);
kernel = gaussianKernel(2.0);
map = conv2(map, kernel, 'same');
end

function edges = gridEdges(grid)
step = median(diff(grid));
edges = [grid(1) - 0.5 * step, 0.5 * (grid(1:end-1) + grid(2:end)), grid(end) + 0.5 * step];
end

function kernel = gaussianKernel(sigma)
radius = max(2, ceil(3 * sigma));
[yy, xx] = ndgrid(-radius:radius, -radius:radius);
kernel = exp(-0.5 * ((xx / sigma) .^ 2 + (yy / sigma) .^ 2));
kernel = kernel / max(sum(kernel(:)), eps);
end

function map = normalizeMap(map)
map = max(map, 0);
values = map(isfinite(map) & map > 0);
if isempty(values)
    return;
end
scale = prctile(values, 97);
if ~isfinite(scale) || scale <= 0
    scale = max(values);
end
map = min(map ./ max(scale, eps), 1.35);
end

function value = fieldOrDefault(data, name, defaultValue)
if isstruct(data) && isfield(data, name) && ~isempty(data.(name))
    value = data.(name);
else
    value = defaultValue;
end
end

function table2 = buildTable2Ablation(summaryTable)
scenario = string(summaryTable.Scenario);
models = ["M0"; "M1"; "M2"; "M3"];
pe = ["Mean-depth approximation"; "Monte Carlo PE interaction volume"; ...
    "Monte Carlo PE interaction volume"; "Monte Carlo PE interaction volume"];
photon = ["No"; "No"; "Yes"; "Yes"];
secondary = ["No"; "No"; "No"; "Yes"];
objective = nan(numel(models), 1);
averageNcc = nan(numel(models), 1);
nrmsePercent = nan(numel(models), 1);
deltaLambdaMqw = nan(numel(models), 1);
deltaFwhmMqw = nan(numel(models), 1);
for idx = 1:numel(models)
    rowIdx = find(strcmpi(scenario, models(idx)), 1);
    if isempty(rowIdx)
        continue;
    end
    objective(idx) = summaryTable.Objective(rowIdx);
    averageNcc(idx) = summaryTable.Average_NCC(rowIdx);
    nrmsePercent(idx) = 100 * summaryTable.nRMSE(rowIdx);
    deltaLambdaMqw(idx) = summaryTable.AbsDeltaLambdaMQW_nm(rowIdx);
    deltaFwhmMqw(idx) = summaryTable.AbsDeltaFWHMMQW_nm(rowIdx);
end
table2 = table(models, pe, photon, secondary, objective, averageNcc, nrmsePercent, deltaLambdaMqw, deltaFwhmMqw, ...
    'VariableNames', {'Model', 'PEInteractionVolume', 'PhotonTransport', 'SecondaryExcitation', ...
    'CompositeObjective', 'AverageNCC', 'nRMSE_percent', 'AbsDeltaLambdaMQW_nm', 'AbsDeltaFWHMMQW_nm'});
end

function files = writeTable2Files(table2, outputDir)
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
files = struct();
files.csv = fullfile(outputDir, 'Table2_forward_ablation_metrics.csv');
files.xlsx = fullfile(outputDir, 'Table2_forward_ablation_metrics.xlsx');
writetable(table2, files.csv);
writetable(table2, files.xlsx, 'Sheet', 'Table 2');
end

function activityPackage = loadActivityPackage(matFile)
if ~isfile(matFile)
    error('precies:workflow:MissingMatFile', 'MAT file does not exist: %s', matFile);
end
loaded = load(matFile);
if isfield(loaded, 'activityPackage')
    activityPackage = loaded.activityPackage;
    return;
end
names = fieldnames(loaded);
for idx = 1:numel(names)
    value = loaded.(names{idx});
    if isstruct(value) && isfield(value, 'points') && isfield(value, 'layers')
        activityPackage = value;
        return;
    end
end
error('precies:workflow:InvalidActivityPackage', 'The selected MAT file does not contain an Activity-Structure package.');
end

function [X, Y, Z, V, W] = extractActivityPoints(activityPackage)
pointTable = activityPackage.points;
names = lower(string(pointTable.Properties.VariableNames));
xCol = findColumn(names, ["x_nm", "x"]);
yCol = findColumn(names, ["y_nm", "y"]);
zCol = findColumn(names, ["z_nm", "z", "depth"]);
vCol = findColumn(names, ["intensity", "value", "signal"]);
wCol = findColumn(names, ["wavelength_nm", "wavelength", "lambda"]);
if any([xCol, yCol, zCol, vCol] == 0)
    error('precies:workflow:InvalidPointTable', 'Activity points must contain X, Y, Z, and Intensity columns.');
end
X = double(pointTable{:, xCol});
Y = double(pointTable{:, yCol});
Z = double(pointTable{:, zCol});
V = double(pointTable{:, vCol});
count = min([numel(X), numel(Y), numel(Z), numel(V)]);
X = X(1:count);
Y = Y(1:count);
Z = Z(1:count);
V = V(1:count);
if wCol > 0
    W = double(pointTable{:, wCol});
    if numel(W) >= count
        W = W(1:count);
    else
        W = zeros(count, 1);
    end
else
    W = zeros(count, 1);
end
if max(W, [], 'omitnan') < 10
    W = W * 1e9;
end
valid = isfinite(X) & isfinite(Y) & isfinite(Z) & isfinite(V);
X = X(valid);
Y = Y(valid);
Z = Z(valid);
V = V(valid);
W = W(valid);
end

function idx = findColumn(names, tokens)
idx = 0;
for tokenIdx = 1:numel(tokens)
    token = lower(string(tokens(tokenIdx)));
    match = find(names == token, 1);
    if ~isempty(match)
        idx = match;
        return;
    end
end
for tokenIdx = 1:numel(tokens)
    token = lower(string(tokens(tokenIdx)));
    match = find(contains(names, token), 1);
    if ~isempty(match)
        idx = match;
        return;
    end
end
end

function volumeData = buildLegacyVolume(activityPackage, gridPoints)
if nargin < 2 || isempty(gridPoints)
    gridPoints = 90;
end
[X, Y, Z, V, W] = extractActivityPoints(activityPackage);
valid = isfinite(X) & isfinite(Y) & isfinite(Z) & isfinite(V) & V > 0;
X = X(valid);
Y = Y(valid);
Z = Z(valid);
V = V(valid);
W = W(valid);
[X, Y, Z, V, W] = limitPointCloudSamples(X, Y, Z, V, W);
[X, Y, Z, V, W, depthRange] = applyActiveStackDepthFilter(activityPackage, X, Y, Z, V, W);
if numel(V) < 4
    error('precies:workflow:InsufficientPointCloud', 'At least four valid source points are required.');
end
[xGrid, yGrid] = resolvePackageLateralGrids(activityPackage, X, Y, gridPoints);
zGrid = linspace(depthRange.zMin_nm, depthRange.zMax_nm, 10 * gridPoints);
[Xg, Yg, Zg] = meshgrid(xGrid, yGrid, zGrid);
intensityInterpolant = scatteredInterpolant(X(:), Y(:), Z(:), V(:), 'natural', 'linear');
Vq = intensityInterpolant(Xg, Yg, Zg);
VForLog = V(:);
positiveValues = VForLog(VForLog > 1.1);
if isempty(positiveValues)
    minPositive = 1.1;
else
    minPositive = min(positiveValues);
end
VForLog(VForLog <= 1.1) = minPositive;
logInterpolant = scatteredInterpolant(X(:), Y(:), Z(:), log10(VForLog), 'natural', 'linear');
logVolume = logInterpolant(Xg, Yg, Zg);
wavelengthVolume = buildPointWavelengthVolume(X, Y, Z, W, Xg, Yg, Zg, Vq);
volumeData = struct( ...
    'xGrid', xGrid, ...
    'yGrid', yGrid, ...
    'zGrid', zGrid, ...
    'volume', Vq, ...
    'intensityVolume', Vq, ...
    'logVolume', logVolume, ...
    'wavelengthVolume', wavelengthVolume, ...
    'Vq_raw', Vq, ...
    'Vq_log', logVolume, ...
    'Vq_wavelength', wavelengthVolume, ...
    'points', table(X(:), Y(:), Z(:), V(:), W(:), 'VariableNames', {'X_nm', 'Y_nm', 'Z_nm', 'Intensity', 'Wavelength_nm'}), ...
    'activityPackage', activityPackage, ...
    'depthRange', depthRange, ...
    'method', 'Legacy Grid');
end

function volumeData = buildActivityStructureVolume(activityPackage, gridPoints)
if nargin < 2 || isempty(gridPoints)
    gridPoints = 90;
end
[X, Y, Z, V, W] = extractActivityPoints(activityPackage);
[X, Y, Z, V, W] = limitPointCloudSamples(X, Y, Z, V, W);
[X, Y, Z, V, W, depthRange] = applyActiveStackDepthFilter(activityPackage, X, Y, Z, V, W);
[xGrid, yGrid] = resolvePackageLateralGrids(activityPackage, X, Y, gridPoints);
zGrid = linspace(depthRange.zMin_nm, depthRange.zMax_nm, 10 * gridPoints);
[layers, layerBoundaries] = layerTableToLayerContext(activityPackage.layers);
layerInfos = buildActivityStructureLayerInfosFromContext(layers, layerBoundaries);
if isempty(layerInfos)
    legacy = buildLegacyVolume(activityPackage, gridPoints);
    volume = legacy.volume;
    logVolume = legacy.logVolume;
    wavelengthVolume = legacy.wavelengthVolume;
    xGrid = legacy.xGrid;
    yGrid = legacy.yGrid;
    zGrid = legacy.zGrid;
else
    hasWavelengthData = ~isempty(W) && any(isfinite(W) & W > 0);
    baseMaps = buildActivityChannelMaps(X, Y, Z, V, W, xGrid, yGrid, layerInfos, hasWavelengthData);
    vPitData = buildVpitDataFromPackage(activityPackage, layers, layerBoundaries);
    [volume, wavelengthVolume] = buildActivityStructureVolumeFromParentLogic(layerInfos, baseMaps, ...
        xGrid, yGrid, zGrid, hasWavelengthData, vPitData);
    logVolume = buildLogVolume(volume);
end
volumeData = struct( ...
    'xGrid', xGrid, ...
    'yGrid', yGrid, ...
    'zGrid', zGrid, ...
    'volume', volume, ...
    'intensityVolume', volume, ...
    'logVolume', logVolume, ...
    'wavelengthVolume', wavelengthVolume, ...
    'Vq_raw', volume, ...
    'Vq_log', logVolume, ...
    'Vq_wavelength', wavelengthVolume, ...
    'activityPackage', activityPackage, ...
    'vPits', getPackageTable(activityPackage, 'vPits'), ...
    'layers', activityPackage.layers, ...
    'layerInfos', layerInfos, ...
    'depthRange', depthRange, ...
    'method', 'Activity-Structure');
end

function [X, Y, Z, V, W] = limitPointCloudSamples(X, Y, Z, V, W)
maxSamplePoints = 50000;
if numel(X) <= maxSamplePoints
    return;
end
sampleIndices = randperm(numel(X), maxSamplePoints);
X = X(sampleIndices);
Y = Y(sampleIndices);
Z = Z(sampleIndices);
V = V(sampleIndices);
if ~isempty(W)
    W = W(sampleIndices);
end
end

function [xGrid, yGrid] = resolvePackageLateralGrids(activityPackage, X, Y, gridPoints)
defaultScanX = 500;
defaultScanY = 500;
xGrid = linspace(-defaultScanX / 2, defaultScanX / 2, gridPoints);
yGrid = linspace(-defaultScanY / 2, defaultScanY / 2, gridPoints);
end

function value = getActivityPackageNumeric(activityPackage, keyName)
value = NaN;
if ~isstruct(activityPackage)
    return;
end
fieldNames = fieldnames(activityPackage);
match = find(strcmpi(fieldNames, keyName), 1);
if ~isempty(match)
    rawValue = activityPackage.(fieldNames{match});
    parsed = parseNumericScalar(rawValue);
    if isfinite(parsed)
        value = parsed;
        return;
    end
end
if ~isfield(activityPackage, 'metadata') || ~istable(activityPackage.metadata) || height(activityPackage.metadata) < 1
    return;
end
metadata = activityPackage.metadata;
names = lower(string(metadata.Properties.VariableNames));
keyCol = findColumn(names, ["key", "name", "parameter"]);
valueCol = findColumn(names, "value");
if keyCol == 0 || valueCol == 0
    return;
end
keys = string(metadata{:, keyCol});
rowIdx = find(strcmpi(keys, keyName), 1);
if isempty(rowIdx)
    return;
end
parsed = parseNumericScalar(metadata{rowIdx, valueCol});
if isfinite(parsed)
    value = parsed;
end
end

function value = parseNumericScalar(rawValue)
value = NaN;
if isnumeric(rawValue) && ~isempty(rawValue)
    rawValue = rawValue(:);
    value = double(rawValue(1));
    return;
end
parsedValues = str2double(string(rawValue));
parsedValues = parsedValues(isfinite(parsedValues));
if ~isempty(parsedValues)
    value = parsedValues(1);
end
end

function [layers, layerBoundaries] = layerTableToLayerContext(layerTable)
layers = {};
layerBoundaries = [];
if isempty(layerTable) || ~istable(layerTable) || height(layerTable) < 1
    return;
end
names = lower(string(layerTable.Properties.VariableNames));
nameCol = findColumn(names, ["layername", "name"]);
bottomCol = findColumn(names, ["zbottom_nm", "bottom"]);
topCol = findColumn(names, ["ztop_nm", "top"]);
thicknessCol = findColumn(names, ["thickness_nm", "thickness"]);
if any([nameCol, bottomCol, topCol] == 0)
    return;
end
layerNames = string(layerTable{:, nameCol});
zBottom = double(layerTable{:, bottomCol});
zTop = double(layerTable{:, topCol});
valid = isfinite(zBottom) & isfinite(zTop) & zTop > zBottom;
layerNames = layerNames(valid);
zBottom = zBottom(valid);
zTop = zTop(valid);
if isempty(zBottom)
    return;
end
[zBottom, order] = sort(zBottom);
zTop = zTop(order);
layerNames = layerNames(order);
if thicknessCol > 0
    thickness = double(layerTable{valid, thicknessCol});
    thickness = thickness(order);
else
    thickness = zTop - zBottom;
end
layers = cell(numel(layerNames), 3);
layers(:, 1) = cellstr(layerNames);
layers(:, 2) = {''};
layers(:, 3) = num2cell(abs(thickness(:)));
layerBoundaries = [zBottom(1); zTop(:)]' * 1e-9;
if ~isinf(layerBoundaries(end))
    layerBoundaries(end + 1) = Inf;
end
end

function layerInfos = buildActivityStructureLayerInfosFromContext(layers, layerBoundaries)
layerInfos = repmat(struct( ...
    'name', '', ...
    'index', 0, ...
    'zBottom_nm', NaN, ...
    'zTop_nm', NaN, ...
    'channel', 'inactive', ...
    'weight', 0, ...
    'nominalWavelength_nm', NaN, ...
    'faceAlpha', 0), 0, 1);
if isempty(layers) || isempty(layerBoundaries)
    return;
end
layerBoundariesNm = layerBoundaries(:) * 1e9;
finiteLayerCount = min(size(layers, 1), numel(layerBoundariesNm) - 1);
for layerIdx = 1:finiteLayerCount
    zBottom = layerBoundariesNm(layerIdx);
    zTop = layerBoundariesNm(layerIdx + 1);
    if ~isfinite(zBottom) || ~isfinite(zTop) || zTop <= zBottom
        continue;
    end
    layerName = layers{layerIdx, 1};
    [channel, weight, nominalWavelength, faceAlpha] = classifyActivityStructureLayer(layerName);
    if weight <= 0
        continue;
    end
    layerInfos(end + 1, 1) = struct( ...
        'name', layerName, ...
        'index', layerIdx, ...
        'zBottom_nm', zBottom, ...
        'zTop_nm', zTop, ...
        'channel', channel, ...
        'weight', weight, ...
        'nominalWavelength_nm', nominalWavelength, ...
        'faceAlpha', faceAlpha);
end
end

function [channel, weight, nominalWavelength, faceAlpha] = classifyActivityStructureLayer(layerName)
text = char(string(layerName));
if contains(text, 'MQW-Well')
    channel = 'mqw';
    weight = 1.00;
    nominalWavelength = 515;
    faceAlpha = 0.34;
elseif contains(text, 'MQW-Barrier')
    channel = 'shoulder';
    weight = 0.18;
    nominalWavelength = 445;
    faceAlpha = 0.13;
elseif contains(text, 'Prestrained')
    channel = 'shoulder';
    weight = 0.09;
    nominalWavelength = 420;
    faceAlpha = 0.10;
elseif contains(text, 'p-EBL')
    channel = 'short';
    weight = 0.06;
    nominalWavelength = 405;
    faceAlpha = 0.08;
elseif contains(text, 'p-GaN')
    channel = 'short';
    weight = 0.05;
    nominalWavelength = 389;
    faceAlpha = 0.07;
elseif contains(text, 'n-GaN')
    channel = 'short';
    weight = 0.08;
    nominalWavelength = 390;
    faceAlpha = 0.08;
else
    channel = 'inactive';
    weight = 0;
    nominalWavelength = NaN;
    faceAlpha = 0;
end
end

function maps = buildActivityChannelMaps(X, Y, Z, V, W, xGrid, yGrid, layerInfos, hasWavelengthData)
layerZMin = min([layerInfos.zBottom_nm]) - 35;
layerZMax = max([layerInfos.zTop_nm]) + 35;
wavelengthNm = W(:);
finiteMask = isfinite(Z) & Z >= layerZMin & Z <= layerZMax & isfinite(V) & V > 0 & isfinite(X) & isfinite(Y);
if hasWavelengthData && ~isempty(W)
    finiteMask = finiteMask & isfinite(wavelengthNm);
else
    wavelengthNm = nan(size(V));
end
allMap = computeActivityMap2D(X, Y, V, xGrid, yGrid, finiteMask);
if hasWavelengthData && ~isempty(W)
    mqwMap = computeActivityMap2D(X, Y, V, xGrid, yGrid, finiteMask & wavelengthNm >= 470 & wavelengthNm <= 570);
    shoulderMap = computeActivityMap2D(X, Y, V, xGrid, yGrid, finiteMask & wavelengthNm >= 405 & wavelengthNm < 470);
    shortMap = computeActivityMap2D(X, Y, V, xGrid, yGrid, finiteMask & wavelengthNm >= 370 & wavelengthNm <= 405);
else
    mqwMap = allMap;
    shoulderMap = 0.25 * allMap;
    shortMap = 0.20 * allMap;
end
mqwMap = fillSparseActivityMap(mqwMap, allMap, 0.65);
shoulderMap = fillSparseActivityMap(shoulderMap, allMap, 0.22);
shortMap = fillSparseActivityMap(shortMap, allMap, 0.20);
maps = struct('all', allMap, 'mqw', mqwMap, 'shoulder', shoulderMap, 'short', shortMap);
end

function map2D = computeActivityMap2D(X, Y, V, xGrid, yGrid, sampleMask)
dims = [numel(yGrid), numel(xGrid)];
map2D = zeros(dims);
if nargin < 6 || isempty(sampleMask) || nnz(sampleMask) < 3
    return;
end
xEdges = buildGridEdges(xGrid);
yEdges = buildGridEdges(yGrid);
ix = discretize(X(sampleMask), xEdges);
iy = discretize(Y(sampleMask), yEdges);
values = V(sampleMask);
valid = isfinite(ix) & isfinite(iy) & isfinite(values) & values > 0;
if nnz(valid) < 3
    return;
end
rawMap = accumarray([iy(valid), ix(valid)], values(valid), dims, @sum, 0);
occupancyMap = accumarray([iy(valid), ix(valid)], ones(nnz(valid), 1), dims, @sum, 0);
xStep = max(median(diff(xGrid)), eps);
yStep = max(median(diff(yGrid)), eps);
kernel = buildGaussianKernel2D(max(1.2, 16 / yStep), max(1.2, 16 / xStep));
supportKernel = buildGaussianKernel2D(max(1.0, 10 / yStep), max(1.0, 10 / xStep));
smoothMap = conv2(rawMap, kernel, 'same');
supportMap = conv2(double(occupancyMap > 0), supportKernel, 'same');
smoothMap(supportMap <= 0.015 * max(supportMap(:), [], 'omitnan')) = 0;
map2D = normalizeActivityMap(smoothMap);
end

function map2D = fillSparseActivityMap(map2D, fallbackMap, fallbackWeight)
if any(map2D(:) > 0)
    map2D = max(map2D, fallbackWeight * fallbackMap);
else
    map2D = fallbackWeight * fallbackMap;
end
map2D = normalizeActivityMap(map2D);
end

function normalizedMap = normalizeActivityMap(map2D)
normalizedMap = max(map2D, 0);
validValues = normalizedMap(isfinite(normalizedMap) & normalizedMap > 0);
if isempty(validValues)
    normalizedMap(:) = 0;
    return;
end
scaleValue = prctile(validValues, 97);
if ~isfinite(scaleValue) || scaleValue <= eps
    scaleValue = max(validValues);
end
normalizedMap = normalizedMap ./ max(scaleValue, eps);
normalizedMap = min(max(normalizedMap, 0), 1.35);
end

function [volume, wavelengthVolume] = buildActivityStructureVolumeFromParentLogic(layerInfos, maps, xGrid, yGrid, zGrid, hasWavelengthData, vPitData)
dims = [numel(yGrid), numel(xGrid), numel(zGrid)];
volume = zeros(dims);
wavelengthNumerator = zeros(dims);
supportMask = false(dims);
for layerIdx = 1:numel(layerInfos)
    info = layerInfos(layerIdx);
    zMask = zGrid >= info.zBottom_nm & zGrid <= info.zTop_nm;
    if ~any(zMask)
        continue;
    end
    layerMap = resolveActivityMapForChannel(maps, info.channel) * info.weight;
    if ~any(layerMap(:) > 0)
        continue;
    end
    layerCenter = mean([info.zBottom_nm, info.zTop_nm]);
    layerHalfWidth = max(0.5 * (info.zTop_nm - info.zBottom_nm), 0.5);
    zIndices = find(zMask);
    for localIdx = 1:numel(zIndices)
        zIdx = zIndices(localIdx);
        depthEnvelope = exp(-0.5 * ((zGrid(zIdx) - layerCenter) / max(layerHalfWidth, 0.7)) .^ 2);
        sliceMap = layerMap * max(depthEnvelope, 0.25);
        sliceMap = applyVpitActivityModulation(sliceMap, xGrid, yGrid, zGrid(zIdx), vPitData, info.channel, maps);
        positiveMask = isfinite(sliceMap) & sliceMap > 0;
        volume(:, :, zIdx) = volume(:, :, zIdx) + sliceMap;
        wavelengthNumerator(:, :, zIdx) = wavelengthNumerator(:, :, zIdx) + sliceMap * info.nominalWavelength_nm;
        supportMask(:, :, zIdx) = supportMask(:, :, zIdx) | positiveMask;
    end
end
volume(~supportMask) = NaN;
if hasWavelengthData
    wavelengthVolume = wavelengthNumerator ./ max(volume, eps);
    wavelengthVolume(~supportMask) = NaN;
else
    wavelengthVolume = [];
end
end

function layerMap = resolveActivityMapForChannel(maps, channel)
switch channel
    case 'mqw'
        layerMap = maps.mqw;
    case 'shoulder'
        layerMap = maps.shoulder;
    case 'short'
        layerMap = maps.short;
    otherwise
        layerMap = maps.all;
end
end

function sliceMap = applyVpitActivityModulation(sliceMap, xGrid, yGrid, zNm, vPitData, channel, maps)
if isempty(vPitData) || ~isstruct(vPitData) || ~isfield(vPitData, 'Vpits') || isempty(vPitData.Vpits)
    return;
end
[X2, Y2] = meshgrid(xGrid, yGrid);
for vpitIdx = 1:numel(vPitData.Vpits)
    vpit = vPitData.Vpits{vpitIdx};
    if ~isstruct(vpit) || ~isfield(vpit, 'surfaceZ_nm') || ~isfield(vpit, 'depth') || vpit.depth <= 0
        continue;
    end
    depthBelowSurface = vpit.surfaceZ_nm - zNm;
    if depthBelowSurface < 0 || depthBelowSurface > vpit.depth
        continue;
    end
    apothem = calculateDisplayApothem(vpit, zNm);
    if ~isfinite(apothem) || apothem <= 0.5
        continue;
    end
    orientationDeg = 0;
    if isfield(vpit, 'orientationDeg') && isfinite(vpit.orientationDeg)
        orientationDeg = vpit.orientationDeg;
    end
    [outerX, outerY] = buildHexagonFromApothem(vpit.center, 1.06 * apothem, orientationDeg);
    [innerX, innerY] = buildHexagonFromApothem(vpit.center, 0.70 * apothem, orientationDeg);
    outerMask = inpolygon(X2, Y2, outerX, outerY);
    innerMask = inpolygon(X2, Y2, innerX, innerY);
    shellMask = outerMask & ~innerMask;
    sliceMap(innerMask) = 0;
    switch channel
        case 'mqw'
            shortContribution = 0.26 * maps.short;
            sliceMap(shellMask) = max(sliceMap(shellMask), shortContribution(shellMask));
        case 'shoulder'
            sliceMap(shellMask) = 1.20 * sliceMap(shellMask);
        case 'short'
            sliceMap(shellMask) = 1.45 * sliceMap(shellMask);
    end
end
end

function vPitData = buildVpitDataFromPackage(activityPackage, layers, layerBoundaries)
vPitData = [];
if ~isstruct(activityPackage) || ~isfield(activityPackage, 'vPits') || isempty(activityPackage.vPits) || ...
        ~istable(activityPackage.vPits) || height(activityPackage.vPits) < 1 || isempty(layers) || isempty(layerBoundaries)
    return;
end
vpitTable = activityPackage.vPits;
names = lower(string(vpitTable.Properties.VariableNames));
centerXCol = findColumn(names, ["centerx_nm", "centerx"]);
centerYCol = findColumn(names, ["centery_nm", "centery"]);
radiusCol = findColumn(names, ["topradius_nm", "radius"]);
depthCol = findColumn(names, ["depth_nm", "depth"]);
orientationCol = findColumn(names, ["orientation_deg", "orientation"]);
if any([centerXCol, centerYCol, radiusCol] == 0)
    return;
end
centerX = double(vpitTable{:, centerXCol});
centerY = double(vpitTable{:, centerYCol});
radii = double(vpitTable{:, radiusCol});
if depthCol > 0
    depths = double(vpitTable{:, depthCol});
else
    depths = nan(size(radii));
end
if orientationCol > 0
    orientations = double(vpitTable{:, orientationCol});
else
    orientations = zeros(size(radii));
end
count = min([numel(centerX), numel(centerY), numel(radii)]);
rawVpits = cell(count, 1);
for idx = 1:count
    if ~isfinite(depths(idx)) || depths(idx) <= 0
        depths(idx) = precies.vpit('calculateDepthFromComposition', radii(idx), layers, layerBoundaries);
    end
    if ~isfinite(orientations(idx))
        orientations(idx) = 0;
    end
    rawVpits{idx} = struct( ...
        'center', [centerX(idx), centerY(idx)], ...
        'topRadius', radii(idx), ...
        'depth', depths(idx), ...
        'orientationDeg', orientations(idx));
end
vPitData = struct('Vpits', {precies.vpit('buildGeometries', rawVpits, layers, layerBoundaries)});
end

function edges = buildGridEdges(gridAxis)
step = median(diff(gridAxis));
edges = [gridAxis(1) - 0.5 * step, 0.5 * (gridAxis(1:end-1) + gridAxis(2:end)), gridAxis(end) + 0.5 * step];
end

function kernel = buildGaussianKernel2D(sigmaY, sigmaX)
radiusY = max(2, ceil(3 * sigmaY));
radiusX = max(2, ceil(3 * sigmaX));
[gridY, gridX] = ndgrid(-radiusY:radiusY, -radiusX:radiusX);
kernel = exp(-0.5 * ((gridY / sigmaY) .^ 2 + (gridX / sigmaX) .^ 2));
kernel = kernel / max(sum(kernel(:)), eps);
end

function apothem = calculateDisplayApothem(vpit, zNm)
depthBelowSurface = vpit.surfaceZ_nm - zNm;
apothem = vpit.topRadius * max(0, 1 - depthBelowSurface / max(vpit.depth, eps));
end

function [xHex, yHex] = buildHexagonFromApothem(center, apothem, orientationDeg)
vertexRadius = apothem / cosd(30);
angles = orientationDeg + (30:60:390);
xHex = center(1) + vertexRadius * cosd(angles);
yHex = center(2) + vertexRadius * sind(angles);
end

function [X, Y, Z, V, W, depthRange] = applyActiveStackDepthFilter(activityPackage, X, Y, Z, V, W)
depthRange = resolveActivityDepthRange(activityPackage, Z);
inDepthRange = Z >= depthRange.zMin_nm & Z <= depthRange.zMax_nm;
minimumPoints = max(8, ceil(0.002 * numel(Z)));
if nnz(inDepthRange) >= minimumPoints
    X = X(inDepthRange);
    Y = Y(inDepthRange);
    Z = Z(inDepthRange);
    V = V(inDepthRange);
    if ~isempty(W)
        W = W(inDepthRange);
    end
    depthRange.applied = true;
    depthRange.pointCount = numel(V);
else
    depthRange = fallbackActivityDepthRange(min(Z, [], 'omitnan'), max(Z, [], 'omitnan'), ...
        sprintf('%s; fallback to point-cloud range because too few points were inside the active-stack window', depthRange.strategy));
    depthRange.applied = false;
    depthRange.pointCount = numel(V);
end
end

function depthRange = resolveActivityDepthRange(activityPackage, Z)
finiteZ = Z(isfinite(Z));
if isempty(finiteZ)
    depthRange = fallbackActivityDepthRange(0, 1, 'empty point-cloud fallback');
    return;
end
dataZMin = min(finiteZ);
dataZMax = max(finiteZ);
marginNm = 10;
depthRange = fallbackActivityDepthRange(dataZMin, dataZMax, 'point-cloud range');
depthRange.margin_nm = marginNm;
if ~isstruct(activityPackage) || ~isfield(activityPackage, 'layers') || ~istable(activityPackage.layers)
    return;
end
layers = activityPackage.layers;
requiredNames = {'LayerName', 'ZBottom_nm', 'ZTop_nm'};
if ~all(ismember(requiredNames, layers.Properties.VariableNames))
    return;
end
[activeBottom, activeTop, activeNames] = resolveActivityStackBounds(layers);
finiteTop = double(layers.ZTop_nm(isfinite(layers.ZTop_nm)));
if isempty(finiteTop) || ~isfinite(activeBottom)
    return;
end
surfaceZ = max(finiteTop);
targetDepth = surfaceZ - activeBottom + marginNm;
if ~isfinite(targetDepth) || targetDepth <= 0
    return;
end
zMin = surfaceZ - targetDepth;
zMax = surfaceZ;
if zMin >= zMax
    return;
end
activeToken = strjoin(cellstr(activeNames), ', ');
depthRange = struct( ...
    'zMin_nm', zMin, ...
    'zMax_nm', zMax, ...
    'surfaceZ_nm', surfaceZ, ...
    'depthLimit_nm', targetDepth, ...
    'activeBottomZ_nm', activeBottom, ...
    'activeTopZ_nm', activeTop, ...
    'margin_nm', marginNm, ...
    'applied', false, ...
    'pointCount', 0, ...
    'strategy', sprintf('active stack (%s) bottom + %.1f nm', activeToken, marginNm));
end

function depthRange = fallbackActivityDepthRange(zMin, zMax, strategy)
if ~isfinite(zMin)
    zMin = 0;
end
if ~isfinite(zMax)
    zMax = zMin + 1;
end
if zMax <= zMin
    zMax = zMin + 1;
end
depthRange = struct( ...
    'zMin_nm', zMin, ...
    'zMax_nm', zMax, ...
    'surfaceZ_nm', zMax, ...
    'depthLimit_nm', zMax - zMin, ...
    'activeBottomZ_nm', NaN, ...
    'activeTopZ_nm', NaN, ...
    'margin_nm', NaN, ...
    'applied', false, ...
    'pointCount', 0, ...
    'strategy', strategy);
end

function [activeBottom, activeTop, activeNames] = resolveActivityStackBounds(layers)
activeBottom = NaN;
activeTop = NaN;
activeNames = strings(0, 1);
layerNames = string(layers.LayerName);
activeMask = contains(layerNames, "MQW-Well", 'IgnoreCase', true) | ...
    contains(layerNames, "MQW-Barrier", 'IgnoreCase', true);
if ~any(activeMask)
    activeMask = contains(layerNames, "MQW", 'IgnoreCase', true);
end
activeIdx = find(activeMask);
if isempty(activeIdx)
    return;
end
bottomCandidates = double(layers.ZBottom_nm(activeIdx));
topCandidates = double(layers.ZTop_nm(activeIdx));
finiteMask = isfinite(bottomCandidates) & isfinite(topCandidates);
bottomCandidates = bottomCandidates(finiteMask);
topCandidates = topCandidates(finiteMask);
activeIdx = activeIdx(finiteMask);
if isempty(bottomCandidates)
    return;
end
activeBottom = min(bottomCandidates);
activeTop = max(topCandidates);
activeNames = unique(layerNames(activeIdx), 'stable');
end

function maps = resolveBaseMaps(activityPackage, X, Y, V, W, xGrid, yGrid)
maps = struct();
if isfield(activityPackage, 'activityMaps') && isstruct(activityPackage.activityMaps)
    maps = activityPackage.activityMaps;
end
if ~isfield(maps, 'all') || isempty(maps.all)
    maps.all = pointMap(X, Y, V, xGrid, yGrid, true(size(V)));
end
if ~isfield(maps, 'mqw') || isempty(maps.mqw)
    maps.mqw = pointMap(X, Y, V, xGrid, yGrid, W >= 470 & W <= 570);
end
if ~isfield(maps, 'shoulder') || isempty(maps.shoulder)
    maps.shoulder = pointMap(X, Y, V, xGrid, yGrid, W >= 405 & W < 470);
end
if ~isfield(maps, 'short') || isempty(maps.short)
    maps.short = pointMap(X, Y, V, xGrid, yGrid, W >= 370 & W < 405);
end
maps.all = resizeMap(maps.all, numel(yGrid), numel(xGrid));
maps.mqw = resizeMap(max(maps.mqw, 0.65 * maps.all), numel(yGrid), numel(xGrid));
maps.shoulder = resizeMap(max(maps.shoulder, 0.22 * maps.all), numel(yGrid), numel(xGrid));
maps.short = resizeMap(max(maps.short, 0.20 * maps.all), numel(yGrid), numel(xGrid));
end

function weight = getLayerWeight(layers, idx)
if ismember('ActivityWeight', layers.Properties.VariableNames)
    weight = double(layers.ActivityWeight(idx));
else
    name = char(string(layers.LayerName(idx)));
    if contains(name, 'MQW-Well')
        weight = 1.0;
    elseif contains(name, 'MQW-Barrier')
        weight = 0.18;
    elseif contains(name, 'Prestrained')
        weight = 0.09;
    elseif contains(name, 'GaN') || contains(name, 'EBL')
        weight = 0.06;
    else
        weight = 0;
    end
end
end

function wavelength = getLayerNominalWavelength(layers, idx, layerName)
if ismember('NominalWavelength_nm', layers.Properties.VariableNames)
    wavelength = double(layers.NominalWavelength_nm(idx));
else
    wavelength = NaN;
end
if isfinite(wavelength) && wavelength > 0
    return;
end
text = char(layerName);
if contains(text, 'MQW-Well')
    wavelength = 515;
elseif contains(text, 'MQW-Barrier')
    wavelength = 445;
elseif contains(text, 'Prestrained')
    wavelength = 420;
elseif contains(text, 'GaN') || contains(text, 'EBL')
    wavelength = 390;
else
    wavelength = NaN;
end
end

function map = selectLayerMap(maps, layerName)
text = char(layerName);
if contains(text, 'MQW-Well')
    map = maps.mqw;
elseif contains(text, 'MQW-Barrier') || contains(text, 'Prestrained')
    map = maps.shoulder;
elseif contains(text, 'GaN') || contains(text, 'EBL')
    map = maps.short;
else
    map = zeros(size(maps.all));
end
end

function map = pointMap(X, Y, V, xGrid, yGrid, mask)
map = zeros(numel(yGrid), numel(xGrid));
mask = mask(:) & isfinite(X(:)) & isfinite(Y(:)) & isfinite(V(:)) & V(:) > 0;
if nnz(mask) < 1
    return;
end
xEdges = gridEdges(xGrid);
yEdges = gridEdges(yGrid);
ix = discretize(X(mask), xEdges);
iy = discretize(Y(mask), yEdges);
values = V(mask);
valid = isfinite(ix) & isfinite(iy) & isfinite(values);
if any(valid)
    map = accumarray([iy(valid), ix(valid)], values(valid), size(map), @sum, 0);
    map = conv2(map, gaussianKernel(2), 'same');
end
map = normalizeMap(map);
end

function map = resizeMap(map, rows, cols)
if size(map, 1) == rows && size(map, 2) == cols
    map = normalizeMap(double(map));
else
    map = normalizeMap(imresize(double(map), [rows, cols]));
end
end

function wavelengthVolume = buildPointWavelengthVolume(X, Y, Z, W, Xg, Yg, Zg, intensityVolume)
wavelengthVolume = nan(size(intensityVolume));
valid = isfinite(X) & isfinite(Y) & isfinite(Z) & isfinite(W) & W > 0;
if nnz(valid) < 4
    return;
end
try
    wavelengthInterpolant = scatteredInterpolant(X(valid), Y(valid), Z(valid), W(valid), 'natural', 'linear');
    wavelengthVolume = wavelengthInterpolant(Xg, Yg, Zg);
    wavelengthVolume(~isfinite(wavelengthVolume) | intensityVolume <= 0) = NaN;
catch
    wavelengthVolume(:) = NaN;
end
end

function logVolume = buildLogVolume(volume)
values = volume(isfinite(volume) & volume > 0);
logVolume = nan(size(volume));
if isempty(values)
    return;
end
floorValue = max(min(values), eps);
logVolume = log10(max(volume, floorValue));
logVolume(~isfinite(volume)) = NaN;
end

function volume = normalizeVolume(volume)
values = volume(isfinite(volume) & volume > 0);
if isempty(values)
    volume(:) = 0;
    return;
end
scale = prctile(values, 98);
if ~isfinite(scale) || scale <= 0
    scale = max(values);
end
volume = min(max(volume ./ max(scale, eps), 0), 1.5);
end

function value = getPackageTable(package, name)
if isfield(package, name)
    value = package.(name);
else
    value = table();
end
end

function showLegacySliceViewer(volumeData)
fig = uifigure('Name', 'Legacy Grid Slice Model', 'Position', [120, 120, 920, 720]);
layout = uigridlayout(fig, [2, 1]);
layout.RowHeight = {'1x', 'fit'};
ax = uiaxes(layout);
ax.Layout.Row = 1;
controls = uigridlayout(layout, [4, 8]);
controls.Layout.Row = 2;
controls.ColumnWidth = {'fit', 120, 'fit', 90, 'fit', '1x', 90, 150};
controls.RowHeight = {'fit', 'fit', 'fit', 'fit'};
modeLabel = uilabel(controls, 'Text', 'Mode');
modeLabel.Layout.Row = 1;
modeLabel.Layout.Column = 1;
modeDrop = uidropdown(controls, 'Items', {'Intensity', 'Log Intensity', 'Wavelength'}, 'Value', 'Intensity');
modeDrop.Layout.Row = 1;
modeDrop.Layout.Column = 2;
directionLabel = uilabel(controls, 'Text', 'Direction');
directionLabel.Layout.Row = 1;
directionLabel.Layout.Column = 3;
directionDrop = uidropdown(controls, 'Items', {'X', 'Y', 'Z'}, 'Value', 'Z');
directionDrop.Layout.Row = 1;
directionDrop.Layout.Column = 4;
positionLabel = uilabel(controls, 'Text', 'Position (nm)');
positionLabel.Layout.Row = 1;
positionLabel.Layout.Column = 5;
positionSlider = uislider(controls);
positionSlider.Layout.Row = 1;
positionSlider.Layout.Column = 6;
positionEdit = uieditfield(controls, 'numeric');
positionEdit.Layout.Row = 1;
positionEdit.Layout.Column = 7;
colormapDrop = uidropdown(controls, 'Items', {'Turbo', 'Parula', 'Hot', 'Jet', 'Gray', 'Spring', 'Cool'}, 'Value', 'Turbo');
colormapDrop.Layout.Row = 1;
colormapDrop.Layout.Column = 8;
minRangeLabel = uilabel(controls, 'Text', 'Min Range');
minRangeLabel.Layout.Row = 2;
minRangeLabel.Layout.Column = 1;
minRangeSlider = uislider(controls);
minRangeSlider.Layout.Row = 2;
minRangeSlider.Layout.Column = [2, 6];
minRangeValue = uilabel(controls, 'Text', '0');
minRangeValue.HorizontalAlignment = 'right';
minRangeValue.Layout.Row = 2;
minRangeValue.Layout.Column = 7;
maxRangeLabel = uilabel(controls, 'Text', 'Max Range');
maxRangeLabel.Layout.Row = 3;
maxRangeLabel.Layout.Column = 1;
maxRangeSlider = uislider(controls);
maxRangeSlider.Layout.Row = 3;
maxRangeSlider.Layout.Column = [2, 6];
maxRangeValue = uilabel(controls, 'Text', '1');
maxRangeValue.HorizontalAlignment = 'right';
maxRangeValue.Layout.Row = 3;
maxRangeValue.Layout.Column = 7;
autoRangeButton = uibutton(controls, 'push', 'Text', 'Auto');
autoRangeButton.Layout.Row = [2, 3];
autoRangeButton.Layout.Column = 8;
statusLabel = uilabel(controls, 'Text', volumeData.method);
statusLabel.Layout.Row = 4;
statusLabel.Layout.Column = [1, 8];
modeDrop.ValueChangedFcn = @(~, ~) resetColorRange();
directionDrop.ValueChangedFcn = @(~, ~) resetPosition();
positionSlider.ValueChangedFcn = @(~, ~) updateFromSlider();
positionEdit.ValueChangedFcn = @(~, ~) updateFromEdit();
colormapDrop.ValueChangedFcn = @(~, ~) updateSlice();
minRangeSlider.ValueChangedFcn = @(~, ~) updateMinRange();
maxRangeSlider.ValueChangedFcn = @(~, ~) updateMaxRange();
autoRangeButton.ButtonPushedFcn = @(~, ~) resetColorRange();
resetColorRange();
resetPosition();

    function resetColorRange()
        modeName = modeDrop.Value;
        [~, colorLimits] = selectSliceVolumeMode(volumeData, modeName);
        if isempty(colorLimits) || ~all(isfinite(colorLimits)) || colorLimits(2) <= colorLimits(1)
            colorLimits = [0, 1];
        end
        minRangeSlider.Limits = colorLimits;
        maxRangeSlider.Limits = colorLimits;
        minRangeSlider.Value = colorLimits(1);
        maxRangeSlider.Value = colorLimits(2);
        updateRangeLabels();
        updateSlice();
    end

    function updateMinRange()
        if minRangeSlider.Value >= maxRangeSlider.Value
            minRangeSlider.Value = maxRangeSlider.Value - max(eps, 0.001 * diff(maxRangeSlider.Limits));
        end
        updateRangeLabels();
        updateSlice();
    end

    function updateMaxRange()
        if maxRangeSlider.Value <= minRangeSlider.Value
            maxRangeSlider.Value = minRangeSlider.Value + max(eps, 0.001 * diff(minRangeSlider.Limits));
        end
        updateRangeLabels();
        updateSlice();
    end

    function updateRangeLabels()
        minRangeValue.Text = formatRangeValue(minRangeSlider.Value, modeDrop.Value);
        maxRangeValue.Text = formatRangeValue(maxRangeSlider.Value, modeDrop.Value);
    end

    function resetPosition()
        switch directionDrop.Value
            case 'X'
                limits = [min(volumeData.xGrid), max(volumeData.xGrid)];
            case 'Y'
                limits = [min(volumeData.yGrid), max(volumeData.yGrid)];
            otherwise
                depthGrid = buildDisplayDepthGrid(volumeData.zGrid);
                limits = [min(depthGrid), max(depthGrid)];
        end
        if limits(1) == limits(2)
            limits = limits + [-1, 1];
        end
        positionSlider.Limits = limits;
        positionSlider.Value = mean(limits);
        positionEdit.Value = positionSlider.Value;
        updateSlice();
    end

    function updateFromSlider()
        positionEdit.Value = positionSlider.Value;
        updateSlice();
    end

    function updateFromEdit()
        positionEdit.Value = min(max(positionEdit.Value, positionSlider.Limits(1)), positionSlider.Limits(2));
        positionSlider.Value = positionEdit.Value;
        updateSlice();
    end

    function updateSlice()
        direction = directionDrop.Value;
        position = positionSlider.Value;
        modeName = modeDrop.Value;
        [V, colorLimits, colorLabel] = selectSliceVolumeMode(volumeData, modeName);
        depthGrid = buildDisplayDepthGrid(volumeData.zGrid);
        V = flip(V, 3);
        switch direction
            case 'X'
                [~, idx] = min(abs(volumeData.xGrid - position));
                imageData = squeeze(V(:, idx, :))';
                xData = volumeData.yGrid;
                yData = depthGrid;
                xlabelText = 'Y (nm)';
                ylabelText = 'Depth (nm)';
                titleText = sprintf('X = %.2f nm', volumeData.xGrid(idx));
            case 'Y'
                [~, idx] = min(abs(volumeData.yGrid - position));
                imageData = squeeze(V(idx, :, :))';
                xData = volumeData.xGrid;
                yData = depthGrid;
                xlabelText = 'X (nm)';
                ylabelText = 'Depth (nm)';
                titleText = sprintf('Y = %.2f nm', volumeData.yGrid(idx));
            otherwise
                [~, idx] = min(abs(depthGrid - position));
                imageData = V(:, :, idx);
                xData = volumeData.xGrid;
                yData = volumeData.yGrid;
                xlabelText = 'X (nm)';
                ylabelText = 'Y (nm)';
                titleText = sprintf('Depth = %.2f nm', depthGrid(idx));
        end
        rangeLimits = sort([minRangeSlider.Value, maxRangeSlider.Value]);
        imageData(imageData < rangeLimits(1) | imageData > rangeLimits(2)) = NaN;
        imageHandle = imagesc(ax, xData, yData, imageData);
        set(imageHandle, 'AlphaData', double(isfinite(imageData)));
        ax.Color = [1, 1, 1];
        axis(ax, 'tight');
        if strcmp(direction, 'Z')
            axis(ax, 'image');
            set(ax, 'YDir', 'normal');
        else
            axis(ax, 'normal');
            pbaspect(ax, [3, 1.6, 1]);
            set(ax, 'YDir', 'reverse');
        end
        xlim(ax, [min(xData), max(xData)]);
        ylim(ax, [min(yData), max(yData)]);
        applyEndpointTicks(ax, 'x', xData);
        applyEndpointTicks(ax, 'y', yData);
        colormap(ax, resolveColormap(colormapDrop.Value, 256));
        cb = colorbar(ax);
        cb.Label.String = colorLabel;
        if all(isfinite(rangeLimits)) && rangeLimits(2) > rangeLimits(1)
            clim(ax, rangeLimits);
        elseif ~isempty(colorLimits)
            clim(ax, colorLimits);
        end
        xlabel(ax, xlabelText);
        ylabel(ax, ylabelText);
        title(ax, ['Legacy Grid Slice, ' titleText]);
        statusLabel.Text = sprintf('Legacy Grid slice: %s, %s direction, position %.2f nm', modeName, direction, position);
    end
end

function showActivityStructureVolume(volumeData)
fig = uifigure('Name', 'Activity-Structure 3D Model', 'Position', [180, 90, 980, 760]);
layout = uigridlayout(fig, [2, 1]);
layout.RowHeight = {'1x', 'fit'};
ax = uiaxes(layout);
ax.Layout.Row = 1;
controls = uigridlayout(layout, [2, 7]);
controls.Layout.Row = 2;
controls.ColumnWidth = {'fit', 140, 'fit', 150, 'fit', '1x', 120};
controls.RowHeight = {'fit', 'fit'};
modeLabel = uilabel(controls, 'Text', 'Mode');
modeLabel.Layout.Row = 1;
modeLabel.Layout.Column = 1;
modeDrop = uidropdown(controls, 'Items', {'Intensity', 'Log Intensity', 'Wavelength'}, 'Value', 'Intensity');
modeDrop.Layout.Row = 1;
modeDrop.Layout.Column = 2;
colormapLabel = uilabel(controls, 'Text', 'Color Map');
colormapLabel.Layout.Row = 1;
colormapLabel.Layout.Column = 3;
colormapDrop = uidropdown(controls, 'Items', {'Turbo', 'Parula', 'Hot', 'Jet', 'Gray', 'Spring', 'Cool'}, 'Value', 'Turbo');
colormapDrop.Layout.Row = 1;
colormapDrop.Layout.Column = 4;
alphaTitle = uilabel(controls, 'Text', 'Alpha');
alphaTitle.Layout.Row = 1;
alphaTitle.Layout.Column = 5;
alphaSlider = uislider(controls, 'Limits', [0.05, 0.85], 'Value', 0.32);
alphaSlider.Layout.Row = 1;
alphaSlider.Layout.Column = 6;
alphaLabel = uilabel(controls, 'Text', '0.32');
alphaLabel.Layout.Row = 1;
alphaLabel.Layout.Column = 7;
statusLabel = uilabel(controls, 'Text', volumeData.method);
statusLabel.Layout.Row = 2;
statusLabel.Layout.Column = [1, 7];
modeDrop.ValueChangedFcn = @(~, ~) renderVolume();
colormapDrop.ValueChangedFcn = @(~, ~) renderVolume();
alphaSlider.ValueChangedFcn = @(~, ~) renderVolume();
renderVolume();

    function renderVolume()
        cla(ax);
        alphaValue = alphaSlider.Value;
        alphaLabel.Text = sprintf('%.2f', alphaValue);
        modeName = modeDrop.Value;
        [surfaceVolume, colorVolume, colorLimits, colorLabel] = selectActivityVolumeMode(volumeData, modeName);
        displayZGrid = buildDisplayDepthGrid(volumeData.zGrid);
        geometryVolume = flip(surfaceVolume, 3);
        if ~isempty(colorVolume)
            colorVolume = flip(colorVolume, 3);
        end
        finiteGeometry = geometryVolume(isfinite(geometryVolume));
        if strcmp(char(modeName), 'Log Intensity')
            values = finiteGeometry(finiteGeometry > min(finiteGeometry));
        else
            values = finiteGeometry(finiteGeometry > 0);
        end
        if isempty(values)
            text(ax, 0.5, 0.5, 'No activity volume data', 'Units', 'normalized', 'HorizontalAlignment', 'center');
            return;
        end
        cmap = resolveColormap(colormapDrop.Value, 256);
        levels = unique(prctile(values, [55, 70, 84, 93]));
        levels = levels(levels > min(values) & levels < max(values));
        if isempty(levels)
            levels = max(values) * [0.35, 0.55, 0.75];
        end
        hold(ax, 'on');
        for idx = 1:numel(levels)
            faceColor = cmap(max(1, min(256, round(1 + 255 * idx / numel(levels)))), :);
            try
                if isempty(colorVolume)
                    patchData = patch(ax, isosurface(volumeData.xGrid, volumeData.yGrid, displayZGrid, geometryVolume, levels(idx)));
                    patchData.FaceColor = faceColor;
                else
                    patchData = patch(ax, isosurface(volumeData.xGrid, volumeData.yGrid, displayZGrid, geometryVolume, levels(idx), colorVolume));
                    patchData.FaceColor = 'interp';
                end
                isonormals(volumeData.xGrid, volumeData.yGrid, displayZGrid, geometryVolume, patchData);
                patchData.EdgeColor = 'none';
                patchData.FaceAlpha = min(0.9, alphaValue * (0.55 + 0.12 * idx));
            catch
            end
        end
        hold(ax, 'off');
        view(ax, 3);
        axis(ax, 'vis3d');
        grid(ax, 'on');
        box(ax, 'on');
        xlabel(ax, 'X (nm)');
        ylabel(ax, 'Y (nm)');
        zlabel(ax, 'Depth (nm)');
        title(ax, 'Activity-Structure 3D Volume');
        xlim(ax, [min(volumeData.xGrid), max(volumeData.xGrid)]);
        ylim(ax, [min(volumeData.yGrid), max(volumeData.yGrid)]);
        zlim(ax, [min(displayZGrid), max(displayZGrid)]);
        applyEndpointTicks(ax, 'x', volumeData.xGrid);
        applyEndpointTicks(ax, 'y', volumeData.yGrid);
        applyEndpointTicks(ax, 'z', displayZGrid);
        set(ax, 'ZDir', 'reverse');
        view(ax, 42, 22);
        camlight(ax, 'headlight');
        camlight(ax, 'right');
        lighting(ax, 'gouraud');
        colormap(ax, cmap);
        cb = colorbar(ax);
        cb.Label.String = colorLabel;
        if ~isempty(colorLimits)
            clim(ax, colorLimits);
        end
        statusLabel.Text = sprintf('Activity-Structure volume: %s, %d isosurfaces, alpha %.2f', modeName, numel(levels), alphaValue);
    end
end

function [surfaceVolume, colorVolume, colorLimits, colorLabel] = selectActivityVolumeMode(volumeData, modeName)
intensityVolume = getVolumeField(volumeData, 'intensityVolume', getVolumeField(volumeData, 'volume', []));
logVolume = getVolumeField(volumeData, 'logVolume', buildLogVolume(intensityVolume));
wavelengthVolume = getVolumeField(volumeData, 'wavelengthVolume', []);
switch char(modeName)
    case 'Log Intensity'
        surfaceVolume = logVolume;
        colorVolume = [];
        colorLimits = finiteColorLimits(logVolume);
        colorLabel = 'Log intensity';
    case 'Wavelength'
        surfaceVolume = intensityVolume;
        colorVolume = wavelengthVolume;
        validWavelength = wavelengthVolume(isfinite(wavelengthVolume) & wavelengthVolume > 0);
        if isempty(validWavelength)
            colorVolume = [];
            colorLimits = [];
            colorLabel = 'Intensity';
        else
            colorLimits = [min(validWavelength), max(validWavelength)];
            if colorLimits(1) == colorLimits(2)
                colorLimits = colorLimits + [-1, 1];
            end
            colorLabel = 'Wavelength (nm)';
        end
    otherwise
        surfaceVolume = intensityVolume;
        colorVolume = [];
        colorLimits = finiteColorLimits(intensityVolume);
        colorLabel = 'Intensity';
end
surfaceVolume = fillMissingVolumeValues(surfaceVolume, modeName);
if ~isempty(colorVolume)
    colorVolume(~isfinite(colorVolume)) = NaN;
end
end

function [displayVolume, colorLimits, colorLabel] = selectSliceVolumeMode(volumeData, modeName)
intensityVolume = getVolumeField(volumeData, 'intensityVolume', getVolumeField(volumeData, 'volume', []));
logVolume = getVolumeField(volumeData, 'logVolume', buildLogVolume(intensityVolume));
wavelengthVolume = getVolumeField(volumeData, 'wavelengthVolume', []);
switch char(modeName)
    case 'Log Intensity'
        displayVolume = logVolume;
        colorLimits = finiteColorLimits(logVolume);
        colorLabel = 'Log intensity';
    case 'Wavelength'
        validWavelength = wavelengthVolume(isfinite(wavelengthVolume) & wavelengthVolume > 0);
        if isempty(validWavelength)
            displayVolume = intensityVolume;
            colorLimits = finiteColorLimits(intensityVolume);
            colorLabel = 'Intensity';
        else
            displayVolume = wavelengthVolume;
            colorLimits = [min(validWavelength), max(validWavelength)];
            if colorLimits(1) == colorLimits(2)
                colorLimits = colorLimits + [-1, 1];
            end
            colorLabel = 'Wavelength (nm)';
        end
    otherwise
        displayVolume = intensityVolume;
        colorLimits = finiteColorLimits(intensityVolume);
        colorLabel = 'Intensity';
end
displayVolume = fillMissingVolumeValues(displayVolume, modeName);
end

function depthGrid = buildDisplayDepthGrid(zGrid)
surfaceZ = max(zGrid);
depthGrid = sort(surfaceZ - zGrid, 'ascend');
end

function values = fillMissingVolumeValues(values, modeName)
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    values = zeros(size(values));
    return;
end
if strcmp(char(modeName), 'Log Intensity')
    valueRange = max(finiteValues) - min(finiteValues);
    fillValue = min(finiteValues) - max(0.02 * valueRange, eps);
else
    fillValue = 0;
end
values(~isfinite(values)) = fillValue;
end

function applyEndpointTicks(ax, axisName, axisData)
axisData = axisData(isfinite(axisData));
if isempty(axisData)
    return;
end
limits = [min(axisData), max(axisData)];
if limits(1) == limits(2)
    return;
end
if max(abs(limits - [-250, 250])) < 1e-6
    ticks = [-250, -125, 0, 125, 250];
else
    ticks = linspace(limits(1), limits(2), 5);
end
switch char(axisName)
    case 'x'
        xticks(ax, ticks);
    case 'y'
        yticks(ax, ticks);
    case 'z'
        zticks(ax, ticks);
end
end

function textValue = formatRangeValue(value, modeName)
if strcmp(char(modeName), 'Wavelength')
    textValue = sprintf('%.1f', value);
elseif strcmp(char(modeName), 'Log Intensity')
    textValue = sprintf('%.3f', value);
else
    textValue = sprintf('%.4g', value);
end
end

function limits = finiteColorLimits(values)
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    limits = [];
    return;
end
limits = [min(finiteValues), max(finiteValues)];
if limits(1) == limits(2)
    limits = limits + [-1, 1];
end
end

function value = getVolumeField(volumeData, fieldName, defaultValue)
if isstruct(volumeData) && isfield(volumeData, fieldName) && ~isempty(volumeData.(fieldName))
    value = volumeData.(fieldName);
else
    value = defaultValue;
end
end

function cmap = resolveColormap(name, n)
if nargin < 2 || isempty(n)
    n = 256;
end
name = char(string(name));
switch lower(name)
    case 'turbo'
        try
            cmap = turbo(n);
        catch
            cmap = parula(n);
        end
    case 'parula'
        cmap = parula(n);
    case 'hot'
        cmap = hot(n);
    case 'jet'
        cmap = jet(n);
    case 'gray'
        cmap = gray(n);
    case 'spring'
        cmap = spring(n);
    case 'cool'
        cmap = cool(n);
    otherwise
        cmap = parula(n);
end
end

function fitResult = fitAdaptiveGaussianMixture(axisNm, observedSpectrum)
    axisNm = axisNm(:);
    observedSpectrum = normalizeSpectrum(observedSpectrum(:));
    maxAdaptivePeaks = 8;
    peakCandidates = detectPeakCandidates(axisNm, observedSpectrum, maxAdaptivePeaks);
    physicsCandidates = detectPhysicsWindowCandidates(axisNm, observedSpectrum);
    peakCandidates = mergePeakCandidates([peakCandidates; physicsCandidates], maxAdaptivePeaks);
    primaryFit = fitMixtureFromCandidates(axisNm, observedSpectrum, peakCandidates);

    residualCandidates = detectResidualCandidates( ...
        axisNm, observedSpectrum, primaryFit.fittedSpectrum, primaryFit.componentsNm, maxAdaptivePeaks);
    if ~isempty(residualCandidates)
        augmentedCandidates = mergePeakCandidates([peakCandidates; residualCandidates], maxAdaptivePeaks);
        augmentedFit = fitMixtureFromCandidates(axisNm, observedSpectrum, augmentedCandidates);
        useAugmentedFit = augmentedFit.fitError <= primaryFit.fitError - 0.004 || ...
            (augmentedFit.fitError <= primaryFit.fitError + 0.001 && augmentedFit.numPeaks > primaryFit.numPeaks);
        if useAugmentedFit
            fitState = augmentedFit;
        else
            fitState = primaryFit;
        end
    else
        fitState = primaryFit;
    end

    fittedSpectrum = fitState.fittedSpectrum;
    componentsNm = fitState.componentsNm;
    componentSpectra = fitState.componentSpectra;
    baseline = fitState.baseline;
    significantComponents = selectSignificantComponents(componentsNm);
    dominantIdx = find(significantComponents(:,1) == max(significantComponents(:,1)), 1, 'first');
    windowMetrics = computeCalibrationWindowMetrics(axisNm, observedSpectrum, fittedSpectrum);
    shortwaveComponent = selectComponentInBand(significantComponents, [376, 396], windowMetrics.shortwavePeakNm);
    mainComponent = selectComponentInBand(significantComponents, [492, 545], windowMetrics.mainPeakNm);

    fitResult = struct( ...
        'fittedSpectrum', fittedSpectrum, ...
        'fitError', fitState.fitError, ...
        'numPeaks', size(significantComponents, 1), ...
        'componentsNm', significantComponents, ...
        'allComponentsNm', componentsNm, ...
        'significantComponentsNm', significantComponents, ...
        'componentSpectra', componentSpectra, ...
        'baseline', baseline, ...
        'dominantComponent', significantComponents(dominantIdx, :), ...
        'shortwaveComponent', shortwaveComponent, ...
        'blueComponent', shortwaveComponent, ...
        'redComponent', mainComponent, ...
        'windowMetrics', windowMetrics, ...
        'shortwavePeakNm', windowMetrics.shortwavePeakNm, ...
        'shoulderPeakNm', windowMetrics.shoulderPeakNm, ...
        'mainPeakNm', windowMetrics.mainPeakNm, ...
        'shortwaveFitError', windowMetrics.shortwaveFitError, ...
        'shoulderFitError', windowMetrics.shoulderFitError, ...
        'shortwaveToMainAreaRatio', windowMetrics.shortwaveToMainAreaRatio, ...
        'shoulderToMainAreaRatio', windowMetrics.shoulderToMainAreaRatio, ...
        'mainPeakArea', windowMetrics.mainPeakArea);
end

function component = selectComponentInBand(componentsNm, bandRangeNm, targetCenterNm)
    bandMask = componentsNm(:, 2) >= bandRangeNm(1) & componentsNm(:, 2) <= bandRangeNm(2);
    if any(bandMask)
        candidateIdx = find(bandMask);
        if isfinite(targetCenterNm)
            [~, bestLocalIdx] = min(abs(componentsNm(candidateIdx, 2) - targetCenterNm));
        else
            [~, bestLocalIdx] = max(componentsNm(candidateIdx, 1));
        end
        component = componentsNm(candidateIdx(bestLocalIdx), :);
        return;
    end

    if isfinite(targetCenterNm)
        [~, bestIdx] = min(abs(componentsNm(:, 2) - targetCenterNm));
    else
        [~, bestIdx] = max(componentsNm(:, 1));
    end
    component = componentsNm(bestIdx, :);
end

function fitState = fitMixtureFromCandidates(axisNm, observedSpectrum, peakCandidates)
    peakCount = size(peakCandidates, 1);
    initialGuess = buildInitialMixtureGuess(peakCandidates);
    fitConstraints = buildCandidateFitConstraints(peakCandidates);
    options = optimset('Display', 'off', 'MaxIter', 6000, 'MaxFunEvals', 12000, 'TolX', 1e-7, 'TolFun', 1e-7);

    bestParams = initialGuess;
    bestError = inf;
    for restartIdx = 1:4
        restartGuess = perturbMixtureGuess(initialGuess, peakCandidates, restartIdx);
        objective = @(p) fitMixtureObjective(p, peakCount, axisNm, observedSpectrum, fitConstraints);
        fittedParams = fminsearch(objective, restartGuess, options);
        currentError = fitMixtureObjective(fittedParams, peakCount, axisNm, observedSpectrum, fitConstraints);
        if currentError < bestError
            bestError = currentError;
            bestParams = fittedParams;
        end
    end

    [fittedSpectrum, componentsNm, componentSpectra, baseline] = evaluateGaussianMixtureModel(bestParams, peakCount, axisNm);
    fittedSpectrum = normalizeSpectrum(fittedSpectrum);
    fitState = struct( ...
        'fittedSpectrum', fittedSpectrum, ...
        'fitError', computeShapeError(observedSpectrum, fittedSpectrum), ...
        'componentsNm', componentsNm, ...
        'componentSpectra', componentSpectra, ...
        'baseline', baseline, ...
        'numPeaks', size(selectSignificantComponents(componentsNm), 1), ...
        'significantComponentsNm', selectSignificantComponents(componentsNm));
end

function peakCandidates = detectPeakCandidates(axisNm, observedSpectrum, maxPeaks)
    smoothedSpectrum = normalizeSpectrum(movmean(observedSpectrum(:), 5));
    localMaxMask = false(size(smoothedSpectrum));
    localMaxMask(2:end-1) = smoothedSpectrum(2:end-1) >= smoothedSpectrum(1:end-2) & ...
        smoothedSpectrum(2:end-1) > smoothedSpectrum(3:end);

    candidateIdx = find(localMaxMask);
    if isempty(candidateIdx)
        [~, candidateIdx] = max(smoothedSpectrum);
    end

    axisStepNm = median(diff(axisNm));
    if ~isfinite(axisStepNm) || axisStepNm <= 0
        axisStepNm = 1;
    end
    spacingSamples = max(2, round(6 / axisStepNm));
    accepted = [];

    for idx = candidateIdx(:)'
        peakHeight = smoothedSpectrum(idx);
        if peakHeight < max(0.018, 0.035 * max(smoothedSpectrum))
            continue;
        end
        leftMin = min(smoothedSpectrum(max(1, idx-25):idx));
        rightMin = min(smoothedSpectrum(idx:min(numel(smoothedSpectrum), idx+25)));
        prominence = peakHeight - max(leftMin, rightMin);
        if prominence < 0.0055
            continue;
        end

        if isempty(accepted) || all(abs(idx - accepted(:,1)) > spacingSamples)
            peakWidthNm = estimatePeakWidthNm(axisNm, smoothedSpectrum, idx);
            accepted = [accepted; idx, peakHeight, prominence, peakWidthNm];
        else
            [nearestDistance, nearestIdx] = min(abs(idx - accepted(:,1)));
            if nearestDistance <= spacingSamples && peakHeight > accepted(nearestIdx, 2)
                peakWidthNm = estimatePeakWidthNm(axisNm, smoothedSpectrum, idx);
                accepted(nearestIdx, :) = [idx, peakHeight, prominence, peakWidthNm];
            end
        end
    end

    accepted = addShoulderCandidates(axisNm, smoothedSpectrum, accepted, spacingSamples, maxPeaks);

    if isempty(accepted)
        [~, maxIdx] = max(smoothedSpectrum);
        accepted = [maxIdx, smoothedSpectrum(maxIdx), smoothedSpectrum(maxIdx), estimatePeakWidthNm(axisNm, smoothedSpectrum, maxIdx)];
    end

    [~, sortOrder] = sort(accepted(:,2) + 0.8 * accepted(:,3), 'descend');
    accepted = accepted(sortOrder(1:min(maxPeaks, size(accepted, 1))), :);
    [~, positionOrder] = sort(accepted(:,1), 'ascend');
    accepted = accepted(positionOrder, :);

    peakCandidates = [ ...
        normalizePositive(accepted(:,2)), ...
        axisNm(accepted(:,1)), ...
        accepted(:,4)];
end

function peakCandidates = detectPhysicsWindowCandidates(axisNm, observedSpectrum)
    axisNm = axisNm(:);
    smoothedSpectrum = normalizeSpectrum(movmean(observedSpectrum(:), 3));
    peakCandidates = zeros(0, 3);
    windows = [ ...
        376, 396, 386, 10, 0.010, 0.0025; ...
        398, 432, 414, 18, 0.012, 0.0030; ...
        442, 475, 459, 20, 0.012, 0.0030; ...
        490, 545, 517, 30, 0.040, 0.0040];

    for windowIdx = 1:size(windows, 1)
        windowRange = windows(windowIdx, 1:2);
        defaultCenterNm = windows(windowIdx, 3);
        defaultWidthNm = windows(windowIdx, 4);
        minHeight = windows(windowIdx, 5);
        minProminence = windows(windowIdx, 6);

        windowMask = axisNm >= windowRange(1) & axisNm <= windowRange(2);
        if ~any(windowMask)
            continue;
        end

        windowAxis = axisNm(windowMask);
        windowValues = smoothedSpectrum(windowMask);
        [peakHeight, localPeakIdx] = max(windowValues);
        centerNm = windowAxis(localPeakIdx);
        expandedMask = axisNm >= windowRange(1) - 18 & axisNm <= windowRange(2) + 18;
        localBaseline = lowerQuantile(smoothedSpectrum(expandedMask), 0.22);
        prominence = max(peakHeight - localBaseline, 0);

        if windowIdx < size(windows, 1)
            accepted = peakHeight >= minHeight && prominence >= minProminence;
        else
            accepted = peakHeight >= minHeight || prominence >= minProminence;
        end
        if ~accepted
            continue;
        end

        widthNm = estimateWindowPeakWidthNm(axisNm, smoothedSpectrum, centerNm, windowRange, defaultWidthNm);
        strength = max([prominence, 0.25 * peakHeight, 0.002]);
        if abs(centerNm - defaultCenterNm) > 0.55 * diff(windowRange)
            centerNm = 0.70 * centerNm + 0.30 * defaultCenterNm;
        end
        peakCandidates(end + 1, :) = [strength, centerNm, widthNm];
    end

    if ~isempty(peakCandidates)
        peakCandidates(:, 1) = normalizePositive(peakCandidates(:, 1));
    end
end

function widthNm = estimateWindowPeakWidthNm(axisNm, spectrum, centerNm, windowRange, defaultWidthNm)
    [~, peakIdx] = min(abs(axisNm - centerNm));
    estimatedWidthNm = estimatePeakWidthNm(axisNm, spectrum, peakIdx);
    windowWidthNm = diff(windowRange);
    widthNm = clampValue(0.65 * estimatedWidthNm + 0.35 * defaultWidthNm, ...
        max(8, 0.25 * defaultWidthNm), max(windowWidthNm, 1.5 * defaultWidthNm));
end

function residualCandidates = detectResidualCandidates(axisNm, observedSpectrum, fittedSpectrum, componentsNm, maxPeaks)
    residualSpectrum = normalizeSpectrum(max(observedSpectrum(:) - fittedSpectrum(:), 0));
    residualCandidates = zeros(0, 3);
    if max(residualSpectrum) < 0.035
        return;
    end

    candidateList = detectPeakCandidates(axisNm, residualSpectrum, maxPeaks);
    if isempty(candidateList)
        return;
    end

    existingCenters = componentsNm(:, 2);
    for candidateIdx = 1:size(candidateList, 1)
        candidate = candidateList(candidateIdx, :);
        if any(abs(candidate(2) - existingCenters) < max(10, 0.35 * candidate(3)))
            continue;
        end
        if candidate(2) > min(existingCenters) + 5 && candidate(2) < max(existingCenters) - 5
            residualCandidates(end + 1, :) = candidate;
        elseif candidate(1) >= 0.05
            residualCandidates(end + 1, :) = candidate;
        end
    end
end

function mergedCandidates = mergePeakCandidates(peakCandidates, maxPeaks)
    if isempty(peakCandidates)
        mergedCandidates = [1, 500, 40];
        return;
    end

    peakCandidates = double(peakCandidates);
    [~, order] = sort(peakCandidates(:, 1), 'descend');
    peakCandidates = peakCandidates(order, :);
    mergedCandidates = zeros(0, 3);
    for idx = 1:size(peakCandidates, 1)
        candidate = peakCandidates(idx, :);
        if isempty(mergedCandidates)
            mergedCandidates = candidate;
            continue;
        end
        if any(abs(candidate(2) - mergedCandidates(:, 2)) < max(8, 0.30 * candidate(3)))
            continue;
        end
        mergedCandidates(end + 1, :) = candidate;
        if size(mergedCandidates, 1) >= maxPeaks
            break;
        end
    end

    mergedCandidates(:, 1) = normalizePositive(mergedCandidates(:, 1));
    [~, positionOrder] = sort(mergedCandidates(:, 2), 'ascend');
    mergedCandidates = mergedCandidates(positionOrder, :);
end

function accepted = addShoulderCandidates(axisNm, smoothedSpectrum, accepted, spacingSamples, maxPeaks)
    if size(accepted, 1) >= maxPeaks
        return;
    end

    axisStepNm = median(diff(axisNm));
    if ~isfinite(axisStepNm) || axisStepNm <= 0
        axisStepNm = 1;
    end
    broadWindow = max(7, round(18 / axisStepNm));
    broadBaseline = movmean(smoothedSpectrum, broadWindow);
    shoulderResidual = max(0, smoothedSpectrum - broadBaseline);
    candidateMask = false(size(shoulderResidual));
    candidateMask(2:end-1) = shoulderResidual(2:end-1) >= shoulderResidual(1:end-2) & ...
        shoulderResidual(2:end-1) > shoulderResidual(3:end);
    candidateIdx = find(candidateMask);

    for idx = candidateIdx(:)'
        if size(accepted, 1) >= maxPeaks
            break;
        end
        residualHeight = shoulderResidual(idx);
        peakHeight = smoothedSpectrum(idx);
        if peakHeight < 0.02 || residualHeight < 0.004
            continue;
        end
        if ~isempty(accepted) && any(abs(idx - accepted(:, 1)) <= max(1, round(spacingSamples / 3)))
            continue;
        end
        accepted(end + 1, :) = [idx, peakHeight, residualHeight, estimatePeakWidthNm(axisNm, smoothedSpectrum, idx)];
    end
end

function initialGuess = buildInitialMixtureGuess(peakCandidates)
    weights = normalizePositive(peakCandidates(:,1));
    centersNm = peakCandidates(:,2);
    widthsNm = max(peakCandidates(:,3), 10);
    etas = 0.72 * ones(size(weights));
    baseline = 0.002;
    initialGuess = [weights(:); centersNm(:); widthsNm(:); etas(:); baseline];
end

function perturbedGuess = perturbMixtureGuess(initialGuess, peakCandidates, restartIdx)
    peakCount = size(peakCandidates, 1);
    perturbedGuess = initialGuess;
    if restartIdx == 1
        return;
    end

    rng(100 + restartIdx);
    weights = initialGuess(1:peakCount);
    centers = initialGuess(peakCount+1:2*peakCount);
    widths = initialGuess(2*peakCount+1:3*peakCount);
    etas = initialGuess(3*peakCount+1:4*peakCount);
    baseline = initialGuess(end);

    weights = max(weights .* (1 + 0.20 * randn(size(weights))), 0.05);
    centers = centers + 4.0 * randn(size(centers));
    widths = max(widths .* (1 + 0.15 * randn(size(widths))), 6);
    etas = clampValue(etas + 0.10 * randn(size(etas)), 0.15, 0.95);
    baseline = max(0, baseline + 0.003 * randn());
    perturbedGuess = [weights; centers; widths; etas; baseline];
end

function fitConstraints = buildCandidateFitConstraints(peakCandidates)
    peakCandidates = double(peakCandidates);
    [~, order] = sort(peakCandidates(:, 2), 'ascend');
    peakCandidates = peakCandidates(order, :);
    centersNm = peakCandidates(:, 2);
    lowerNm = centersNm - 45;
    upperNm = centersNm + 45;
    targetNm = centersNm;
    stiffness = 0.0008 * ones(size(centersNm));

    shortMask = centersNm >= 372 & centersNm <= 398;
    lowerNm(shortMask) = 376;
    upperNm(shortMask) = 396;
    targetNm(shortMask) = clampValue(centersNm(shortMask), 380, 391);
    stiffness(shortMask) = 0.0065;

    shoulderMask = centersNm >= 397 & centersNm <= 434;
    lowerNm(shoulderMask) = 398;
    upperNm(shoulderMask) = 434;
    targetNm(shoulderMask) = clampValue(centersNm(shoulderMask), 405, 425);
    stiffness(shoulderMask) = 0.0040;

    barrierMask = centersNm >= 438 & centersNm <= 480;
    lowerNm(barrierMask) = 438;
    upperNm(barrierMask) = 482;
    targetNm(barrierMask) = clampValue(centersNm(barrierMask), 448, 468);
    stiffness(barrierMask) = 0.0020;

    mainMask = centersNm >= 485 & centersNm <= 560;
    lowerNm(mainMask) = 488;
    upperNm(mainMask) = 560;
    stiffness(mainMask) = 0.0004;

    fitConstraints = struct( ...
        'lowerNm', lowerNm(:), ...
        'upperNm', upperNm(:), ...
        'targetNm', targetNm(:), ...
        'stiffness', stiffness(:));
end

function penalty = computeCenterConstraintPenalty(componentsNm, fitConstraints)
    if isempty(componentsNm) || numel(fitConstraints.targetNm) ~= size(componentsNm, 1)
        penalty = 0;
        return;
    end

    centersNm = componentsNm(:, 2);
    lowerExcess = max(fitConstraints.lowerNm - centersNm, 0);
    upperExcess = max(centersNm - fitConstraints.upperNm, 0);
    outsideDistance = lowerExcess + upperExcess;
    targetDistance = centersNm - fitConstraints.targetNm;
    penalty = 0.12 * sum(outsideDistance .^ 2) + ...
        sum(fitConstraints.stiffness .* targetDistance .^ 2);
end

function errorValue = fitMixtureObjective(params, peakCount, axisNm, observedSpectrum, fitConstraints)
    [modelSpectrum, componentsNm, ~, baseline] = evaluateGaussianMixtureModel(params, peakCount, axisNm);
    normalizedModel = normalizeSpectrum(modelSpectrum);
    errorValue = computeWeightedFitError(observedSpectrum, normalizedModel, axisNm);

    centerSpacing = diff(componentsNm(:,2));
    if any(centerSpacing < 6)
        errorValue = errorValue + 0.30 * sum(6 - centerSpacing(centerSpacing < 6));
    end
    if any(componentsNm(:,3) < 5)
        errorValue = errorValue + 0.10 * sum(5 - componentsNm(componentsNm(:,3) < 5, 3));
    end
    if size(componentsNm, 2) >= 4 && any(componentsNm(:, 4) < 0.05 | componentsNm(:, 4) > 0.98)
        errorValue = errorValue + 0.05 * sum(abs(componentsNm(:, 4) - clampValue(componentsNm(:, 4), 0.05, 0.98)));
    end
    errorValue = errorValue + 1.8 * baseline;
    if baseline > 0.045
        errorValue = errorValue + 8 * (baseline - 0.045);
    end

    backgroundMask = observedSpectrum < 0.025;
    if any(backgroundMask)
        backgroundOvershoot = max(normalizedModel(backgroundMask) - 0.04, 0);
        if any(backgroundOvershoot > 0)
            errorValue = errorValue + 0.75 * mean(backgroundOvershoot .^ 2) / max(mean(backgroundMask), eps);
        end
    end

    valleyMask = axisNm >= 405 & axisNm <= 490;
    if any(valleyMask)
        valleyResidual = abs(observedSpectrum(valleyMask) - normalizedModel(valleyMask));
        errorValue = errorValue + 0.45 * mean(valleyResidual .^ 2);
    end

    targetResidual = normalizeSpectrum(max(observedSpectrum - normalizedModel, 0));
    if any(targetResidual > 0.03)
        errorValue = errorValue + 0.30 * mean(targetResidual(targetResidual > 0.03));
    end

    errorValue = errorValue + computePhysicsWindowRatioPenalty(axisNm, observedSpectrum, normalizedModel);
    errorValue = errorValue + computeMainWindowShapePenalty(axisNm, observedSpectrum, normalizedModel);
    if nargin >= 5 && ~isempty(fitConstraints)
        errorValue = errorValue + computeCenterConstraintPenalty(componentsNm, fitConstraints);
    end
end

function [modelSpectrum, componentsNm, componentSpectra, baseline] = evaluateGaussianMixtureModel(params, peakCount, axisNm)
    rawWeights = abs(params(1:peakCount));
    if ~any(rawWeights)
        rawWeights = ones(peakCount, 1);
    end
    weights = rawWeights / sum(rawWeights);
    centersNm = params(peakCount+1:2*peakCount);
    fwhmNm = abs(params(2*peakCount+1:3*peakCount));
    etaValues = clampValue(abs(params(3*peakCount+1:4*peakCount)), 0.05, 0.95);
    baseline = clampValue(abs(params(end)), 0, 0.06);

    [centersNm, order] = sort(centersNm(:), 'ascend');
    weights = weights(order);
    fwhmNm = clampValue(fwhmNm(order), 5, 180);
    etaValues = etaValues(order);

    componentSpectra = zeros(numel(axisNm), peakCount);
    for peakIdx = 1:peakCount
        componentSpectra(:, peakIdx) = weights(peakIdx) .* evaluatePseudoVoigtProfile( ...
            axisNm, centersNm(peakIdx), fwhmNm(peakIdx), etaValues(peakIdx));
    end
    modelSpectrum = sum(componentSpectra, 2) + baseline;
    scaleValue = max(modelSpectrum);
    if scaleValue > 0
        componentSpectra = componentSpectra ./ scaleValue;
    end
    componentsNm = [weights(:), centersNm(:), fwhmNm(:), etaValues(:)];
end

function profileValues = evaluatePseudoVoigtProfile(axisNm, centerNm, fwhmNm, eta)
    sigmaNm = fwhmToSigma(fwhmNm);
    gaussianPart = exp(-0.5 * ((axisNm - centerNm) ./ sigmaNm) .^ 2);
    lorentzianPart = 1 ./ (1 + 4 * ((axisNm - centerNm) ./ fwhmNm) .^ 2);
    profileValues = eta .* gaussianPart + (1 - eta) .* lorentzianPart;
end

function peakWidthNm = estimatePeakWidthNm(axisNm, spectrum, peakIdx)
    peakHeight = spectrum(peakIdx);
    halfHeight = peakHeight * 0.5;
    leftIdx = peakIdx;
    rightIdx = peakIdx;
    while leftIdx > 1 && spectrum(leftIdx) > halfHeight
        leftIdx = leftIdx - 1;
    end
    while rightIdx < numel(spectrum) && spectrum(rightIdx) > halfHeight
        rightIdx = rightIdx + 1;
    end
    peakWidthNm = max(8, axisNm(rightIdx) - axisNm(leftIdx));
end

function significantComponents = selectSignificantComponents(componentsNm)
    maxWeight = max(componentsNm(:,1));
    physicalMask = componentsNm(:, 2) >= 360 & componentsNm(:, 2) <= 680;
    broadSignificantMask = componentsNm(:,1) >= max(0.020, 0.06 * maxWeight);
    shortwaveWindowMask = componentsNm(:, 2) >= 370 & componentsNm(:, 2) <= 432;
    shortwaveMask = shortwaveWindowMask & componentsNm(:,1) >= max(0.004, 0.010 * maxWeight);
    shoulderWindowMask = componentsNm(:, 2) >= 395 & componentsNm(:, 2) <= 505;
    shoulderMask = shoulderWindowMask & componentsNm(:,1) >= max(0.012, 0.035 * maxWeight);
    significantMask = physicalMask & (broadSignificantMask | shoulderMask | shortwaveMask);
    significantComponents = componentsNm(significantMask, :);
    if isempty(significantComponents)
        significantComponents = componentsNm;
    end
end

function weights = normalizePositive(values)
    weights = max(double(values(:)), eps);
    weights = weights / sum(weights);
end

function spectrumModel = buildGaussianMixtureSpectrum(componentsNm)
    etaColumnNm = ensureEtaColumnNm(componentsNm);
    spectrumModel = struct( ...
        'modelType', 'pseudo_voigt_mixture', ...
        'components', [componentsNm(:,1), componentsNm(:,2) * 1e-9, componentsNm(:,3) * 1e-9, etaColumnNm]);
end

function assignments = assignComponentsToEmissionFamilies(componentsNm, familyDefaults, observedMainPeakNm)
    if nargin < 3 || ~isfinite(observedMainPeakNm)
        observedMainPeakNm = NaN;
    end
    componentsNm = normalizeComponentWeights(componentsNm);
    componentCount = size(componentsNm, 1);
    assignedMask = false(componentCount, 1);

    nGaNRef = dominantComponentFromGroup(familyDefaults.nGaN);
    prestrainedRef = dominantComponentFromGroup(familyDefaults.prestrained);
    barrierRef = dominantComponentFromGroup(familyDefaults.barrier);
    wellRef = dominantComponentFromGroup(familyDefaults.well);
    pTypeRef = dominantComponentFromGroup(familyDefaults.pType);

    nGaNCandidates = find(componentsNm(:, 2) >= 372 & componentsNm(:, 2) <= 397);
    if isempty(nGaNCandidates)
        nGaNIdx = selectClosestUnassignedComponent(componentsNm, assignedMask, nGaNRef(2), [372, 397]);
    else
        [~, strongestIdx] = max(componentsNm(nGaNCandidates, 1));
        nGaNIdx = nGaNCandidates(strongestIdx);
    end
    if ~isempty(nGaNIdx)
        assignedMask(nGaNIdx) = true;
    end

    wellCandidates = find(~assignedMask & componentsNm(:, 2) >= 488 & componentsNm(:, 2) <= 560);
    if isfinite(observedMainPeakNm)
        wellIdx = selectClosestUnassignedComponent(componentsNm, assignedMask, observedMainPeakNm, [488, 560]);
    elseif isempty(wellCandidates)
        wellIdx = selectClosestUnassignedComponent(componentsNm, assignedMask, wellRef(2), [488, 590]);
    else
        [~, strongestWellIdx] = max(componentsNm(wellCandidates, 1));
        wellIdx = wellCandidates(strongestWellIdx);
    end
    if isempty(wellIdx)
        wellIdx = componentCount;
    end
    assignedMask(wellIdx) = true;

    semipolarIdx = [];
    if componentCount >= 4
        wellCenterNm = componentsNm(wellIdx, 2);
        redTailCandidates = find(~assignedMask & componentsNm(:, 2) >= wellCenterNm + 8 & componentsNm(:, 2) <= 570);
        if ~isempty(redTailCandidates)
            [~, strongestTailIdx] = max(componentsNm(redTailCandidates, 1));
            semipolarIdx = redTailCandidates(strongestTailIdx);
        else
            semipolarTarget = max(barrierRef(2) + 22, wellCenterNm - 30);
            semipolarIdx = selectClosestUnassignedComponent(componentsNm, assignedMask, semipolarTarget, [470, wellCenterNm - 6]);
        end
    end
    if ~isempty(semipolarIdx)
        assignedMask(semipolarIdx) = true;
    end

    prestrainedIdx = selectClosestUnassignedComponent(componentsNm, assignedMask, prestrainedRef(2), [398, 455]);
    if ~isempty(prestrainedIdx)
        assignedMask(prestrainedIdx) = true;
    end

    barrierCandidates = find(~assignedMask & componentsNm(:, 2) >= 430 & componentsNm(:, 2) <= max(505, wellRef(2) - 18));
    if isempty(barrierCandidates)
        barrierIdx = selectClosestUnassignedComponent(componentsNm, assignedMask, barrierRef(2), [430, 515]);
    else
        barrierIdx = barrierCandidates(:)';
    end
    assignedMask(barrierIdx) = true;

    pTypeIdx = selectClosestUnassignedComponent(componentsNm, assignedMask, pTypeRef(2), [376, 397]);

    nGaNComponentsNm = [];
    if ~isempty(nGaNIdx)
        nGaNComponentsNm = normalizeComponentWeights(componentsNm(nGaNIdx, :));
    end
    if isempty(nGaNComponentsNm)
        nGaNComponentsNm = normalizeComponentWeights(familyDefaults.nGaN);
    end

    if isempty(prestrainedIdx)
        prestrainedComponentsNm = derivePrestrainedComponents(nGaNComponentsNm, barrierRef);
    else
        prestrainedComponentsNm = normalizeComponentWeights(componentsNm(prestrainedIdx, :));
    end

    wellComponentsNm = normalizeComponentWeights(componentsNm(wellIdx, :));
    if isempty(wellComponentsNm)
        wellComponentsNm = normalizeComponentWeights(familyDefaults.well);
    end

    if isempty(barrierIdx)
        barrierComponentsNm = deriveFallbackBarrierComponents(familyDefaults.barrier, prestrainedComponentsNm, wellComponentsNm);
    else
        barrierComponentsNm = normalizeComponentWeights(componentsNm(barrierIdx, :));
    end

    if isempty(semipolarIdx)
        semipolarComponentsNm = deriveSemipolarComponents(wellComponentsNm, NaN);
    else
        semipolarComponentsNm = normalizeComponentWeights(componentsNm(semipolarIdx, :));
    end

    if isempty(pTypeIdx)
        pTypeComponentsNm = derivePTypeComponents(prestrainedComponentsNm, barrierComponentsNm, familyDefaults.pType);
    else
        pTypeComponentsNm = normalizeComponentWeights(componentsNm(pTypeIdx, :));
    end

    assignments = struct( ...
        'nGaNComponentsNm', nGaNComponentsNm, ...
        'prestrainedComponentsNm', prestrainedComponentsNm, ...
        'barrierComponentsNm', barrierComponentsNm, ...
        'wellComponentsNm', wellComponentsNm, ...
        'semipolarComponentsNm', semipolarComponentsNm, ...
        'pTypeComponentsNm', pTypeComponentsNm);
end

function componentsNm = normalizeComponentWeights(componentsNm)
    if isempty(componentsNm)
        return;
    end
    componentsNm = double(componentsNm);
    componentsNm(:, 1) = normalizePositive(componentsNm(:, 1));
    [~, order] = sort(componentsNm(:, 2), 'ascend');
    componentsNm = componentsNm(order, :);
end

function component = dominantComponentFromGroup(componentsNm)
    componentsNm = normalizeComponentWeights(componentsNm);
    [~, dominantIdx] = max(componentsNm(:, 1));
    component = componentsNm(dominantIdx, :);
end

function barrierComponentsNm = deriveFallbackBarrierComponents(defaultBarrierSpectrum, prestrainedComponentsNm, wellComponentsNm)
    defaultBarrierComponents = normalizeComponentWeights(defaultBarrierSpectrum);
    prestrainedPeak = dominantComponentFromGroup(prestrainedComponentsNm);
    wellPeak = dominantComponentFromGroup(wellComponentsNm);
    barrierComponentsNm = defaultBarrierComponents(1:min(2, size(defaultBarrierComponents, 1)), :);
    centerNm = max(prestrainedPeak(2) + 35, wellPeak(2) - 58);
    widthNm = max(14, 0.85 * wellPeak(3));
    if isempty(barrierComponentsNm)
        barrierComponentsNm = [1, centerNm, widthNm, 0.72];
    else
        barrierComponentsNm(1, 1:3) = [0.75, centerNm, widthNm];
        if size(barrierComponentsNm, 1) >= 2
            barrierComponentsNm(2, 1:3) = [0.25, centerNm + 18, max(widthNm, barrierComponentsNm(2, 3))];
        end
    end
    barrierComponentsNm = normalizeComponentWeights(barrierComponentsNm);
end

function semipolarComponentsNm = deriveSemipolarComponents(wellComponentsNm, wellInComposition)
    semipolarComponentsNm = normalizeComponentWeights(wellComponentsNm);
    if ~isfinite(wellInComposition)
        compositionShift_eV = 0.065;
    else
        semipolarInComposition = clampValue(0.85 * wellInComposition, 0.02, 0.95);
        compositionShift_eV = max(0, bandgapFromComposition(semipolarInComposition) - bandgapFromComposition(wellInComposition));
    end
    confinementShift_eV = 0.025;
    totalShift_eV = compositionShift_eV + confinementShift_eV;

    photonEnergy_eV = 1240 ./ semipolarComponentsNm(:, 2);
    shiftedCentersNm = 1240 ./ max(photonEnergy_eV + totalShift_eV, 0.05);
    semipolarComponentsNm(:, 2) = shiftedCentersNm;
    semipolarComponentsNm(:, 3) = max(8, 0.9 * semipolarComponentsNm(:, 3));
    semipolarComponentsNm = normalizeComponentWeights(semipolarComponentsNm);
end

function prestrainedComponentsNm = derivePrestrainedComponents(nGaNComponentsNm, barrierRef)
    baseComponent = dominantComponentFromGroup(nGaNComponentsNm);
    upperBoundNm = max(405, barrierRef(2) - 18);
    centerNm = min(upperBoundNm, max(392, baseComponent(2) + 24));
    widthNm = max(12, 1.8 * baseComponent(3));
    etaValue = getEtaFromComponent(baseComponent);
    prestrainedComponentsNm = normalizeComponentWeights([1, centerNm, widthNm, etaValue]);
end

function pTypeComponentsNm = derivePTypeComponents(~, ~, defaultPTypeSpectrum)
    defaultPTypeComponents = normalizeComponentWeights(defaultPTypeSpectrum);
    if isempty(defaultPTypeComponents)
        centerNm = 388;
        widthNm = 12;
        etaValue = 0.78;
    else
        defaultPeak = dominantComponentFromGroup(defaultPTypeComponents);
        centerNm = clampValue(defaultPeak(2), 380, 398);
        widthNm = clampValue(defaultPeak(3), 8, 18);
        etaValue = max(getEtaFromComponent(defaultPeak), 0.70);
    end
    pTypeComponentsNm = defaultPTypeComponents(1:min(2, size(defaultPTypeComponents, 1)), :);
    if isempty(pTypeComponentsNm)
        pTypeComponentsNm = [1, centerNm, widthNm, etaValue];
    else
        pTypeComponentsNm(1, 1:4) = [0.76, centerNm, widthNm, etaValue];
        if size(pTypeComponentsNm, 1) >= 2
            pTypeComponentsNm(2, 1:4) = [0.24, centerNm + 9, max(widthNm, pTypeComponentsNm(2, 3)), etaValue];
        end
    end
    pTypeComponentsNm = normalizeComponentWeights(pTypeComponentsNm);
end

function selectedIdx = selectClosestUnassignedComponent(componentsNm, assignedMask, targetCenterNm, validRangeNm)
    candidateMask = ~assignedMask;
    if nargin >= 4 && ~isempty(validRangeNm)
        candidateMask = candidateMask & componentsNm(:, 2) >= validRangeNm(1) & componentsNm(:, 2) <= validRangeNm(2);
    end
    candidateIdx = find(candidateMask);
    if isempty(candidateIdx)
        selectedIdx = [];
        return;
    end
    [~, bestIdx] = min(abs(componentsNm(candidateIdx, 2) - targetCenterNm));
    selectedIdx = candidateIdx(bestIdx);
end

function Eg = bandgapFromComposition(inComposition)
    EgGaN = 3.3032;
    EgInN = 0.6086;
    Eg = EgInN * inComposition + EgGaN * (1 - inComposition) - 1.43 * inComposition * (1 - inComposition);
end

function barrierSpectrum = deriveBarrierSpectrum(defaultBarrierSpectrum, fittedComponentsNm)
    legacyComponents = convertSpectrumModelToComponentsNm(defaultBarrierSpectrum);
    selectedCount = min(2, size(fittedComponentsNm, 1));
    selectedComponents = fittedComponentsNm(1:selectedCount, :);
    barrierComponents = legacyComponents(1:min(size(legacyComponents, 1), selectedCount), :);
    if isempty(barrierComponents)
        barrierComponents = selectedComponents;
    end

    for idx = 1:size(barrierComponents, 1)
        fitComp = selectedComponents(min(idx, size(selectedComponents, 1)), :);
        barrierComponents(idx, 1) = max(0.05, fitComp(1));
        barrierComponents(idx, 2) = max(200, 0.65 * fitComp(2) + 0.35 * barrierComponents(idx, 2));
        barrierComponents(idx, 3) = max(barrierComponents(idx, 3), 0.85 * fitComp(3));
    end

    barrierSpectrum = componentsNmToLegacySpectrum(barrierComponents);
end

function [meanSpectrum, validCount] = buildMeanNormalizedSpectrum(spectraCell, roi, commonAxisNm)
    spectraSum = zeros(numel(commonAxisNm), 1);
    validCount = 0;

    for rowIdx = roi(1):roi(2)
        for colIdx = roi(3):roi(4)
            thisSpectrum = spectraCell{rowIdx, colIdx};
            if isempty(thisSpectrum)
                continue;
            end

            normalizedSpectrum = preprocessSpectrum(thisSpectrum, commonAxisNm);
            if ~any(normalizedSpectrum)
                continue;
            end

            spectraSum = spectraSum + normalizedSpectrum;
            validCount = validCount + 1;
        end
    end

    if validCount < 1
        meanSpectrum = zeros(numel(commonAxisNm), 1);
        return;
    end

    meanSpectrum = spectraSum / validCount;
    meanSpectrum = normalizeSpectrum(meanSpectrum);
end

function [meanSpectrum, validCount] = buildMeanNormalizedSpectrumFromMask(spectraCell, pixelMask, commonAxisNm)
    spectraSum = zeros(numel(commonAxisNm), 1);
    validCount = 0;
    [rows, cols] = size(spectraCell);
    pixelMask = logical(pixelMask);
    for rowIdx = 1:rows
        for colIdx = 1:cols
            if rowIdx > size(pixelMask, 1) || colIdx > size(pixelMask, 2) || ~pixelMask(rowIdx, colIdx)
                continue;
            end
            thisSpectrum = spectraCell{rowIdx, colIdx};
            if isempty(thisSpectrum)
                continue;
            end

            normalizedSpectrum = preprocessSpectrum(thisSpectrum, commonAxisNm);
            if ~any(normalizedSpectrum)
                continue;
            end

            spectraSum = spectraSum + normalizedSpectrum;
            validCount = validCount + 1;
        end
    end

    if validCount < 1
        meanSpectrum = zeros(numel(commonAxisNm), 1);
        return;
    end

    meanSpectrum = spectraSum / validCount;
    meanSpectrum = normalizeSpectrum(meanSpectrum);
end

function selection = buildCalibrationSelection(experimentalData, roi)
    commonAxisNm = getCommonAxisNm(experimentalData);
    [meanSpectrum, validCount] = buildMeanNormalizedSpectrum(experimentalData.totalSpectra, roi, commonAxisNm);
    selection = struct( ...
        'roi', roi, ...
        'axisNm', commonAxisNm, ...
        'meanSpectrum', meanSpectrum, ...
        'validCount', validCount, ...
        'sourceFile', getFieldOrDefault(experimentalData, 'sourceFile', ''), ...
        'sourceName', getFieldOrDefault(experimentalData, 'sourceName', 'Imported TIFF'));
end

function normalizedSpectrum = preprocessSpectrum(spectrum, commonAxisNm)
    normalizedSpectrum = normalizeSpectrum(preprocessSpectrumForIntensity(spectrum, commonAxisNm));
end

function processedSpectrum = preprocessSpectrumForIntensity(spectrum, commonAxisNm)
    processedSpectrum = zeros(numel(commonAxisNm), 1);
    if isempty(spectrum) || size(spectrum, 2) < 2
        return;
    end

    wavelengthNm = double(spectrum(:, 1));
    intensity = double(spectrum(:, 2));
    validMask = isfinite(wavelengthNm) & isfinite(intensity);
    wavelengthNm = wavelengthNm(validMask);
    intensity = intensity(validMask);
    if isempty(wavelengthNm)
        return;
    end

    [wavelengthNm, sortIdx] = sort(wavelengthNm);
    intensity = intensity(sortIdx);
    [uniqueWavelengths, ~, groupIdx] = unique(wavelengthNm);
    if numel(uniqueWavelengths) ~= numel(wavelengthNm)
        intensity = accumarray(groupIdx, intensity, [], @mean);
        wavelengthNm = uniqueWavelengths;
    end

    if isscalar(wavelengthNm)
        [~, closestIdx] = min(abs(commonAxisNm - wavelengthNm));
        interpolated = zeros(size(commonAxisNm));
        interpolated(closestIdx) = max(intensity, 0);
    else
        interpolated = interp1(wavelengthNm, intensity, commonAxisNm, 'linear', 0);
    end

    interpolated = max(interpolated, 0);
    if numel(interpolated) >= 5
        interpolated = movmean(interpolated, 5);
    end
    baseline = prctile(interpolated, 5);
    processedSpectrum = max(interpolated - baseline, 0);
end

function commonAxisNm = getCommonAxisNm(data)
    if isfield(data, 'wavelengthAxis') && ~isempty(data.wavelengthAxis)
        commonAxisNm = clipCalibrationAxisNm(double(data.wavelengthAxis(:)));
        return;
    end

    commonAxisNm = [];
    if isfield(data, 'totalSpectra') && ~isempty(data.totalSpectra)
        for idx = 1:numel(data.totalSpectra)
            spectrum = data.totalSpectra{idx};
            if ~isempty(spectrum) && size(spectrum, 2) >= 2
                commonAxisNm = clipCalibrationAxisNm(double(spectrum(:, 1)));
                break;
            end
        end
    end

    if isempty(commonAxisNm)
        commonAxisNm = linspace(350, 700, 512)';
    end
end

function axisNm = clipCalibrationAxisNm(axisNm)
    axisNm = double(axisNm(:));
    axisNm = axisNm(isfinite(axisNm));
    if isempty(axisNm)
        axisNm = linspace(350, 700, 512)';
        return;
    end

    validMask = axisNm >= 350 & axisNm <= 700;
    if nnz(validMask) >= 64
        axisNm = axisNm(validMask);
    end

    axisNm = unique(axisNm, 'stable');
    if numel(axisNm) < 64
        axisNm = linspace(max(350, min(axisNm)), min(700, max(axisNm)), 512)';
    end
end

function refinement = buildLocalRefinement(experimentalData, simulationData, commonAxisNm, roi, profile)
    expSize = size(experimentalData.totalSpectra);
    simSize = size(simulationData.totalSpectra);
    rows = min(expSize(1), simSize(1));
    cols = min(expSize(2), simSize(2));
    nGaNGainMap = ones(rows, cols);
    prestrainedGainMap = ones(rows, cols);
    barrierGainMap = ones(rows, cols);
    wellGainMap = ones(rows, cols);
    semipolarGainMap = ones(rows, cols);
    pTypeGainMap = ones(rows, cols);

    nGaNBand = buildPeakBand(getFieldOrDefault(profile, 'shortwavePeakNm', getFieldOrDefault(profile, 'nGaNPeakNm', 390)), 11, [376, 396]);
    prestrainedBand = buildPeakBand(getFieldOrDefault(profile, 'shoulderPeakNm', getFieldOrDefault(profile, 'prestrainedPeakNm', 413)), 18, [398, 432]);
    barrierBand = buildPeakBand(getFieldOrDefault(profile, 'barrierPeakNm', 459), 18, [442, 475]);
    semipolarBand = buildPeakBand(getFieldOrDefault(profile, 'semipolarPeakNm', 510), 24, [470, 545]);
    wellBand = buildPeakBand(getFieldOrDefault(profile, 'wellPeakNm', getFieldOrDefault(profile, 'redPeakNm', 540)), 28, [490, 650]);
    pTypeBand = buildPeakBand(getFieldOrDefault(profile, 'pTypePeakNm', 388), 11, [376, 396]);

    for rowIdx = roi(1):roi(2)
        for colIdx = roi(3):roi(4)
            experimentalSpectrum = experimentalData.totalSpectra{rowIdx, colIdx};
            simulationSpectrum = simulationData.totalSpectra{rowIdx, colIdx};
            if isempty(experimentalSpectrum) || isempty(simulationSpectrum)
                continue;
            end

            yExperimental = preprocessSpectrum(experimentalSpectrum, commonAxisNm);
            ySimulation = preprocessSpectrum(simulationSpectrum, commonAxisNm);
            if ~any(yExperimental) || ~any(ySimulation)
                continue;
            end

            nGaNGainMap(rowIdx, colIdx) = computeBandRatioGain(commonAxisNm, yExperimental, ySimulation, nGaNBand, 0.55, 1.45);
            prestrainedGainMap(rowIdx, colIdx) = computeBandRatioGain(commonAxisNm, yExperimental, ySimulation, prestrainedBand, 0.50, 1.50);
            barrierGainMap(rowIdx, colIdx) = computeBandRatioGain(commonAxisNm, yExperimental, ySimulation, barrierBand, 0.40, 1.35);
            semipolarGainMap(rowIdx, colIdx) = computeBandRatioGain(commonAxisNm, yExperimental, ySimulation, semipolarBand, 0.35, 1.20);
            wellGainMap(rowIdx, colIdx) = computeBandRatioGain(commonAxisNm, yExperimental, ySimulation, wellBand, 0.70, 1.35);
            pTypeGainMap(rowIdx, colIdx) = computeBandRatioGain(commonAxisNm, yExperimental, ySimulation, pTypeBand, 0.45, 1.45);
        end
    end

    refinement = struct( ...
        'version', 1, ...
        'createdAt', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
        'roi', roi, ...
        'nGaNGainMap', nGaNGainMap, ...
        'gaNGainMap', nGaNGainMap, ...
        'prestrainedGainMap', prestrainedGainMap, ...
        'barrierGainMap', barrierGainMap, ...
        'wellGainMap', wellGainMap, ...
        'semipolarGainMap', semipolarGainMap, ...
        'pTypeGainMap', pTypeGainMap, ...
        'nGaNBand', nGaNBand, ...
        'gaNBand', nGaNBand, ...
        'prestrainedBand', prestrainedBand, ...
        'barrierBand', barrierBand, ...
        'wellBand', wellBand, ...
        'semipolarBand', semipolarBand, ...
        'pTypeBand', pTypeBand);
end

function localField = buildExperimentalLocalParameterField(experimentalData, commonAxisNm, profile)
    dataSize = size(experimentalData.totalSpectra);
    rows = dataSize(1);
    cols = dataSize(2);
    globalWellIn = profile.wellInComposition;
    globalBarrierIn = profile.barrierInComposition;
    globalSemipolarIn = invertInComposition(getFieldOrDefault(profile, 'semipolarPeakNm', profile.barrierPeakNm), 0.85 * globalWellIn);

    shortwaveBand = buildPeakBand(getFieldOrDefault(profile, 'shortwavePeakNm', getFieldOrDefault(profile, 'nGaNPeakNm', 390)), 11, [376, 396]);
    shoulderBand = buildPeakBand(getFieldOrDefault(profile, 'shoulderPeakNm', getFieldOrDefault(profile, 'prestrainedPeakNm', 413)), 18, [398, 432]);
    barrierBand = buildPeakBand(getFieldOrDefault(profile, 'barrierPeakNm', 459), 18, [442, 475]);
    semipolarBand = buildPeakBand(getFieldOrDefault(profile, 'semipolarPeakNm', 510), 24, [470, 545]);
    wellBand = buildPeakBand(getFieldOrDefault(profile, 'wellPeakNm', getFieldOrDefault(profile, 'redPeakNm', 540)), 30, [490, 650]);
    pTypeBand = shortwaveBand;

    nGaNAreaMap = nan(rows, cols);
    prestrainedAreaMap = nan(rows, cols);
    barrierAreaMap = nan(rows, cols);
    semipolarAreaMap = nan(rows, cols);
    wellAreaMap = nan(rows, cols);
    pTypeAreaMap = nan(rows, cols);
    totalIntensityMap = nan(rows, cols);
    nGaNPeakMap_nm = nan(rows, cols);
    prestrainedPeakMap_nm = nan(rows, cols);
    barrierPeakMap_nm = nan(rows, cols);
    semipolarPeakMap_nm = nan(rows, cols);
    wellPeakMap_nm = nan(rows, cols);
    pTypePeakMap_nm = nan(rows, cols);

    for rowIdx = 1:rows
        for colIdx = 1:cols
            spectrum = experimentalData.totalSpectra{rowIdx, colIdx};
            intensitySpectrum = preprocessSpectrumForIntensity(spectrum, commonAxisNm);
            normalizedSpectrum = normalizeSpectrum(intensitySpectrum);
            if ~any(normalizedSpectrum)
                continue;
            end

            totalIntensityMap(rowIdx, colIdx) = trapz(commonAxisNm, intensitySpectrum);
            nGaNAreaMap(rowIdx, colIdx) = integrateBandArea(commonAxisNm, normalizedSpectrum, shortwaveBand);
            prestrainedAreaMap(rowIdx, colIdx) = integrateBandArea(commonAxisNm, normalizedSpectrum, shoulderBand);
            barrierAreaMap(rowIdx, colIdx) = integrateBandArea(commonAxisNm, normalizedSpectrum, barrierBand);
            semipolarAreaMap(rowIdx, colIdx) = integrateBandArea(commonAxisNm, normalizedSpectrum, semipolarBand);
            wellAreaMap(rowIdx, colIdx) = integrateBandArea(commonAxisNm, normalizedSpectrum, wellBand);
            pTypeAreaMap(rowIdx, colIdx) = integrateBandArea(commonAxisNm, normalizedSpectrum, pTypeBand);
            nGaNPeakMap_nm(rowIdx, colIdx) = computeBandCentroid(commonAxisNm, normalizedSpectrum, shortwaveBand, getFieldOrDefault(profile, 'shortwavePeakNm', getFieldOrDefault(profile, 'nGaNPeakNm', 390)));
            prestrainedPeakMap_nm(rowIdx, colIdx) = computeBandCentroid(commonAxisNm, normalizedSpectrum, shoulderBand, getFieldOrDefault(profile, 'shoulderPeakNm', getFieldOrDefault(profile, 'prestrainedPeakNm', 413)));
            barrierPeakMap_nm(rowIdx, colIdx) = computeBandCentroid(commonAxisNm, normalizedSpectrum, barrierBand, getFieldOrDefault(profile, 'barrierPeakNm', 405));
            semipolarPeakMap_nm(rowIdx, colIdx) = computeBandCentroid(commonAxisNm, normalizedSpectrum, semipolarBand, getFieldOrDefault(profile, 'semipolarPeakNm', 510));
            wellPeakMap_nm(rowIdx, colIdx) = computeBandCentroid(commonAxisNm, normalizedSpectrum, wellBand, getFieldOrDefault(profile, 'wellPeakNm', getFieldOrDefault(profile, 'redPeakNm', 540)));
            pTypePeakMap_nm(rowIdx, colIdx) = computeBandCentroid(commonAxisNm, normalizedSpectrum, pTypeBand, getFieldOrDefault(profile, 'pTypePeakNm', 388));
        end
    end

    shortwaveAreaMap = nGaNAreaMap;
    rawShoulderAreaMap = prestrainedAreaMap;
    shoulderAreaMap = max(rawShoulderAreaMap - 0.28 * shortwaveAreaMap - 0.035 * wellAreaMap, 0);
    prestrainedAreaMap = shoulderAreaMap;
    [vpitInfluenceMap, vpitCoreMask, vpitRimMask, regionCategoryMap] = ...
        buildVpitRegionMaps(experimentalData, rows, cols, shortwaveAreaMap, shoulderAreaMap, wellAreaMap);
    shortMainRatioMap = shortwaveAreaMap ./ max(wellAreaMap, eps);
    shoulderMainRatioMap = shoulderAreaMap ./ max(wellAreaMap, eps);
    pTypeShareMap = clampValue( ...
        0.34 + 0.30 * vpitInfluenceMap + 0.10 * normalizeFiniteMap(shoulderMainRatioMap) - ...
        0.08 * normalizeFiniteMap(wellAreaMap), 0.24, 0.72);
    pTypeAreaMap = shortwaveAreaMap .* pTypeShareMap;
    nGaNAreaMap = shortwaveAreaMap .* (1 - pTypeShareMap);
    pTypeReferencePeak_nm = getFieldOrDefault(profile, 'pTypePeakNm', 388.5);
    pTypePeakMap_nm = min(nGaNPeakMap_nm, pTypeReferencePeak_nm);

    nGaNRef = medianPositive(nGaNAreaMap, 0.08);
    prestrainedRef = medianPositive(prestrainedAreaMap, 0.08);
    barrierRef = medianPositive(barrierAreaMap, 0.12);
    semipolarRef = medianPositive(semipolarAreaMap, 0.16);
    wellRef = medianPositive(wellAreaMap, 0.22);
    pTypeRef = medianPositive(pTypeAreaMap, 0.10);
    totalIntensityRef = medianPositive(totalIntensityMap, 1.0);
    semipolarToWellRef = medianPositive(semipolarAreaMap ./ max(wellAreaMap, eps), 0.75);

    wellInMap = globalWellIn * ones(rows, cols);
    barrierInMap = globalBarrierIn * ones(rows, cols);
    semipolarInMap = globalSemipolarIn * ones(rows, cols);
    shellThicknessScaleMap = ones(rows, cols);
    totalIntensityGainMap = ones(rows, cols);
    vpitCollectionGainMap = ones(rows, cols);
    centerFillGainMap = ones(rows, cols);
    nGaNGainMap = ones(rows, cols);
    prestrainedGainMap = ones(rows, cols);
    barrierGainMap = ones(rows, cols);
    wellGainMap = ones(rows, cols);
    semipolarGainMap = ones(rows, cols);
    pTypeGainMap = ones(rows, cols);

    localIntensityContext = computeLocalIntensityContext(totalIntensityMap);

    for rowIdx = 1:rows
        for colIdx = 1:cols
            if isfinite(wellPeakMap_nm(rowIdx, colIdx))
                localWellIn = invertInComposition(wellPeakMap_nm(rowIdx, colIdx), globalWellIn);
                wellInMap(rowIdx, colIdx) = clampValue(localWellIn, globalWellIn - 0.08, globalWellIn + 0.08);
            end

            if isfinite(barrierPeakMap_nm(rowIdx, colIdx))
                localBarrierIn = invertInComposition(barrierPeakMap_nm(rowIdx, colIdx), globalBarrierIn);
                barrierInMap(rowIdx, colIdx) = clampValue(localBarrierIn, 0.01, min(wellInMap(rowIdx, colIdx) - 0.02, globalBarrierIn + 0.06));
            end

            if isfinite(semipolarPeakMap_nm(rowIdx, colIdx))
                localSemipolarIn = invertInComposition(semipolarPeakMap_nm(rowIdx, colIdx), globalSemipolarIn);
                semipolarInMap(rowIdx, colIdx) = clampValue(localSemipolarIn, 0.02, max(0.03, wellInMap(rowIdx, colIdx) - 0.015));
            end

            semipolarRatio = semipolarAreaMap(rowIdx, colIdx) / max(wellAreaMap(rowIdx, colIdx), eps);
            if isfinite(semipolarRatio) && semipolarRatio > 0
                shellThicknessScaleMap(rowIdx, colIdx) = clampValue((semipolarRatio / max(semipolarToWellRef, eps)) ^ 0.40, 0.75, 1.35);
            end

            if isfinite(totalIntensityMap(rowIdx, colIdx)) && totalIntensityMap(rowIdx, colIdx) > 0
                totalIntensityGainMap(rowIdx, colIdx) = clampValue( ...
                    (totalIntensityMap(rowIdx, colIdx) / max(totalIntensityRef, eps)) ^ 0.72, 0.60, 2.20);
            end

            centerFillGainMap(rowIdx, colIdx) = clampValue(localIntensityContext.centerFillGainMap(rowIdx, colIdx), 0.85, 1.35);

            nGaNGainMap(rowIdx, colIdx) = clampValue(nGaNAreaMap(rowIdx, colIdx) / max(nGaNRef, eps), 0.55, 1.45);
            prestrainedGainMap(rowIdx, colIdx) = clampValue(prestrainedAreaMap(rowIdx, colIdx) / max(prestrainedRef, eps), 0.50, 1.50);
            barrierGainMap(rowIdx, colIdx) = clampValue(barrierAreaMap(rowIdx, colIdx) / max(barrierRef, eps), 0.55, 1.45);
            wellGainMap(rowIdx, colIdx) = clampValue(wellAreaMap(rowIdx, colIdx) / max(wellRef, eps), 0.60, 1.55);
            semipolarGainMap(rowIdx, colIdx) = clampValue(semipolarAreaMap(rowIdx, colIdx) / max(semipolarRef, eps), 0.50, 1.60);
            pTypeGainMap(rowIdx, colIdx) = clampValue(pTypeAreaMap(rowIdx, colIdx) / max(pTypeRef, eps), 0.45, 1.50);
            vpitCollectionGainMap(rowIdx, colIdx) = clampValue( ...
                1 + 0.35 * max(0, shellThicknessScaleMap(rowIdx, colIdx) - 1) + ...
                0.22 * max(0, semipolarGainMap(rowIdx, colIdx) - 1) + ...
                0.18 * max(0, totalIntensityGainMap(rowIdx, colIdx) - 1) + ...
                0.14 * max(0, centerFillGainMap(rowIdx, colIdx) - 1), ...
                0.90, 1.55);
        end
    end

    localField = struct( ...
        'version', 1, ...
        'createdAt', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
        'rows', rows, ...
        'cols', cols, ...
        'nGaNAreaMap', nGaNAreaMap, ...
        'prestrainedAreaMap', prestrainedAreaMap, ...
        'barrierAreaMap', barrierAreaMap, ...
        'wellAreaMap', wellAreaMap, ...
        'semipolarAreaMap', semipolarAreaMap, ...
        'pTypeAreaMap', pTypeAreaMap, ...
        'shortwaveAreaMap', shortwaveAreaMap, ...
        'shoulderAreaMap', shoulderAreaMap, ...
        'rawShoulderAreaMap', rawShoulderAreaMap, ...
        'shortMainRatioMap', shortMainRatioMap, ...
        'shoulderMainRatioMap', shoulderMainRatioMap, ...
        'pTypeShareMap', pTypeShareMap, ...
        'vpitInfluenceMap', vpitInfluenceMap, ...
        'vpitCoreMask', vpitCoreMask, ...
        'vpitRimMask', vpitRimMask, ...
        'regionCategoryMap', regionCategoryMap, ...
        'totalIntensityMap', totalIntensityMap, ...
        'nGaNPeakMap_nm', nGaNPeakMap_nm, ...
        'gaNPeakMap_nm', nGaNPeakMap_nm, ...
        'prestrainedPeakMap_nm', prestrainedPeakMap_nm, ...
        'barrierPeakMap_nm', barrierPeakMap_nm, ...
        'wellPeakMap_nm', wellPeakMap_nm, ...
        'semipolarPeakMap_nm', semipolarPeakMap_nm, ...
        'pTypePeakMap_nm', pTypePeakMap_nm, ...
        'wellInMap', wellInMap, ...
        'barrierInMap', barrierInMap, ...
        'semipolarInMap', semipolarInMap, ...
        'shellThicknessScaleMap', shellThicknessScaleMap, ...
        'totalIntensityGainMap', totalIntensityGainMap, ...
        'vpitCollectionGainMap', vpitCollectionGainMap, ...
        'centerFillGainMap', centerFillGainMap, ...
        'nGaNGainMap', nGaNGainMap, ...
        'gaNGainMap', nGaNGainMap, ...
        'prestrainedGainMap', prestrainedGainMap, ...
        'barrierGainMap', barrierGainMap, ...
        'wellGainMap', wellGainMap, ...
        'semipolarGainMap', semipolarGainMap, ...
        'pTypeGainMap', pTypeGainMap, ...
        'globalNGaNPeakNm', getFieldOrDefault(profile, 'shortwavePeakNm', getFieldOrDefault(profile, 'nGaNPeakNm', getFieldOrDefault(profile, 'gaNPeakNm', 390))), ...
        'globalPrestrainedPeakNm', getFieldOrDefault(profile, 'shoulderPeakNm', getFieldOrDefault(profile, 'prestrainedPeakNm', 413)), ...
        'globalWellIn', globalWellIn, ...
        'globalBarrierIn', globalBarrierIn, ...
        'globalSemipolarIn', globalSemipolarIn, ...
        'globalPTypePeakNm', getFieldOrDefault(profile, 'pTypePeakNm', 388));
end

function intensityContext = computeLocalIntensityContext(totalIntensityMap)
    rows = size(totalIntensityMap, 1);
    cols = size(totalIntensityMap, 2);
    centerFillGainMap = ones(rows, cols);
    ringNeighborhoodMean = nan(rows, cols);

    [gridX, gridY] = meshgrid(-2:2, -2:2);
    centerMask = hypot(gridX, gridY) <= 1.1;
    ringMask = hypot(gridX, gridY) > 1.1 & hypot(gridX, gridY) <= 2.5;

    for rowIdx = 1:rows
        for colIdx = 1:cols
            if ~isfinite(totalIntensityMap(rowIdx, colIdx)) || totalIntensityMap(rowIdx, colIdx) <= 0
                continue;
            end

            rowRange = max(1, rowIdx - 2):min(rows, rowIdx + 2);
            colRange = max(1, colIdx - 2):min(cols, colIdx + 2);
            patch = totalIntensityMap(rowRange, colRange);

            localCenterMask = centerMask((rowRange - rowIdx) + 3, (colRange - colIdx) + 3);
            localRingMask = ringMask((rowRange - rowIdx) + 3, (colRange - colIdx) + 3);
            ringValues = patch(localRingMask & isfinite(patch) & patch > 0);
            centerValues = patch(localCenterMask & isfinite(patch) & patch > 0);
            if isempty(ringValues) || isempty(centerValues)
                continue;
            end

            ringNeighborhoodMean(rowIdx, colIdx) = mean(ringValues);
            centerFillRatio = mean(centerValues) / max(mean(ringValues), eps);
            centerFillGainMap(rowIdx, colIdx) = clampValue(centerFillRatio ^ 0.35, 0.85, 1.35);
        end
    end

    intensityContext = struct( ...
        'centerFillGainMap', centerFillGainMap, ...
        'ringNeighborhoodMean', ringNeighborhoodMean);
end

function [vpitInfluenceMap, vpitCoreMask, vpitRimMask, regionCategoryMap] = ...
        buildVpitRegionMaps(experimentalData, rows, cols, shortwaveAreaMap, shoulderAreaMap, wellAreaMap)
    vpitInfluenceMap = zeros(rows, cols);
    vpitCoreMask = false(rows, cols);
    vpitRimMask = false(rows, cols);

    hasGeometry = isstruct(experimentalData) && isfield(experimentalData, 'vPitCentroids_nm') && ...
        ~isempty(experimentalData.vPitCentroids_nm) && isfield(experimentalData, 'vPitRadii_nm') && ...
        ~isempty(experimentalData.vPitRadii_nm) && isfield(experimentalData, 'scanX_nm') && ...
        isfield(experimentalData, 'scanY_nm');
    if hasGeometry
        centers_nm = double(experimentalData.vPitCentroids_nm);
        radii_nm = double(experimentalData.vPitRadii_nm(:));
        if size(centers_nm, 2) ~= 2
            centers_nm = zeros(0, 2);
            radii_nm = zeros(0, 1);
        else
            entryCount = min(size(centers_nm, 1), numel(radii_nm));
            centers_nm = centers_nm(1:entryCount, :);
            radii_nm = radii_nm(1:entryCount);
            validMask = all(isfinite(centers_nm), 2) & isfinite(radii_nm) & radii_nm > 0;
            centers_nm = centers_nm(validMask, :);
            radii_nm = radii_nm(validMask);
        end

        if ~isempty(centers_nm)
            scanX_nm = double(experimentalData.scanX_nm);
            scanY_nm = double(experimentalData.scanY_nm);
            pixelSizeX_nm = scanX_nm / max(cols, 1);
            pixelSizeY_nm = scanY_nm / max(rows, 1);
            [colGrid, rowGrid] = meshgrid(1:cols, 1:rows);
            xGrid_nm = colGrid * pixelSizeX_nm - scanX_nm / 2;
            yGrid_nm = rowGrid * pixelSizeY_nm - scanY_nm / 2;

            for pitIdx = 1:size(centers_nm, 1)
                distanceNorm = hypot(xGrid_nm - centers_nm(pitIdx, 1), yGrid_nm - centers_nm(pitIdx, 2)) ./ ...
                    max(radii_nm(pitIdx), eps);
                pitInfluence = exp(-0.5 * (distanceNorm / 1.05) .^ 2);
                vpitInfluenceMap = max(vpitInfluenceMap, pitInfluence);
                vpitCoreMask = vpitCoreMask | distanceNorm <= 0.72;
                vpitRimMask = vpitRimMask | (distanceNorm > 0.72 & distanceNorm <= 1.45);
            end
        end
    end

    spectralVpitMap = estimateSpectralVpitInfluence(shortwaveAreaMap, shoulderAreaMap, wellAreaMap);
    vpitInfluenceMap = max(vpitInfluenceMap, spectralVpitMap);
    vpitInfluenceMap = clampValue(vpitInfluenceMap, 0, 1);
    if ~any(vpitCoreMask(:))
        vpitCoreMask = vpitInfluenceMap >= 0.72;
    end
    if ~any(vpitRimMask(:))
        vpitRimMask = vpitInfluenceMap >= 0.38 & vpitInfluenceMap < 0.72;
    end

    regionCategoryMap = ones(rows, cols);
    regionCategoryMap(vpitRimMask) = 2;
    regionCategoryMap(vpitCoreMask) = 3;
end

function influenceMap = estimateSpectralVpitInfluence(shortwaveAreaMap, shoulderAreaMap, wellAreaMap)
    shortRatio = shortwaveAreaMap ./ max(wellAreaMap, eps);
    shoulderRatio = shoulderAreaMap ./ max(wellAreaMap, eps);
    shortNorm = normalizeFiniteMap(shortRatio);
    shoulderNorm = normalizeFiniteMap(shoulderRatio);
    influenceMap = clampValue(0.78 * shortNorm + 0.22 * shoulderNorm, 0, 1);
end

function normalizedMap = normalizeFiniteMap(values)
    normalizedMap = zeros(size(values));
    validValues = values(isfinite(values));
    if isempty(validValues)
        return;
    end
    lowValue = prctile(validValues, 12);
    highValue = prctile(validValues, 92);
    if highValue <= lowValue
        highValue = max(validValues);
        lowValue = min(validValues);
    end
    if highValue <= lowValue
        return;
    end
    normalizedMap = (values - lowValue) ./ (highValue - lowValue);
    normalizedMap(~isfinite(normalizedMap)) = 0;
    normalizedMap = clampValue(normalizedMap, 0, 1);
end

function bandRange = buildPeakBand(centerNm, halfWidthNm, limits)
    bandRange = [centerNm - halfWidthNm, centerNm + halfWidthNm];
    bandRange(1) = max(bandRange(1), limits(1));
    bandRange(2) = min(bandRange(2), limits(2));
end

function gain = computeBandRatioGain(axisNm, experimentalSpectrum, simulationSpectrum, bandRange, minGain, maxGain)
    bandMask = axisNm >= bandRange(1) & axisNm <= bandRange(2);
    if ~any(bandMask)
        gain = 1;
        return;
    end

    expBand = trapz(axisNm(bandMask), experimentalSpectrum(bandMask));
    simBand = trapz(axisNm(bandMask), simulationSpectrum(bandMask));
    if simBand <= eps || expBand <= 0
        gain = 1;
        return;
    end

    gain = expBand / simBand;
    gain = clampValue(gain, minGain, maxGain);
end

function area = integrateBandArea(axisNm, normalizedSpectrum, bandRange)
    mask = axisNm >= bandRange(1) & axisNm <= bandRange(2);
    if ~any(mask)
        area = 0;
        return;
    end
    area = trapz(axisNm(mask), normalizedSpectrum(mask));
end

function centroidNm = computeBandCentroid(axisNm, normalizedSpectrum, bandRange, defaultCenterNm)
    mask = axisNm >= bandRange(1) & axisNm <= bandRange(2);
    if ~any(mask)
        centroidNm = defaultCenterNm;
        return;
    end
    weights = normalizedSpectrum(mask) .^ 1.3;
    if sum(weights) <= eps
        centroidNm = defaultCenterNm;
        return;
    end
    centroidNm = sum(axisNm(mask) .* weights) / sum(weights);
end

function value = medianPositive(values, defaultValue)
    values = values(isfinite(values) & values > 0);
    if isempty(values)
        value = defaultValue;
    else
        value = median(values);
    end
end

function value = lowerQuantile(values, fraction)
    values = sort(double(values(:)));
    values = values(isfinite(values));
    if isempty(values)
        value = 0;
        return;
    end
    index = max(1, min(numel(values), round(1 + fraction * (numel(values) - 1))));
    value = values(index);
end


function errorValue = computeShapeError(referenceSpectrum, testSpectrum)
    difference = referenceSpectrum(:) - testSpectrum(:);
    errorValue = sqrt(mean(difference .^ 2));
end

function errorValue = computeWeightedFitError(referenceSpectrum, testSpectrum, axisNm)
    referenceSpectrum = normalizeSpectrum(referenceSpectrum);
    testSpectrum = normalizeSpectrum(testSpectrum);
    weights = 0.18 + 1.25 * sqrt(max(referenceSpectrum, 0));
    if nargin >= 3 && ~isempty(axisNm)
        weights = applyPhysicsWindowWeights(axisNm(:), referenceSpectrum, weights);
    end
    weights = weights / mean(weights);
    difference = referenceSpectrum(:) - testSpectrum(:);
    errorValue = sqrt(mean(weights(:) .* (difference .^ 2)));
end

function weights = applyPhysicsWindowWeights(axisNm, referenceSpectrum, weights)
    weightedWindows = [ ...
        376, 396, 0.010, 5.5; ...
        398, 432, 0.012, 3.8; ...
        442, 475, 0.012, 2.5; ...
        492, 545, 0.018, 3.2];
    for windowIdx = 1:size(weightedWindows, 1)
        windowRange = weightedWindows(windowIdx, 1:2);
        minSignal = weightedWindows(windowIdx, 3);
        gain = weightedWindows(windowIdx, 4);
        mask = axisNm >= windowRange(1) & axisNm <= windowRange(2);
        if ~any(mask) || max(referenceSpectrum(mask)) < minSignal
            continue;
        end
        localShape = normalizeSpectrum(referenceSpectrum(mask));
        weights(mask) = weights(mask) .* (1 + gain * (0.35 + localShape));
    end
end

function penalty = computePhysicsWindowRatioPenalty(axisNm, referenceSpectrum, testSpectrum)
    referenceSpectrum = normalizeSpectrum(referenceSpectrum);
    testSpectrum = normalizeSpectrum(testSpectrum);
    mainBand = [492, 545];
    mainReferenceArea = max(integrateBandArea(axisNm, referenceSpectrum, mainBand), eps);
    mainTestArea = max(integrateBandArea(axisNm, testSpectrum, mainBand), eps);
    shortBands = [376, 396; 398, 432; 442, 475];
    penalty = 0;
    for bandIdx = 1:size(shortBands, 1)
        referenceArea = integrateBandArea(axisNm, referenceSpectrum, shortBands(bandIdx, :));
        testArea = integrateBandArea(axisNm, testSpectrum, shortBands(bandIdx, :));
        referenceRatio = referenceArea / mainReferenceArea;
        testRatio = testArea / mainTestArea;
        if referenceRatio < 0.006
            continue;
        end
        logReference = log1p(28 * referenceRatio);
        logTest = log1p(28 * testRatio);
        bandWeight = 2.5 + 1.5 * (bandIdx == 1);
        penalty = penalty + bandWeight * (logReference - logTest) .^ 2;
    end
end

function penalty = computeMainWindowShapePenalty(axisNm, referenceSpectrum, testSpectrum)
    referenceSpectrum = normalizeSpectrum(referenceSpectrum);
    testSpectrum = normalizeSpectrum(testSpectrum);
    mainBand = [492, 545];
    mainMask = axisNm >= mainBand(1) & axisNm <= mainBand(2);
    penalty = 0;
    if ~any(mainMask) || max(referenceSpectrum(mainMask)) < 0.018
        return;
    end

    mainAxis = axisNm(mainMask);
    referenceMain = referenceSpectrum(mainMask);
    testMain = testSpectrum(mainMask);
    [~, refPeakIdx] = max(referenceMain);
    [~, testPeakIdx] = max(testMain);
    referencePeakNm = mainAxis(refPeakIdx);
    testPeakNm = mainAxis(testPeakIdx);
    peakShiftNm = testPeakNm - referencePeakNm;
    penalty = penalty + 0.00018 * peakShiftNm .^ 2;

    referenceMainArea = max(integrateBandArea(axisNm, referenceSpectrum, mainBand), eps);
    testMainArea = max(integrateBandArea(axisNm, testSpectrum, mainBand), eps);
    redTailBand = [min(max(referencePeakNm + 16, 522), 540), 570];
    referenceRedRatio = integrateBandArea(axisNm, referenceSpectrum, redTailBand) / referenceMainArea;
    testRedRatio = integrateBandArea(axisNm, testSpectrum, redTailBand) / testMainArea;
    allowedRedRatio = referenceRedRatio * 1.10 + 0.012;
    redTailExcess = max(testRedRatio - allowedRedRatio, 0);
    if redTailExcess > 0
        penalty = penalty + 2.4 * log1p(6.0 * redTailExcess) .^ 2;
    end
end

function metrics = computeCalibrationWindowMetrics(axisNm, referenceSpectrum, fittedSpectrum)
    referenceSpectrum = normalizeSpectrum(referenceSpectrum);
    fittedSpectrum = normalizeSpectrum(fittedSpectrum);
    shortwaveBand = [376, 396];
    shoulderBand = [398, 432];
    mainBand = [492, 545];
    shortwaveMask = axisNm >= shortwaveBand(1) & axisNm <= shortwaveBand(2);
    if any(shortwaveMask)
        [~, localIdx] = max(referenceSpectrum(shortwaveMask));
        shortAxis = axisNm(shortwaveMask);
        shortwavePeakNm = shortAxis(localIdx);
        shortwaveFitError = sqrt(mean((referenceSpectrum(shortwaveMask) - fittedSpectrum(shortwaveMask)) .^ 2));
    else
        shortwavePeakNm = NaN;
        shortwaveFitError = NaN;
    end
    shoulderMask = axisNm >= shoulderBand(1) & axisNm <= shoulderBand(2);
    if any(shoulderMask)
        [~, localIdx] = max(referenceSpectrum(shoulderMask));
        shoulderAxis = axisNm(shoulderMask);
        shoulderPeakNm = shoulderAxis(localIdx);
        shoulderFitError = sqrt(mean((referenceSpectrum(shoulderMask) - fittedSpectrum(shoulderMask)) .^ 2));
    else
        shoulderPeakNm = NaN;
        shoulderFitError = NaN;
    end
    mainMask = axisNm >= mainBand(1) & axisNm <= mainBand(2);
    if any(mainMask)
        [~, localIdx] = max(referenceSpectrum(mainMask));
        mainAxis = axisNm(mainMask);
        mainPeakNm = mainAxis(localIdx);
    else
        mainPeakNm = NaN;
    end
    shortwaveArea = integrateBandArea(axisNm, referenceSpectrum, shortwaveBand);
    shoulderArea = integrateBandArea(axisNm, referenceSpectrum, shoulderBand);
    mainPeakArea = integrateBandArea(axisNm, referenceSpectrum, mainBand);
    fittedShortwaveArea = integrateBandArea(axisNm, fittedSpectrum, shortwaveBand);
    fittedShoulderArea = integrateBandArea(axisNm, fittedSpectrum, shoulderBand);
    fittedMainPeakArea = integrateBandArea(axisNm, fittedSpectrum, mainBand);
    metrics = struct( ...
        'shortwaveBandNm', shortwaveBand, ...
        'shoulderBandNm', shoulderBand, ...
        'mainBandNm', mainBand, ...
        'shortwavePeakNm', shortwavePeakNm, ...
        'shoulderPeakNm', shoulderPeakNm, ...
        'mainPeakNm', mainPeakNm, ...
        'shortwaveFitError', shortwaveFitError, ...
        'shoulderFitError', shoulderFitError, ...
        'shortwaveArea', shortwaveArea, ...
        'shoulderArea', shoulderArea, ...
        'mainPeakArea', mainPeakArea, ...
        'fittedShortwaveArea', fittedShortwaveArea, ...
        'fittedShoulderArea', fittedShoulderArea, ...
        'fittedMainPeakArea', fittedMainPeakArea, ...
        'shortwaveToMainAreaRatio', shortwaveArea / max(mainPeakArea, eps), ...
        'shoulderToMainAreaRatio', shoulderArea / max(mainPeakArea, eps), ...
        'fittedShortwaveToMainAreaRatio', fittedShortwaveArea / max(fittedMainPeakArea, eps));
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

function sigma = fwhmToSigma(fwhm)
    sigma = fwhm / (2 * sqrt(2 * log(2)));
end

function inComposition = invertInComposition(peakNm, defaultValue)
    EgGaN = 3.3032;
    EgInN = 0.6086;
    targetEnergy = 1240 / peakNm;
    coefficients = [1.43, EgInN - EgGaN - 1.43, EgGaN - targetEnergy];
    rootsValue = roots(coefficients);
    rootsValue = rootsValue(imag(rootsValue) == 0);
    rootsValue = real(rootsValue);
    rootsValue = rootsValue(rootsValue >= 0 & rootsValue <= 1);
    if isempty(rootsValue)
        inComposition = defaultValue;
        return;
    end
    [~, bestIdx] = min(abs(rootsValue - defaultValue));
    inComposition = rootsValue(bestIdx);
end


function componentsNm = convertSpectrumModelToComponentsNm(spectrumModel)
    if isstruct(spectrumModel) && isfield(spectrumModel, 'components') && ~isempty(spectrumModel.components)
        rawComponents = double(spectrumModel.components);
        componentsNm = [rawComponents(:,1), rawComponents(:,2) * 1e9, rawComponents(:,3) * 1e9];
        if size(rawComponents, 2) >= 4
            componentsNm(:, 4) = rawComponents(:, 4);
        end
    elseif isnumeric(spectrumModel) && numel(spectrumModel) >= 6
        rawValues = double(spectrumModel(:)');
        componentCount = floor(numel(rawValues) / 3);
        rawValues = rawValues(1:3*componentCount);
        componentsNm = reshape(rawValues, 3, [])';
        componentsNm(:,2:3) = componentsNm(:,2:3) * 1e9;
    else
        componentsNm = [1, 500, 40];
    end

    componentsNm(:,4) = ensureEtaColumnNm(componentsNm);
    componentsNm(:,1) = normalizePositive(componentsNm(:,1));
    [~, order] = sort(componentsNm(:,2), 'ascend');
    componentsNm = componentsNm(order, :);
end

function legacySpectrum = componentsNmToLegacySpectrum(componentsNm)
    componentsNm = componentsNm(1:min(2, size(componentsNm, 1)), :);
    if size(componentsNm, 1) == 1
        componentsNm = [componentsNm; componentsNm];
    end
    legacySpectrum = [ ...
        max(1, 100 * componentsNm(1,1)), componentsNm(1,2) * 1e-9, componentsNm(1,3) * 1e-9, ...
        max(1, 100 * componentsNm(2,1)), componentsNm(2,2) * 1e-9, componentsNm(2,3) * 1e-9];
end

function etaColumnNm = ensureEtaColumnNm(componentsNm)
    if size(componentsNm, 2) >= 4
        etaColumnNm = clampValue(double(componentsNm(:, 4)), 0.05, 0.95);
    else
        etaColumnNm = 0.72 * ones(size(componentsNm, 1), 1);
    end
end

function eta = getEtaFromComponent(component)
    if numel(component) >= 4 && isfinite(component(4))
        eta = clampValue(component(4), 0.05, 0.95);
    else
        eta = 0.72;
    end
end

function value = getLayerSpectrum(layerParams, layerName)
    rowIdx = find(strcmp(layerParams(:, 1), layerName), 1);
    value = layerParams{rowIdx, 8};
end

function value = getLayerComposition(layerParams, layerName)
    rowIdx = find(strcmp(layerParams(:, 1), layerName), 1);
    value = layerParams{rowIdx, 2};
end

function value = getFieldOrDefault(data, fieldName, defaultValue)
    if isstruct(data) && isfield(data, fieldName) && ~isempty(data.(fieldName))
        value = data.(fieldName);
    else
        value = defaultValue;
    end
end

function value = clampValue(value, lowerBound, upperBound)
    value = max(lowerBound, min(upperBound, value));
end

function mapData = loadTIFFForVpitsFile(fullpath)
    if nargin < 1 || isempty(fullpath)
        error('precies:vpit:MissingFile', 'A TIFF file path is required.');
    end
    if ~isfile(fullpath)
        error('precies:vpit:MissingFile', 'TIFF file does not exist: %s', fullpath);
    end

    warnState = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
    cleanupObj = onCleanup(@() warning(warnState));
    tiffInfo = imfinfo(fullpath);
    numLayers = numel(tiffInfo);

    t = Tiff(fullpath, 'r');
    omeXML = '';
    try
        omeXML = t.getTag('ImageDescription');
    catch
    end
    t.close();

    [imageTypes, ~] = precies.tiffops('analyzeTiffStructure', tiffInfo);
    spectralLayerIdx = find(contains(imageTypes, 'Spectrum with Spectrometer', 'IgnoreCase', true), 1);
    if isempty(spectralLayerIdx)
        spectralLayerIdx = find(contains(imageTypes, 'spectral', 'IgnoreCase', true), 1);
    end
    if isempty(spectralLayerIdx)
        error('precies:vpit:NoSpectralLayers', 'No spectral layers found in TIFF file: %s', fullpath);
    end

    spectralHeight = tiffInfo(spectralLayerIdx).Height;
    spectralWidth = tiffInfo(spectralLayerIdx).Width;
    [pixelSizeX_um, pixelSizeY_um] = readPixelSizeFromTiffInfo(tiffInfo(spectralLayerIdx));
    scanX_nm = spectralWidth * pixelSizeX_um * 1000;
    scanY_nm = spectralHeight * pixelSizeY_um * 1000;

    surveyLayerIdx = find(contains(imageTypes, 'survey', 'IgnoreCase', true), 1);
    if isempty(surveyLayerIdx)
        surveyLayerIdx = min(2, numLayers);
    end
    surveyHeight = tiffInfo(surveyLayerIdx).Height;
    surveyWidth = tiffInfo(surveyLayerIdx).Width;

    wavelengthLayerCount = numLayers - spectralLayerIdx + 1;
    dataCube = zeros(surveyHeight, surveyWidth, wavelengthLayerCount, 'uint16');
    problematicLayers = zeros(wavelengthLayerCount, 1);
    problematicLayerCount = 0;
    for layerIdx = spectralLayerIdx:numLayers
        layerNum = layerIdx - spectralLayerIdx + 1;
        try
            wavelengthData = imread(fullpath, layerIdx);
            if size(wavelengthData, 1) ~= surveyHeight || size(wavelengthData, 2) ~= surveyWidth
                wavelengthData = imresize(wavelengthData, [surveyHeight, surveyWidth]);
            end
            dataCube(:, :, layerNum) = wavelengthData;
        catch
            problematicLayerCount = problematicLayerCount + 1;
            problematicLayers(problematicLayerCount) = layerIdx;
        end
    end
    problematicLayers = problematicLayers(1:problematicLayerCount);

    wavelengthAxis = precies.tiffops('parseOmeXmlForWavelength', omeXML);
    if isempty(wavelengthAxis) || numel(wavelengthAxis) ~= wavelengthLayerCount
        wavelengthAxis = precies.tiffops('tryAlternativeWavelengthSources', omeXML, tiffInfo, wavelengthLayerCount);
    end
    if isempty(wavelengthAxis) || numel(wavelengthAxis) ~= wavelengthLayerCount
        wavelengthAxis = linspace(400, 800, wavelengthLayerCount)';
    end
    wavelengthAxis = double(wavelengthAxis(:));

    spectraData = cell(spectralHeight, spectralWidth);
    scaleX = spectralWidth / surveyWidth;
    scaleY = spectralHeight / surveyHeight;
    for rowIdx = 1:spectralHeight
        for colIdx = 1:spectralWidth
            surveyRow = round(rowIdx / scaleY);
            surveyCol = round(colIdx / scaleX);
            if surveyRow >= 1 && surveyRow <= surveyHeight && surveyCol >= 1 && surveyCol <= surveyWidth
                intensities = double(squeeze(dataCube(surveyRow, surveyCol, :)));
            else
                intensities = zeros(size(wavelengthAxis));
            end
            spectraData{rowIdx, colIdx} = [wavelengthAxis, intensities(:)];
        end
    end

    detectionBands = deriveVpitsDetectionBands(spectraData);
    [vPitCentroids_pixels, vPitRadii_nm, vPitProperties] = ...
        detectVpitsFromSpectralData(spectraData, detectionBands, pixelSizeX_um, pixelSizeY_um);
    if isempty(vPitCentroids_pixels)
        vPitCentroids_pixels = zeros(0, 2);
        vPitRadii_nm = zeros(0, 1);
    end

    vPitCentroids_nm = vPitCentroids_pixels .* [pixelSizeX_um * 1000, pixelSizeY_um * 1000] - ...
        [scanX_nm / 2, scanY_nm / 2];
    [~, filename, ext] = fileparts(fullpath);
    mapData = struct( ...
        'spectraData', {spectraData}, ...
        'spectralRows', spectralHeight, ...
        'spectralCols', spectralWidth, ...
        'pixelSizeX_um', pixelSizeX_um, ...
        'pixelSizeY_um', pixelSizeY_um, ...
        'scanX_nm', scanX_nm, ...
        'scanY_nm', scanY_nm, ...
        'vPitCentroids_nm', vPitCentroids_nm, ...
        'vPitRadii_nm', vPitRadii_nm, ...
        'vPitOrientations_deg', extractPropertyVector(vPitProperties, 'OrientationDeg', size(vPitCentroids_nm, 1), 0), ...
        'vPitDetectionScores', extractPropertyVector(vPitProperties, 'DetectionScore', size(vPitCentroids_nm, 1), 1), ...
        'wavelengthAxis', wavelengthAxis, ...
        'sourceFile', fullpath, ...
        'sourceName', [filename ext], ...
        'tiffInfo', tiffInfo, ...
        'imageTypes', {imageTypes}, ...
        'detectionBands', detectionBands, ...
        'problematicLayers', problematicLayers);
end

function [pixelSizeX_um, pixelSizeY_um] = readPixelSizeFromTiffInfo(info)
    if isfield(info, 'XResolution') && ~isempty(info.XResolution)
        xResInch = info.XResolution;
        yResInch = info.YResolution;
        if isfield(info, 'ResolutionUnit') && strcmpi(info.ResolutionUnit, 'Centimeter')
            xResInch = xResInch * 2.54;
            yResInch = yResInch * 2.54;
        end
        pixelSizeX_um = 1 / xResInch * 25400;
        pixelSizeY_um = 1 / yResInch * 25400;
    else
        pixelSizeX_um = 0.024132355839844;
        pixelSizeY_um = 0.024132355839844;
    end
end

function [vPitCentroids_pixels, vPitRadii_pixels, vPitProperties] = detectVpitsFromSpectralData(spectraData, detectionBands, pixelSizeX_um, pixelSizeY_um)
    [featureMaps, enhancedMap] = createMultiFeatureMaps(spectraData, detectionBands, pixelSizeX_um, pixelSizeY_um);

    [vPitCentroids_pixels, vPitRadii_pixels, vPitProperties] = improvedVPitDetection( ...
        enhancedMap, featureMaps, pixelSizeX_um, pixelSizeY_um);
    [vPitCentroids_pixels, vPitRadii_pixels, vPitProperties] = postProcessDetections(...
        vPitCentroids_pixels, vPitRadii_pixels, vPitProperties, enhancedMap);
end

function detectionBands = deriveVpitsDetectionBands(spectraData)
    wavelengthAxis = [];
    accumulatedSpectrum = [];
    validCount = 0;
    [rows, cols] = size(spectraData);

    for rowIdx = 1:rows
        for colIdx = 1:cols
            spectrum = spectraData{rowIdx, colIdx};
            if isempty(spectrum) || size(spectrum, 2) < 2
                continue;
            end
            if isempty(wavelengthAxis)
                wavelengthAxis = spectrum(:, 1);
                accumulatedSpectrum = zeros(size(wavelengthAxis));
            end
            intensities = double(spectrum(:, 2));
            if numel(intensities) ~= numel(wavelengthAxis)
                continue;
            end
            accumulatedSpectrum = accumulatedSpectrum + intensities;
            validCount = validCount + 1;
        end
    end

    if isempty(wavelengthAxis) || validCount == 0
        detectionBands = struct( ...
            'shortBand_nm', [372, 418], ...
            'longBand_nm', [490, 550], ...
            'contextBand_nm', [440, 590]);
        return;
    end

    meanSpectrum = accumulatedSpectrum / validCount;
    shortMask = wavelengthAxis >= 340 & wavelengthAxis <= 450;
    longMask = wavelengthAxis >= 460 & wavelengthAxis <= 620;

    if any(shortMask)
        [~, shortLocalIdx] = max(meanSpectrum(shortMask));
        shortAxis = wavelengthAxis(shortMask);
        shortPeak_nm = shortAxis(shortLocalIdx);
    else
        shortPeak_nm = 392;
    end

    if any(longMask)
        [~, longLocalIdx] = max(meanSpectrum(longMask));
        longAxis = wavelengthAxis(longMask);
        longPeak_nm = longAxis(longLocalIdx);
    else
        longPeak_nm = 520;
    end

    detectionBands = struct( ...
        'shortBand_nm', [max(330, shortPeak_nm - 18), min(460, shortPeak_nm + 20)], ...
        'longBand_nm', [max(450, longPeak_nm - 28), min(650, longPeak_nm + 35)], ...
        'contextBand_nm', [max(420, longPeak_nm - 50), min(650, longPeak_nm + 60)]);
end

function [featureMaps, enhancedMap] = createMultiFeatureMaps(spectraData, detectionBands, pixelSizeX_um, pixelSizeY_um)
    [rows, cols] = size(spectraData);
    shortMap = zeros(rows, cols);
    longMap = zeros(rows, cols);
    totalMap = zeros(rows, cols);

    for i = 1:rows
        for j = 1:cols
            spectrum = spectraData{i, j};
            if isempty(spectrum)
                continue;
            end
            shortMask = spectrum(:, 1) >= detectionBands.shortBand_nm(1) & spectrum(:, 1) <= detectionBands.shortBand_nm(2);
            longMask = spectrum(:, 1) >= detectionBands.longBand_nm(1) & spectrum(:, 1) <= detectionBands.longBand_nm(2);
            contextMask = spectrum(:, 1) >= detectionBands.contextBand_nm(1) & spectrum(:, 1) <= detectionBands.contextBand_nm(2);

            shortMap(i, j) = sum(spectrum(shortMask, 2));
            longMap(i, j) = sum(spectrum(longMask, 2));
            totalMap(i, j) = sum(spectrum(contextMask, 2));
        end
    end

    avgPixelSize_nm = ((pixelSizeX_um + pixelSizeY_um) / 2) * 1000;
    minRadiusPx = max(2, round(35 / avgPixelSize_nm));
    shortNorm = mat2gray(log1p(shortMap));
    totalNorm = mat2gray(log1p(totalMap));
    ratioMap = shortMap ./ max(longMap, 1);
    ratioNorm = mat2gray(log1p(ratioMap));
    fractionMap = shortMap ./ max(totalMap, 1);
    fractionNorm = mat2gray(fractionMap);
    localBaseline = imgaussfilt(fractionNorm, max(1.2, 1.4 * minRadiusPx));
    localContrast = max(0, fractionNorm - localBaseline);
    [gx, gy] = gradient(imgaussfilt(fractionNorm, 0.8));
    gradientNorm = mat2gray(sqrt(gx.^2 + gy.^2));
    dogSmall = imgaussfilt(fractionNorm, max(0.8, 0.45 * minRadiusPx));
    dogLarge = imgaussfilt(fractionNorm, max(1.5, 0.85 * minRadiusPx));
    dogNorm = mat2gray(max(0, dogSmall - dogLarge));

    enhancedMap = 0.34 * fractionNorm + ...
                  0.26 * ratioNorm + ...
                  0.18 * localContrast + ...
                  0.12 * gradientNorm + ...
                  0.10 * dogNorm + ...
                  0.06 * totalNorm + ...
                  0.04 * shortNorm;
    enhancedMap = imadjust(mat2gray(enhancedMap), stretchlim(enhancedMap, [0.02 0.995]), [0 1]);
    edgeMap = mat2gray(imgradient(imgaussfilt(totalNorm + 1.4 * fractionNorm, 0.9)));

    featureMaps = struct( ...
        'shortMap', shortMap, ...
        'longMap', longMap, ...
        'totalMap', totalMap, ...
        'ratioMap', ratioMap, ...
        'fractionMap', fractionMap, ...
        'edgeMap', edgeMap, ...
        'enhancedMap', enhancedMap);
end

function [centroids, radii, properties] = improvedVPitDetection(enhancedMap, featureMaps, pixelSizeX_um, pixelSizeY_um)
    avgPixelSize_nm = ((pixelSizeX_um + pixelSizeY_um) / 2) * 1000;
    minRadiusPx = max(2, round(35 / avgPixelSize_nm));
    maxRadiusPx = max(minRadiusPx + 2, round(180 / avgPixelSize_nm));
    circleSearchMinPx = max(6, minRadiusPx);
    circleSearchMaxPx = max(circleSearchMinPx + 2, maxRadiusPx);

    candidateMap = imgaussfilt(enhancedMap, 0.9);
    ringResponseMap = mat2gray(imgradient(candidateMap));

    [centersA, radiiA, metricA] = imfindcircles(candidateMap, [circleSearchMinPx, circleSearchMaxPx], ...
        'ObjectPolarity', 'bright', 'Sensitivity', 0.93, 'EdgeThreshold', 0.03);
    [centersB, radiiB, metricB] = imfindcircles(ringResponseMap, [circleSearchMinPx, circleSearchMaxPx], ...
        'ObjectPolarity', 'bright', 'Sensitivity', 0.90, 'EdgeThreshold', 0.02);
    [centersC, radiiC, metricC] = proposeSplitCandidates(candidateMap, minRadiusPx, maxRadiusPx);

    [centersA, radiiA, metricA] = normalizeCircleCandidates(centersA, radiiA, metricA);
    [centersB, radiiB, metricB] = normalizeCircleCandidates(centersB, radiiB, metricB);
    [centersC, radiiC, metricC] = normalizeCircleCandidates(centersC, radiiC, metricC);

    candidateCenters = [centersA; centersB; centersC];
    candidateRadiiPx = [radiiA; radiiB; radiiC];
    candidateScores = [metricA; 0.86 * metricB; 0.84 * metricC];

    if isempty(candidateCenters)
        centroids = zeros(0, 2);
        radii = zeros(0, 1);
        properties = struct([]);
        return;
    end

    centroids = zeros(0, 2);
    radii = zeros(0, 1);
    properties = struct([]);
    for idx = 1:size(candidateCenters, 1)
        [isValid, prop] = evaluateVpitCandidate(candidateCenters(idx, :), candidateRadiiPx(idx), ...
            candidateScores(idx), featureMaps, avgPixelSize_nm);
        if ~isValid
            continue;
        end

        centroids(end + 1, :) = prop.Centroid;
        radii(end + 1, 1) = prop.Radius_nm;
        properties = [properties; prop];
    end
end

function [centers, radiiPx, metrics] = normalizeCircleCandidates(centers, radiiPx, metrics)
    if isempty(centers)
        centers = zeros(0, 2);
    elseif size(centers, 2) ~= 2
        centers = reshape(centers, [], 2);
    end

    radiiPx = radiiPx(:);
    metrics = metrics(:);
    count = min([size(centers, 1), numel(radiiPx), numel(metrics)]);

    if count < 1
        centers = zeros(0, 2);
        radiiPx = zeros(0, 1);
        metrics = zeros(0, 1);
        return;
    end

    centers = centers(1:count, :);
    radiiPx = radiiPx(1:count);
    metrics = metrics(1:count);
end

function [centers, radiiPx, metrics] = proposeSplitCandidates(candidateMap, minRadiusPx, maxRadiusPx)
    smoothedMap = imgaussfilt(candidateMap, 0.8);
    localMaxMask = imregionalmax(smoothedMap);
    finiteValues = smoothedMap(isfinite(smoothedMap));
    if isempty(finiteValues)
        finiteValues = 0;
    end
    intensityThreshold = max(prctile(finiteValues, 90), mean(finiteValues, 'omitnan') + 0.45 * std(finiteValues, 0, 'omitnan'));
    localMaxMask = localMaxMask & smoothedMap >= intensityThreshold;
    [rowIdx, colIdx] = find(localMaxMask);

    if isempty(rowIdx)
        centers = zeros(0, 2);
        radiiPx = zeros(0, 1);
        metrics = zeros(0, 1);
        return;
    end

    candidateStrength = smoothedMap(localMaxMask);
    candidateStrength = candidateStrength(:);
    candidateStrength = candidateStrength / max(candidateStrength, eps);
    baseRadius = 0.55 * minRadiusPx + 0.25 * maxRadiusPx;
    radiiPx = min(maxRadiusPx, max(minRadiusPx, baseRadius * (0.82 + 0.36 * candidateStrength)));
    centers = [colIdx, rowIdx];
    metrics = candidateStrength;
end

function [isValid, prop] = evaluateVpitCandidate(centerPx, radiusPx, detectorScore, featureMaps, avgPixelSize_nm)
    rows = size(featureMaps.enhancedMap, 1);
    cols = size(featureMaps.enhancedMap, 2);
    [gridX, gridY] = meshgrid(1:cols, 1:rows);
    distMap = hypot(gridX - centerPx(1), gridY - centerPx(2));

    centerMask = distMap <= max(1.2, 0.55 * radiusPx);
    ringMask = distMap > 0.70 * radiusPx & distMap <= 1.25 * radiusPx;
    outerMask = distMap > 1.35 * radiusPx & distMap <= 2.05 * radiusPx;

    centerEnhanced = mean(featureMaps.enhancedMap(centerMask), 'omitnan');
    ringEnhanced = mean(featureMaps.enhancedMap(ringMask), 'omitnan');
    outerEnhanced = mean(featureMaps.enhancedMap(outerMask), 'omitnan');
    centerFraction = mean(featureMaps.fractionMap(centerMask), 'omitnan');
    ringFraction = mean(featureMaps.fractionMap(ringMask), 'omitnan');
    outerFraction = mean(featureMaps.fractionMap(outerMask), 'omitnan');
    centerRatio = mean(featureMaps.ratioMap(centerMask), 'omitnan');
    ringRatio = mean(featureMaps.ratioMap(ringMask), 'omitnan');
    outerRatio = mean(featureMaps.ratioMap(outerMask), 'omitnan');
    localTotal = mean(featureMaps.totalMap(distMap <= 1.4 * radiusPx), 'omitnan');
    edgeClearancePx = min([centerPx(1) - 1, centerPx(2) - 1, cols - centerPx(1), rows - centerPx(2)]);

    pitCoreSignal = max(centerEnhanced, ringEnhanced);
    fractionalContrast = max(centerFraction, ringFraction) - outerFraction;
    ratioContrast = max(centerRatio, ringRatio) - outerRatio;
    edgePenalty = min(1, max(0.52, edgeClearancePx / max(radiusPx, eps)));
    qualityScore = edgePenalty * (0.55 * detectorScore + ...
                   0.25 * max(0, pitCoreSignal - outerEnhanced) + ...
                   0.12 * max(0, fractionalContrast) + ...
                   0.08 * max(0, ratioContrast / max(outerRatio, 1e-6)));
    [effectiveRadiusPx, radiusInfo] = estimateVpitOuterRadiusPx(centerPx, radiusPx, featureMaps, distMap);
    edgeIsUsable = edgeClearancePx >= 0.20 * effectiveRadiusPx || qualityScore > 0.32;
    orientationDeg = estimateHexagonOrientationDeg(featureMaps, centerPx, effectiveRadiusPx);

    isValid = isfinite(localTotal) && localTotal > 0 && ...
        pitCoreSignal > outerEnhanced + 0.01 && ...
        fractionalContrast > 0.005 && ...
        qualityScore > 0.12 && edgeIsUsable;

    prop = struct( ...
        'Centroid', centerPx, ...
        'Area', pi * effectiveRadiusPx^2, ...
        'Radius_pixels', effectiveRadiusPx, ...
        'Radius_nm', effectiveRadiusPx * avgPixelSize_nm, ...
        'CoreRadius_pixels', radiusPx, ...
        'CoreRadius_nm', radiusPx * avgPixelSize_nm, ...
        'OuterRadius_pixels', radiusInfo.outerRadiusPx, ...
        'OuterRadius_nm', radiusInfo.outerRadiusPx * avgPixelSize_nm, ...
        'EdgeRadius_pixels', radiusInfo.edgeRadiusPx, ...
        'ContrastRadius_pixels', radiusInfo.contrastRadiusPx, ...
        'CenterEnhanced', centerEnhanced, ...
        'RingEnhanced', ringEnhanced, ...
        'OuterEnhanced', outerEnhanced, ...
        'CenterFraction', centerFraction, ...
        'RingFraction', ringFraction, ...
        'OuterFraction', outerFraction, ...
        'CenterRatio', centerRatio, ...
        'RingRatio', ringRatio, ...
        'OuterRatio', outerRatio, ...
        'FractionalContrast', fractionalContrast, ...
        'RatioContrast', ratioContrast, ...
        'EdgeClearancePx', edgeClearancePx, ...
        'OrientationDeg', orientationDeg, ...
        'MeanIntensity', localTotal, ...
        'IntensityStd', std(featureMaps.totalMap(distMap <= radiusPx), 0, 'omitnan'), ...
        'Circularity', 1, ...
        'Solidity', 1, ...
        'AspectRatio', 1, ...
        'BoundingBox', [centerPx(1) - effectiveRadiusPx, centerPx(2) - effectiveRadiusPx, ...
        2 * effectiveRadiusPx, 2 * effectiveRadiusPx], ...
        'DetectionScore', qualityScore);
end

function [effectiveRadiusPx, radiusInfo] = estimateVpitOuterRadiusPx(centerPx, coreRadiusPx, featureMaps, distMap)
    rows = size(featureMaps.enhancedMap, 1);
    cols = size(featureMaps.enhancedMap, 2);
    edgeClearancePx = min([centerPx(1) - 1, centerPx(2) - 1, cols - centerPx(1), rows - centerPx(2)]);
    searchMinPx = max(1.0, 0.55 * coreRadiusPx);
    searchMaxPx = min([max(searchMinPx + 0.5, 2.10 * coreRadiusPx), edgeClearancePx + 0.75, 0.48 * min(rows, cols)]);
    if ~isfinite(searchMaxPx) || searchMaxPx <= searchMinPx
        effectiveRadiusPx = coreRadiusPx;
        radiusInfo = buildRadiusInfo(coreRadiusPx, coreRadiusPx, coreRadiusPx, coreRadiusPx);
        return;
    end

    radialSamplesPx = linspace(searchMinPx, searchMaxPx, 36);
    fractionProfile = sampleAnnularMean(featureMaps.fractionMap, distMap, radialSamplesPx);
    edgeProfile = sampleAnnularMean(featureMaps.edgeMap, distMap, radialSamplesPx);

    outerMask = distMap > 1.55 * coreRadiusPx & distMap <= min(searchMaxPx, 2.25 * coreRadiusPx);
    if ~any(outerMask(:))
        outerMask = distMap > 1.25 * coreRadiusPx & distMap <= searchMaxPx;
    end
    outerFraction = mean(featureMaps.fractionMap(outerMask), 'omitnan');
    innerMask = distMap <= max(1.1, 0.80 * coreRadiusPx);
    innerFraction = mean(featureMaps.fractionMap(innerMask), 'omitnan');
    fractionContrast = max(innerFraction - outerFraction, 0);

    if ~isfinite(fractionContrast) || fractionContrast <= 1e-8
        contrastRadiusPx = coreRadiusPx;
    else
        contrastThreshold = outerFraction + 0.22 * fractionContrast;
        validContrastMask = fractionProfile >= contrastThreshold;
        if any(validContrastMask)
            contrastRadiusPx = max(radialSamplesPx(validContrastMask));
        else
            contrastRadiusPx = coreRadiusPx;
        end
    end

    edgeSearchMask = radialSamplesPx >= 0.78 * coreRadiusPx & radialSamplesPx <= searchMaxPx;
    if any(edgeSearchMask) && any(isfinite(edgeProfile(edgeSearchMask)))
        candidateRadii = radialSamplesPx(edgeSearchMask);
        candidateEdges = edgeProfile(edgeSearchMask);
        [~, edgeIdx] = max(candidateEdges);
        edgeRadiusPx = candidateRadii(edgeIdx);
    else
        edgeRadiusPx = coreRadiusPx;
    end

    if ~isfinite(edgeRadiusPx)
        edgeRadiusPx = coreRadiusPx;
    end
    if ~isfinite(contrastRadiusPx)
        contrastRadiusPx = coreRadiusPx;
    end

    outerRadiusPx = max([coreRadiusPx, 0.65 * edgeRadiusPx + 0.35 * contrastRadiusPx, contrastRadiusPx]);
    upperBoundPx = min(searchMaxPx, max(coreRadiusPx + 0.75, 1.45 * coreRadiusPx));
    effectiveRadiusPx = min(max(coreRadiusPx, outerRadiusPx), upperBoundPx);
    radiusInfo = buildRadiusInfo(coreRadiusPx, effectiveRadiusPx, edgeRadiusPx, contrastRadiusPx);
end

function radiusInfo = buildRadiusInfo(coreRadiusPx, outerRadiusPx, edgeRadiusPx, contrastRadiusPx)
    radiusInfo = struct( ...
        'coreRadiusPx', coreRadiusPx, ...
        'outerRadiusPx', outerRadiusPx, ...
        'edgeRadiusPx', edgeRadiusPx, ...
        'contrastRadiusPx', contrastRadiusPx);
end

function profile = sampleAnnularMean(mapData, distMap, radialSamplesPx)
    profile = nan(size(radialSamplesPx));
    if isempty(mapData)
        return;
    end
    for sampleIdx = 1:numel(radialSamplesPx)
        radiusPx = radialSamplesPx(sampleIdx);
        halfWidthPx = max(0.65, 0.055 * max(radiusPx, 1));
        annulusMask = abs(distMap - radiusPx) <= halfWidthPx;
        if any(annulusMask(:))
            profile(sampleIdx) = mean(mapData(annulusMask), 'omitnan');
        end
    end
end

function [centroids, radii, properties] = postProcessDetections(centroids, radii, properties, enhancedMap)
    if isempty(centroids)
        return;
    end
    
    i = 1;
    while i <= size(centroids, 1)
        j = i + 1;
        while j <= size(centroids, 1)
            distance = norm(centroids(i, :) - centroids(j, :));
            mergeDistance = 0.45 * max(properties(i).Radius_pixels, properties(j).Radius_pixels);
            if distance < mergeDistance
                score_i = getFieldOrFallback(properties(i), 'DetectionScore', 0) + ...
                    0.15 * getFieldOrFallback(properties(i), 'Circularity', 0) * getFieldOrFallback(properties(i), 'Solidity', 0);
                score_j = getFieldOrFallback(properties(j), 'DetectionScore', 0) + ...
                    0.15 * getFieldOrFallback(properties(j), 'Circularity', 0) * getFieldOrFallback(properties(j), 'Solidity', 0);
                
                if score_i >= score_j
                    centroids(j, :) = [];
                    radii(j) = [];
                    properties(j) = [];
                else
                    centroids(i, :) = centroids(j, :);
                    radii(i) = radii(j);
                    properties(i) = properties(j);
                    centroids(j, :) = [];
                    radii(j) = [];
                    properties(j) = [];
                end
            else
                j = j + 1;
            end
        end
        i = i + 1;
    end
    validIndices = true(size(radii));
    for i = 1:length(radii)
        centroid = round(centroids(i, :));
        if centroid(1) >= 1 && centroid(1) <= size(enhancedMap, 2) && ...
           centroid(2) >= 1 && centroid(2) <= size(enhancedMap, 1)
            localIntensity = enhancedMap(centroid(2), centroid(1));
            if localIntensity < 0.05
                validIndices(i) = false;
            end
        end
    end

    centroids = centroids(validIndices, :);
    radii = radii(validIndices);
    properties = properties(validIndices);
    if isempty(centroids)
        return;
    end

    [centroids, radii, properties] = suppressOverlappingVpits(centroids, radii, properties);
    if numel(radii) <= 1
        return;
    end

    detectionScores = extractDetectionValues(properties, 'DetectionScore', 0);
    fractionContrasts = extractDetectionValues(properties, 'FractionalContrast', 0);
    edgeClearance = extractDetectionValues(properties, 'EdgeClearancePx', inf);
    radiusPixels = extractDetectionValues(properties, 'Radius_pixels', 1);
    maxScore = max(detectionScores);
    maxContrast = max(fractionContrasts);

    strongScore = detectionScores >= 0.50 * maxScore;
    strongContrast = fractionContrasts >= 0.28 * max(maxContrast, eps);
    usableEdge = edgeClearance >= 0.25 * radiusPixels | detectionScores >= 0.78 * maxScore;
    keepMask = (strongScore & strongContrast & usableEdge) | detectionScores >= 0.82 * maxScore;
    if ~any(keepMask)
        [~, bestIdx] = max(detectionScores);
        keepMask(bestIdx) = true;
    end

    centroids = centroids(keepMask, :);
    radii = radii(keepMask);
    properties = properties(keepMask);
end

function [centroids, radii, properties] = suppressOverlappingVpits(centroids, radii, properties)
    if numel(radii) <= 1
        return;
    end
    detectionScores = extractDetectionValues(properties, 'DetectionScore', 0);
    radiusPixels = extractDetectionValues(properties, 'Radius_pixels', 1);
    [~, order] = sort(detectionScores, 'descend');
    keepMask = false(size(radii));

    for orderIdx = 1:numel(order)
        candidateIdx = order(orderIdx);
        keptIdx = find(keepMask);
        if isempty(keptIdx)
            keepMask(candidateIdx) = true;
            continue;
        end
        distances = sqrt(sum((centroids(keptIdx, :) - centroids(candidateIdx, :)).^2, 2));
        minSeparation = 0.72 * (radiusPixels(keptIdx) + radiusPixels(candidateIdx));
        if all(distances > minSeparation)
            keepMask(candidateIdx) = true;
        end
    end

    centroids = centroids(keepMask, :);
    radii = radii(keepMask);
    properties = properties(keepMask);
end

function values = extractDetectionValues(properties, fieldName, defaultValue)
    values = defaultValue * ones(numel(properties), 1);
    for idx = 1:numel(properties)
        if isfield(properties(idx), fieldName) && ~isempty(properties(idx).(fieldName))
            values(idx) = properties(idx).(fieldName);
        end
    end
end

function value = getFieldOrFallback(structValue, fieldName, fallbackValue)
    if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
        value = structValue.(fieldName);
    else
        value = fallbackValue;
    end
end

function orientationDeg = estimateHexagonOrientationDeg(featureMaps, centerPx, radiusPx)
    if ~isfield(featureMaps, 'edgeMap') || isempty(featureMaps.edgeMap)
        orientationDeg = 0;
        return;
    end

    edgeMap = featureMaps.edgeMap;
    rows = size(edgeMap, 1);
    cols = size(edgeMap, 2);
    apothemPx = max(1.5, 0.92 * radiusPx);
    samplesPerEdge = 7;
    phaseCandidates = 0:2:58;
    scores = -inf(size(phaseCandidates));

    for phaseIdx = 1:numel(phaseCandidates)
        phaseDeg = phaseCandidates(phaseIdx);
        edgeScore = 0;
        vertexPenalty = 0;
        for facetIdx = 0:5
            normalDeg = phaseDeg + 60 * facetIdx;
            tangentDeg = normalDeg + 90;
            sideCenter = centerPx + apothemPx * [cosd(normalDeg), sind(normalDeg)];
            halfSidePx = apothemPx * tand(30);
            edgeValues = zeros(1, samplesPerEdge);
            for sampleIdx = 1:samplesPerEdge
                sampleOffset = -halfSidePx + 2 * halfSidePx * (sampleIdx - 1) / max(samplesPerEdge - 1, 1);
                samplePoint = sideCenter + sampleOffset * [cosd(tangentDeg), sind(tangentDeg)];
                edgeValues(sampleIdx) = bilinearSample(edgeMap, samplePoint(1), samplePoint(2), cols, rows);
            end
            edgeScore = edgeScore + mean(edgeValues, 'omitnan');

            vertexDeg = normalDeg + 30;
            vertexRadiusPx = apothemPx / cosd(30);
            vertexPoint = centerPx + vertexRadiusPx * [cosd(vertexDeg), sind(vertexDeg)];
            vertexPenalty = vertexPenalty + bilinearSample(edgeMap, vertexPoint(1), vertexPoint(2), cols, rows);
        end
        scores(phaseIdx) = edgeScore - 0.20 * vertexPenalty;
    end

    [~, bestIdx] = max(scores);
    orientationDeg = phaseCandidates(bestIdx);
end

function sampleValue = bilinearSample(imageData, xCoord, yCoord, maxX, maxY)
    if nargin < 4
        maxY = size(imageData, 1);
        maxX = size(imageData, 2);
    end
    if xCoord < 1 || yCoord < 1 || xCoord > maxX || yCoord > maxY
        sampleValue = 0;
        return;
    end

    x1 = floor(xCoord);
    x2 = min(maxX, ceil(xCoord));
    y1 = floor(yCoord);
    y2 = min(maxY, ceil(yCoord));
    dx = xCoord - x1;
    dy = yCoord - y1;

    if x1 < 1 || y1 < 1
        sampleValue = 0;
        return;
    end

    q11 = imageData(y1, x1);
    q21 = imageData(y1, x2);
    q12 = imageData(y2, x1);
    q22 = imageData(y2, x2);

    sampleValue = (1 - dx) * (1 - dy) * q11 + ...
                  dx * (1 - dy) * q21 + ...
                  (1 - dx) * dy * q12 + ...
                  dx * dy * q22;
end

function values = extractPropertyVector(properties, fieldName, expectedCount, defaultValue)
    if nargin < 4
        defaultValue = 0;
    end
    values = defaultValue * ones(expectedCount, 1);
    if isempty(properties)
        return;
    end
    count = min(numel(properties), expectedCount);
    for idx = 1:count
        if isfield(properties(idx), fieldName) && ~isempty(properties(idx).(fieldName))
            values(idx) = properties(idx).(fieldName);
        end
    end
end

