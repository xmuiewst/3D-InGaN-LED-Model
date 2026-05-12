function cfg = config(entryMode, params)
    if nargin < 1 || isempty(entryMode)
        entryMode = 'wavelength';
    end
    if nargin < 2
        params = struct();
    end

    entryMode = lower(string(entryMode));
    if entryMode == "wav"
        entryMode = "wavelength";
    elseif entryMode ~= "log"
        entryMode = "wavelength";
    end

    cfg = struct();
    cfg.entryMode = char(entryMode);
    cfg.default3DMode = default3DMode(entryMode);
    cfg.available3DModes = {'intensity', 'log_intensity', 'wavelength'};
    cfg.layerParams = buildLayerParams(entryMode, params);
end

function modeName = default3DMode(entryMode)
    switch entryMode
        case "log"
            modeName = 'log_intensity';
        otherwise
            modeName = 'wavelength';
    end
end

function layerParams = buildLayerParams(~, params)
    calibrationProfile = getFieldOrDefault(params, 'calibrationProfile', struct());

    barrierThick = getProfileOrDefault(calibrationProfile, 'barrierThick', ...
        getFieldOrDefault(params, 'barrierThick', 8));
    wellThick = getProfileOrDefault(calibrationProfile, 'wellThick', ...
        getFieldOrDefault(params, 'wellThick', 5));

    nGaNSpectrum = [684, 390.90e-9, 10.00e-9, 15, 376.71e-9, 3.74e-9];
    prestrainedSpectrum = [12, 414.0e-9, 34.0e-9, 4, 430.0e-9, 45.0e-9];
    barrierSpectrum = [42, 478.0e-9, 28.0e-9, 24, 503.0e-9, 30.0e-9];
    wellSpectrum = [76, 514.0e-9, 28.0e-9, 36, 540.0e-9, 24.0e-9];
    defaultPTypeSpectrum = [72, 388.5e-9, 11.0e-9, 22, 398.0e-9, 16.0e-9];
    pTypeSpectrum = defaultPTypeSpectrum;

    nGaNSpectrum = getProfileOrDefault(calibrationProfile, 'nGaNSpectrum', ...
        getProfileOrDefault(calibrationProfile, 'gaNSpectrum', nGaNSpectrum));
    prestrainedSpectrum = getProfileOrDefault(calibrationProfile, 'prestrainedSpectrum', prestrainedSpectrum);
    barrierSpectrum = getProfileOrDefault(calibrationProfile, 'barrierSpectrum', barrierSpectrum);
    wellSpectrum = getProfileOrDefault(calibrationProfile, 'wellSpectrum', wellSpectrum);
    pTypeSpectrum = getProfileOrDefault(calibrationProfile, 'pTypeSpectrum', pTypeSpectrum);
    pTypeSpectrum = enforcePTypeNearBandEdgeSpectrum(pTypeSpectrum, defaultPTypeSpectrum);

    prestrainedIn = getProfileOrDefault(calibrationProfile, 'prestrainedInComposition', 0.0530);
    barrierIn = getProfileOrDefault(calibrationProfile, 'barrierInComposition', 0.1458);
    wellIn = getProfileOrDefault(calibrationProfile, 'wellInComposition', 0.6647);
    eblIn = getProfileOrDefault(calibrationProfile, 'eblInComposition', 0.2130);

    layerParams = {
        'Substrate',   0.0000, 10,           'none',   [], [], [], [0, 0, 0, 0, 0, 0], 0.4083
        'n-GaN',       0.0000, 2000,         'n-type', [], [], [], nGaNSpectrum,         3
        'Prestrained', prestrainedIn, 30,           'none',   [], [], [], prestrainedSpectrum, 0.08
        'MQW-Barrier', barrierIn, barrierThick, 'none',   [], [], [], barrierSpectrum,      0.18
        'MQW-Well',    wellIn, wellThick,    'none',   [], [], [], wellSpectrum,         0.72
        'p-EBL',       eblIn, 0,            'p-type', [], [], [], pTypeSpectrum,        0.16
        'p-GaN',       0.0000, 10,           'p-type', [], [], [], pTypeSpectrum,        0.30
    };
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
    if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end

function value = getProfileOrDefault(profile, fieldName, defaultValue)
    if isstruct(profile) && isfield(profile, fieldName) && ~isempty(profile.(fieldName))
        value = profile.(fieldName);
    else
        value = defaultValue;
    end
end

function spectrumOut = enforcePTypeNearBandEdgeSpectrum(spectrumIn, fallbackSpectrum)
    spectrumOut = spectrumIn;
    components = [];
    if isstruct(spectrumIn) && isfield(spectrumIn, 'components') && ~isempty(spectrumIn.components)
        components = double(spectrumIn.components);
    elseif isnumeric(spectrumIn) && numel(spectrumIn) >= 6
        rawValues = double(spectrumIn(:)');
        componentCount = floor(numel(rawValues) / 3);
        components = reshape(rawValues(1:3 * componentCount), 3, [])';
    end

    if isempty(components) || size(components, 2) < 3
        spectrumOut = fallbackSpectrum;
        return;
    end

    [~, dominantIdx] = max(abs(components(:, 1)));
    dominantCenterNm = components(dominantIdx, 2) * 1e9;
    if isfinite(dominantCenterNm) && dominantCenterNm >= 365 && dominantCenterNm <= 402
        return;
    end

    spectrumOut = fallbackSpectrum;
end
