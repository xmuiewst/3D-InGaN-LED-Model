function varargout = run_3d_modeling(options)
if nargin < 1 || isempty(options)
    options = struct();
end

rootDir = fileparts(mfilename('fullpath'));
opts = defaultOptions(rootDir);
opts = mergeOptions(opts, options);

if ~isfile(opts.matFile)
    error('PointCloudVPitMQW:MissingMatFile', 'MAT file does not exist: %s', opts.matFile);
end
if ~isfolder(opts.outputDir)
    mkdir(opts.outputDir);
end

activityPackage = loadActivityPackageLocal(opts.matFile);
surfaceZNm = resolveSurfaceZ(activityPackage);
vPits = extractVpitGeometry(activityPackage, surfaceZNm);
[xGridNm, yGridNm, depthGridNm, scanRangeNm] = buildVoxelAxes( ...
    activityPackage, surfaceZNm, opts);
mqwLayers = extractMqwLayers(activityPackage, surfaceZNm, max(depthGridNm));

[voxelVolume, voxelCounts, pointCount] = voxelizePointCloud( ...
    activityPackage.points, xGridNm, yGridNm, depthGridNm, surfaceZNm, opts);
[voxelVolume, cavityMask, shellMask] = applyVpitGeometry( ...
    voxelVolume, xGridNm, yGridNm, depthGridNm, vPits);

displayVolume = voxelVolume .^ 0.72;
result = struct( ...
    'method', 'Point-cloud voxel block with V-pit and MQW geometry', ...
    'sourceMatFile', opts.matFile, ...
    'surfaceZNm', surfaceZNm, ...
    'scanRangeNm', scanRangeNm, ...
    'xGridNm', xGridNm, ...
    'yGridNm', yGridNm, ...
    'depthGridNm', depthGridNm, ...
    'voxelVolume', voxelVolume, ...
    'displayVolume', displayVolume, ...
    'voxelCounts', voxelCounts, ...
    'pointCountUsed', pointCount, ...
    'cavityMask', cavityMask, ...
    'shellMask', shellMask, ...
    'vPits', vPits, ...
    'mqwLayers', mqwLayers, ...
    'options', opts);

result.outputMat = fullfile(opts.outputDir, '3d_modeling_volume.mat');
save(result.outputMat, 'result', '-v7.3');

visibility = figureVisibility(opts.showFigures);
figBlock = renderBlockFigure(result, visibility);
result.blockFigure = fullfile(opts.outputDir, '3d_modeling_block.png');
exportgraphics(figBlock, result.blockFigure, 'Resolution', 320);

figSection = renderCrossSectionFigure(result, visibility);
result.crossSectionFigure = fullfile(opts.outputDir, '3d_modeling_cross_section.png');
exportgraphics(figSection, result.crossSectionFigure, 'Resolution', 320);

if ~opts.showFigures
    close([figBlock, figSection]);
end

fprintf('Point-cloud V-pit/MQW voxel model saved to %s\n', result.outputMat);
fprintf('Block rendering: %s\n', result.blockFigure);
fprintf('Cross section: %s\n', result.crossSectionFigure);
if nargout > 0
    varargout{1} = result;
end
end

function opts = defaultOptions(rootDir)
opts = struct( ...
    'matFile', fullfile(rootDir, 'results', 'three_dim_model_package.mat'), ...
    'outputDir', fullfile(rootDir, 'output', '3d_modeling'), ...
    'gridSize', [110, 92, 110], ...
    'depthMaxNm', [], ...
    'maxPoints', 250000, ...
    'smoothVolume', true, ...
    'externalAlpha', 0.62, ...
    'vPitAlpha', 0.88, ...
    'showFigures', true);
end

function opts = mergeOptions(opts, options)
if ~isstruct(options)
    error('PointCloudVPitMQW:InvalidOptions', 'Options must be a struct.');
end
names = fieldnames(options);
for idx = 1:numel(names)
    opts.(names{idx}) = options.(names{idx});
end
end

function activityPackage = loadActivityPackageLocal(matFile)
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
error('PointCloudVPitMQW:InvalidPackage', ...
    'The MAT file does not contain an Activity-Structure package.');
end

function surfaceZNm = resolveSurfaceZ(activityPackage)
surfaceZNm = NaN;
if isfield(activityPackage, 'vPits') && istable(activityPackage.vPits) && ...
        ismember('SurfaceZ_nm', activityPackage.vPits.Properties.VariableNames)
    values = double(activityPackage.vPits.SurfaceZ_nm);
    values = values(isfinite(values));
    if ~isempty(values)
        surfaceZNm = max(values);
    end
end
if ~isfinite(surfaceZNm) && isfield(activityPackage, 'layers') && istable(activityPackage.layers)
    values = double(activityPackage.layers.ZTop_nm);
    values = values(isfinite(values));
    if ~isempty(values)
        surfaceZNm = max(values);
    end
end
if ~isfinite(surfaceZNm)
    surfaceZNm = max(double(activityPackage.points.Z_nm), [], 'omitnan');
end
end

function [xGridNm, yGridNm, depthGridNm, scanRangeNm] = buildVoxelAxes( ...
    activityPackage, surfaceZNm, opts)
gridSize = round(double(opts.gridSize(:)'));
if numel(gridSize) ~= 3 || any(gridSize < 16)
    error('PointCloudVPitMQW:InvalidGridSize', 'gridSize must be [Nx Ny Nz], with each value at least 16.');
end

if isfield(activityPackage, 'xGrid_nm') && ~isempty(activityPackage.xGrid_nm)
    xLimits = [min(activityPackage.xGrid_nm), max(activityPackage.xGrid_nm)];
else
    xLimits = prctile(double(activityPackage.points.X_nm), [2, 98]);
end
if isfield(activityPackage, 'yGrid_nm') && ~isempty(activityPackage.yGrid_nm)
    yLimits = [min(activityPackage.yGrid_nm), max(activityPackage.yGrid_nm)];
else
    yLimits = prctile(double(activityPackage.points.Y_nm), [2, 98]);
end

xGridNm = linspace(xLimits(1), xLimits(2), gridSize(1));
yGridNm = linspace(yLimits(1), yLimits(2), gridSize(2));

if isempty(opts.depthMaxNm)
    points = activityPackage.points;
    insideScan = double(points.X_nm) >= xLimits(1) & double(points.X_nm) <= xLimits(2) & ...
        double(points.Y_nm) >= yLimits(1) & double(points.Y_nm) <= yLimits(2) & ...
        isfinite(double(points.Z_nm));
    depths = surfaceZNm - double(points.Z_nm(insideScan));
    depths = depths(isfinite(depths) & depths >= 0);
    if isempty(depths)
        depthMaxNm = surfaceZNm - min(double(points.Z_nm), [], 'omitnan');
    else
        depthMaxNm = max(depths);
    end
else
    depthMaxNm = double(opts.depthMaxNm);
end
depthGridNm = linspace(0, depthMaxNm, gridSize(3));
scanRangeNm = struct( ...
    'xMin', xLimits(1), ...
    'xMax', xLimits(2), ...
    'yMin', yLimits(1), ...
    'yMax', yLimits(2), ...
    'depthMin', 0, ...
    'depthMax', depthMaxNm);
end

function [volume, counts, pointCount] = voxelizePointCloud( ...
    pointTable, xGridNm, yGridNm, depthGridNm, surfaceZNm, opts)
required = {'X_nm', 'Y_nm', 'Z_nm', 'Intensity'};
if ~all(ismember(required, pointTable.Properties.VariableNames))
    error('PointCloudVPitMQW:InvalidPointTable', ...
        'Point table must contain X_nm, Y_nm, Z_nm, and Intensity.');
end

x = double(pointTable.X_nm);
y = double(pointTable.Y_nm);
depth = surfaceZNm - double(pointTable.Z_nm);
intensity = double(pointTable.Intensity);
valid = isfinite(x) & isfinite(y) & isfinite(depth) & isfinite(intensity) & intensity > 0 & ...
    x >= min(xGridNm) & x <= max(xGridNm) & ...
    y >= min(yGridNm) & y <= max(yGridNm) & ...
    depth >= min(depthGridNm) & depth <= max(depthGridNm);
x = x(valid);
y = y(valid);
depth = depth(valid);
intensity = intensity(valid);

maxPoints = max(1000, round(double(opts.maxPoints)));
if numel(x) > maxPoints
    indices = unique(round(linspace(1, numel(x), maxPoints)));
    x = x(indices);
    y = y(indices);
    depth = depth(indices);
    intensity = intensity(indices);
end
pointCount = numel(x);
if pointCount < 10
    error('PointCloudVPitMQW:InsufficientPoints', ...
        'Too few point-cloud samples remain inside the voxel domain.');
end

xEdges = gridEdges(xGridNm);
yEdges = gridEdges(yGridNm);
depthEdges = gridEdges(depthGridNm);
ix = discretize(x, xEdges);
iy = discretize(y, yEdges);
iz = discretize(depth, depthEdges);
validBins = isfinite(ix) & isfinite(iy) & isfinite(iz);

positiveIntensity = intensity(validBins);
referenceIntensity = prctile(positiveIntensity, 35);
if ~isfinite(referenceIntensity) || referenceIntensity <= eps
    referenceIntensity = median(positiveIntensity, 'omitnan');
end
logIntensity = log1p(positiveIntensity ./ max(referenceIntensity, eps));
dims = [numel(yGridNm), numel(xGridNm), numel(depthGridNm)];
subs = [iy(validBins), ix(validBins), iz(validBins)];
volume = accumarray(subs, logIntensity, dims, @mean, 0);
counts = accumarray(subs, 1, dims, @sum, 0);

if opts.smoothVolume
    volume = smooth3(volume, 'gaussian', [5, 5, 5], 1.0);
    support = smooth3(double(counts > 0), 'box', [3, 3, 3]);
    volume(support <= 0.01) = 0;
end
volume = normalizeVolume(volume);
end

function edges = gridEdges(gridValues)
step = median(diff(gridValues));
edges = [gridValues(1) - 0.5 * step, ...
    0.5 * (gridValues(1:end-1) + gridValues(2:end)), ...
    gridValues(end) + 0.5 * step];
end

function volume = normalizeVolume(volume)
volume(~isfinite(volume) | volume < 0) = 0;
values = volume(volume > 0);
if isempty(values)
    return;
end
lowValue = prctile(values, 2);
highValue = prctile(values, 99);
if ~isfinite(highValue) || highValue <= lowValue
    highValue = max(values);
end
volume = (volume - lowValue) ./ max(highValue - lowValue, eps);
volume = min(max(volume, 0), 1);
end

function vPits = extractVpitGeometry(activityPackage, surfaceZNm)
vPits = repmat(emptyVpit(), 0, 1);
if ~isfield(activityPackage, 'vPits') || ~istable(activityPackage.vPits)
    return;
end
tbl = activityPackage.vPits;
required = {'CenterX_nm', 'CenterY_nm', 'TopRadius_nm', 'Depth_nm'};
if ~all(ismember(required, tbl.Properties.VariableNames))
    return;
end
for idx = 1:height(tbl)
    pit = emptyVpit();
    pit.index = idx;
    pit.centerNm = [double(tbl.CenterX_nm(idx)), double(tbl.CenterY_nm(idx))];
    pit.topApothemNm = double(tbl.TopRadius_nm(idx));
    pit.depthNm = double(tbl.Depth_nm(idx));
    pit.surfaceZNm = surfaceZNm;
    if ismember('Orientation_deg', tbl.Properties.VariableNames)
        pit.orientationDeg = double(tbl.Orientation_deg(idx));
    end
    if ~isfinite(pit.orientationDeg)
        pit.orientationDeg = 0;
    end
    pit.topVerticesNm = hexagonFromApothem(pit.centerNm, pit.topApothemNm, pit.orientationDeg);
    vPits(end + 1, 1) = pit;
end
end

function pit = emptyVpit()
pit = struct( ...
    'index', 0, ...
    'centerNm', [NaN, NaN], ...
    'topApothemNm', NaN, ...
    'depthNm', NaN, ...
    'orientationDeg', 0, ...
    'surfaceZNm', NaN, ...
    'topVerticesNm', zeros(0, 2));
end

function layers = extractMqwLayers(activityPackage, surfaceZNm, depthMaxNm)
layers = repmat(struct( ...
    'name', '', ...
    'zBottomNm', NaN, ...
    'zTopNm', NaN, ...
    'depthTopNm', NaN, ...
    'depthBottomNm', NaN, ...
    'depthCenterNm', NaN), 0, 1);
if ~isfield(activityPackage, 'layers') || ~istable(activityPackage.layers)
    return;
end
tbl = activityPackage.layers;
required = {'LayerName', 'ZBottom_nm', 'ZTop_nm'};
if ~all(ismember(required, tbl.Properties.VariableNames))
    return;
end
mask = contains(string(tbl.LayerName), 'MQW-Well', 'IgnoreCase', true);
rows = find(mask);
for idx = 1:numel(rows)
    row = rows(idx);
    depthTop = surfaceZNm - double(tbl.ZTop_nm(row));
    depthBottom = surfaceZNm - double(tbl.ZBottom_nm(row));
    if depthBottom < 0 || depthTop > depthMaxNm
        continue;
    end
    layers(end + 1, 1) = struct( ...
        'name', char(string(tbl.LayerName(row))), ...
        'zBottomNm', double(tbl.ZBottom_nm(row)), ...
        'zTopNm', double(tbl.ZTop_nm(row)), ...
        'depthTopNm', depthTop, ...
        'depthBottomNm', depthBottom, ...
        'depthCenterNm', 0.5 * (depthTop + depthBottom));
end
end

function [volume, cavityMask, shellMask] = applyVpitGeometry( ...
    volume, xGridNm, yGridNm, depthGridNm, vPits)
cavityMask = false(size(volume));
shellMask = false(size(volume));
if isempty(vPits)
    return;
end

[X, Y] = meshgrid(xGridNm, yGridNm);
for pitIdx = 1:numel(vPits)
    pit = vPits(pitIdx);
    rho = hexagonalRadius(X, Y, pit);
    for depthIdx = 1:numel(depthGridNm)
        depthNm = depthGridNm(depthIdx);
        if depthNm > pit.depthNm
            continue;
        end
        cavityRadius = max(0, 1 - depthNm / max(pit.depthNm, eps));
        innerMask = rho <= 0.82 * cavityRadius;
        outerMask = rho <= min(1, 1.10 * cavityRadius);
        localShell = outerMask & ~innerMask;
        cavityMask(:, :, depthIdx) = cavityMask(:, :, depthIdx) | innerMask;
        shellMask(:, :, depthIdx) = shellMask(:, :, depthIdx) | localShell;

        slice = volume(:, :, depthIdx);
        shellValues = slice(localShell & slice > 0);
        if isempty(shellValues)
            shellLevel = 0.78;
        else
            shellLevel = max(0.70, prctile(shellValues, 75));
        end
        slice(localShell) = max(slice(localShell), shellLevel);
        slice(innerMask) = 0;
        volume(:, :, depthIdx) = slice;
    end
end
end

function rho = hexagonalRadius(X, Y, pit)
dx = X - pit.centerNm(1);
dy = Y - pit.centerNm(2);
facetAngles = pit.orientationDeg + (0:60:300);
rho = -inf(size(X));
for idx = 1:numel(facetAngles)
    projected = dx * cosd(facetAngles(idx)) + dy * sind(facetAngles(idx));
    rho = max(rho, projected ./ max(pit.topApothemNm, eps));
end
end

function vertices = hexagonFromApothem(centerNm, apothemNm, orientationDeg)
vertexRadiusNm = apothemNm / cosd(30);
angles = orientationDeg + (30:60:330);
vertices = [ ...
    centerNm(1) + vertexRadiusNm * cosd(angles(:)), ...
    centerNm(2) + vertexRadiusNm * sind(angles(:))];
end

function visibility = figureVisibility(showFigures)
if showFigures
    visibility = 'on';
else
    visibility = 'off';
end
end

function fig = renderBlockFigure(result, visibility)
fig = figure('Color', 'w', 'Visible', visibility, ...
    'Position', [80, 60, 1280, 820], 'Name', 'Point-cloud V-pit and MQW voxel block');
plotPosition = [0.06, 0.12, 0.68, 0.74];
ax = axes(fig, 'Position', plotPosition);
hold(ax, 'on');

x = result.xGridNm;
y = result.yGridNm;
depth = result.depthGridNm;
displayVolume = 10 * result.displayVolume;
externalAlpha = clampAlpha(result.options.externalAlpha);
vPitAlpha = clampAlpha(result.options.vPitAlpha);

[Xtop, Ytop] = meshgrid(x, y);
surfaceDepth = buildPitSurfaceDepth(Xtop, Ytop, result.vPits);
surfaceSample = interp3(x, y, depth, displayVolume, Xtop, Ytop, surfaceDepth, 'linear', 0);
shallowMask = depth <= min(105, max(depth));
topProjection = max(displayVolume(:, :, shallowMask), [], 3);
topColor = min(10, 0.58 * topProjection + 0.42 * surfaceSample);

topAlpha = externalAlpha * ones(size(surfaceDepth));
topAlpha(surfaceDepth > 0.5) = vPitAlpha;
topSurface = surf(ax, Xtop, Ytop, surfaceDepth, topColor, ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 'interp');
topSurface.AlphaData = topAlpha;
topSurface.AlphaDataMapping = 'none';

[Xfront, Dfront] = meshgrid(x, depth);
Yfront = min(y) * ones(size(Xfront));
frontColor = squeeze(displayVolume(1, :, :))';
surf(ax, Xfront, Yfront, Dfront, frontColor, ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', externalAlpha);

[Xback, Dback] = meshgrid(x, depth);
Yback = max(y) * ones(size(Xback));
backColor = squeeze(displayVolume(end, :, :))';
surf(ax, Xback, Yback, Dback, backColor, ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', externalAlpha);

[Yright, Dright] = meshgrid(y, depth);
Xright = max(x) * ones(size(Yright));
rightColor = squeeze(displayVolume(:, end, :))';
surf(ax, Xright, Yright, Dright, rightColor, ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', externalAlpha);

[Yleft, Dleft] = meshgrid(y, depth);
Xleft = min(x) * ones(size(Yleft));
leftColor = squeeze(displayVolume(:, 1, :))';
surf(ax, Xleft, Yleft, Dleft, leftColor, ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', externalAlpha);

[Xbottom, Ybottom] = meshgrid(x, y);
Dbottom = max(depth) * ones(size(Xbottom));
bottomColor = displayVolume(:, :, end);
surf(ax, Xbottom, Ybottom, Dbottom, bottomColor, ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', externalAlpha);

colormap(ax, turbo(256));
clim(ax, [0, 10]);
cb = colorbar(ax);
cb.Position = [0.72, 0.21, 0.022, 0.58];
cb.Label.String = 'Normalized point-cloud intensity';
ax.Position = plotPosition;

xlabel(ax, 'X (nm)');
ylabel(ax, 'Y (nm)');
zlabel(ax, 'Depth below surface (nm)');
set(ax, 'ZDir', 'reverse', 'Projection', 'orthographic');
axis(ax, 'vis3d');
ax.SortMethod = 'depth';
daspect(ax, [1, 1, 0.72]);
xlim(ax, [min(x), max(x)]);
ylim(ax, [min(y), max(y)]);
zlim(ax, [min(depth), max(depth)]);
view(ax, 38, 27);
fitOrthographicCamera(ax, 1.22);
grid(ax, 'off');
box(ax, 'off');
ax.Clipping = 'on';
camlight(ax, 'headlight');
camlight(ax, 'right');
lighting(ax, 'gouraud');
rotation = rotate3d(fig);
rotation.ActionPostCallback = @(~, event) fitOrthographicCamera(event.Axes, 1.22);
rotation.Enable = 'on';
fig.SizeChangedFcn = @(~, ~) fitOrthographicCamera(ax, 1.22);
end

function fitOrthographicCamera(ax, margin)
ax.CameraViewAngleMode = 'auto';
drawnow;
viewAngle = ax.CameraViewAngle;
ax.CameraViewAngleMode = 'manual';
ax.CameraViewAngle = min(179, viewAngle * margin);
end

function value = clampAlpha(value)
value = double(value);
if isempty(value) || ~isfinite(value(1))
    value = 0.2;
else
    value = min(max(value(1), 0), 1);
end
end

function surfaceDepth = buildPitSurfaceDepth(X, Y, vPits)
surfaceDepth = zeros(size(X));
for pitIdx = 1:numel(vPits)
    pit = vPits(pitIdx);
    rho = hexagonalRadius(X, Y, pit);
    localDepth = pit.depthNm * max(0, 1 - rho);
    localDepth(rho > 1) = 0;
    surfaceDepth = max(surfaceDepth, localDepth);
end
end

function fig = renderCrossSectionFigure(result, visibility)
fig = figure('Color', 'w', 'Visible', visibility, ...
    'Position', [110, 80, 1180, 680], 'Name', 'Point-cloud V-pit and MQW cross section');
plotPosition = [0.08, 0.11, 0.70, 0.80];
ax = axes(fig, 'Position', plotPosition);

if isempty(result.vPits)
    targetY = mean(result.yGridNm);
else
    targetY = result.vPits(1).centerNm(2);
end
[~, yIdx] = min(abs(result.yGridNm - targetY));
section = 10 * squeeze(result.displayVolume(yIdx, :, :))';
imagesc(ax, result.xGridNm, result.depthGridNm, section);
set(ax, 'YDir', 'reverse');
colormap(ax, turbo(256));
clim(ax, [0, 10]);
hold(ax, 'on');

cb = colorbar(ax);
cb.Position = [0.84, 0.19, 0.024, 0.64];
cb.Label.String = 'Normalized point-cloud intensity';
ax.Position = plotPosition;
xlabel(ax, 'X (nm)');
ylabel(ax, 'Depth below surface (nm)');
title(ax, sprintf('V-pit and MQW cross section at Y = %.1f nm', result.yGridNm(yIdx)));
axis(ax, 'tight');
end
