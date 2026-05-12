function varargout = vpit(action, varargin)
switch lower(string(action))
    case "getdefaults"
        varargout{1} = getVpitsDefaults();
    case "buildgeometries"
        varargout{1} = buildVpits(varargin{:});
    case "resolvepointstate"
        varargout{1} = resolvePointState(varargin{:});
    case "firstboundarycrossing"
        varargout{1} = firstBoundaryCrossing(varargin{:});
    case "projectcarriersource"
        varargout{1} = projectCarrierSource(varargin{:});
    case "calculatedepthfromcomposition"
        varargout{1} = calculateVPitDepthFromComposition(varargin{:});
    otherwise
        error('precies:vpit:InvalidAction', 'Unsupported action: %s', action);
end
end
function depth = calculateVPitDepthFromComposition(topRadius, layers, layerBoundaries)
  
    surfaceZ = layerBoundaries(end-1) * 1e9;
    maxDepth = min(500, surfaceZ - layerBoundaries(1)*1e9);
    numSamples = 100;
    depthSamples = linspace(0, maxDepth, numSamples);
    
    totalInComposition = 0;
    validSamples = 0;
    
    for i = 1:numSamples
        currentDepth = depthSamples(i);
        currentZ = surfaceZ - currentDepth;
        layerIdx = find(currentZ >= layerBoundaries(1:end-1)*1e9 & currentZ < layerBoundaries(2:end)*1e9, 1);
        
        if ~isempty(layerIdx) && layerIdx <= size(layers, 1)
            layerInComp = layers{layerIdx, 2};
            weight = exp(-currentDepth / (maxDepth * 0.3));
            totalInComposition = totalInComposition + layerInComp * weight;
            validSamples = validSamples + weight;
        end
    end
    
    if validSamples > 0
        avgInComposition = totalInComposition / validSamples;
    else
        mqwLayerIdx = find(contains(layers(:, 1), 'MQW-Well'), 1);
        if ~isempty(mqwLayerIdx)
            avgInComposition = layers{mqwLayerIdx, 2};
        else
            avgInComposition = 0.15;
        end
    end
    materialSystem = struct();
    materialSystem.latticeConstant_a_GaN = 0.3189; 
    materialSystem.latticeConstant_c_GaN = 0.5185;  
    materialSystem.latticeConstant_a_InN = 0.3545;
    materialSystem.latticeConstant_c_InN = 0.5703;

    if avgInComposition < 0.15
        materialSystem.facetType = '10-11'; 
    elseif avgInComposition < 0.25  
        materialSystem.facetType = '11-22'; 
    else
        materialSystem.facetType = '10-11'; 
    end

    a_alloy = avgInComposition * materialSystem.latticeConstant_a_InN + ...
              (1 - avgInComposition) * materialSystem.latticeConstant_a_GaN;
    c_alloy = avgInComposition * materialSystem.latticeConstant_c_InN + ...
              (1 - avgInComposition) * materialSystem.latticeConstant_c_GaN;

    switch materialSystem.facetType
        case '10-11'
            facetAngle = atand(c_alloy / (sqrt(3) * a_alloy));
        case '11-22' 
            facetAngle = atand(c_alloy / (sqrt(8/3) * a_alloy));
        case '10-12'
            facetAngle = atand(2*c_alloy / (sqrt(3)*a_alloy));
        case '10-13'
            facetAngle = atand(3*c_alloy / (sqrt(3)*a_alloy));
        otherwise
            facetAngle = atand(c_alloy / (sqrt(3) * a_alloy));
    end

    depth = topRadius / tand(facetAngle);

    maxReasonableDepth = 200; 
    depth = min(depth, maxReasonableDepth);
end

function defaults = getVpitsDefaults()
    defaults = struct( ...
        'facetType', '10-11', ...
        'semipolarInScale', 0.85, ...
        'semipolarThicknessScale', 0.70, ...
        'carrierCaptureDistance_nm', 100, ...
        'carrierLateralDiffusionLength_nm', 800, ...
        'displayShellAlpha', 0.28, ...
        'displayFacetAlpha', 0.16, ...
        'defaultFacetAngleDeg', 60, ...
        'defaultAirAbsorption', 0, ...
        'stateEpsilon_nm', 1e-3);
end

function builtVpits = buildVpits(rawVpits, layers, layerBoundaries)
    defaults = getVpitsDefaults();
    if nargin < 1 || isempty(rawVpits)
        builtVpits = cell(0, 1);
        return;
    end
    if nargin < 2
        layers = {};
    end
    if nargin < 3
        layerBoundaries = [];
    end

    builtVpits = rawVpits;
    avgInComposition = estimateActiveInComposition(layers, layerBoundaries);
    facetAngleDeg = calculateFacetAngleForComposition(avgInComposition, defaults.facetType, defaults.defaultFacetAngleDeg);
    surfaceZ_nm = inferSurfaceZ(layerBoundaries, rawVpits);

    for i = 1:numel(rawVpits)
        vpit = rawVpits{i};
        if isempty(vpit)
            continue;
        end

        if ~isfield(vpit, 'depth') || isempty(vpit.depth) || ~isfinite(vpit.depth) || vpit.depth <= 0
            if ~isempty(layers) && ~isempty(layerBoundaries)
                vpit.depth = calculateVPitDepthFromComposition(vpit.topRadius, layers, layerBoundaries);
            else
                vpit.depth = max(vpit.topRadius / tand(facetAngleDeg), 1);
            end
        end

        localFacetAngleDeg = facetAngleDeg;
        if isfinite(vpit.depth) && vpit.depth > 0
            localFacetAngleDeg = atand(vpit.topRadius / vpit.depth);
        end
        if ~isfield(vpit, 'orientationDeg') || isempty(vpit.orientationDeg) || ~isfinite(vpit.orientationDeg)
            vpit.orientationDeg = 0;
        end
        vpit.orientationDeg = mod(vpit.orientationDeg, 60);

        shellThicknessByLayer_nm = containers.Map('KeyType', 'double', 'ValueType', 'double');
        if ~isempty(layers)
            for layerIdx = 1:size(layers, 1)
                layerName = layers{layerIdx, 1};
                if contains(layerName, 'MQW-Well') || contains(layerName, 'MQW-Barrier')
                    shellThicknessByLayer_nm(layerIdx) = max(0.5, layers{layerIdx, 3} * defaults.semipolarThicknessScale);
                end
            end
        end

        apothemToVertex = 1 / cosd(30);
        vertexAnglesDeg = vpit.orientationDeg + (30:60:390);
        vertexRadius_nm = vpit.topRadius * apothemToVertex;
        topVertices = [ ...
            vpit.center(1) + vertexRadius_nm * cosd(vertexAnglesDeg(:)), ...
            vpit.center(2) + vertexRadius_nm * sind(vertexAnglesDeg(:)), ...
            surfaceZ_nm * ones(numel(vertexAnglesDeg), 1)];
        apex = [vpit.center(1), vpit.center(2), surfaceZ_nm - vpit.depth];

        facets = repmat(struct( ...
            'index', 0, ...
            'normalXY', [0, 0], ...
            'planeNormalOutward', [0, 0, 1], ...
            'planeOffset', 0, ...
            'planeNorm', 1, ...
            'vertexStart', [0, 0, 0], ...
            'vertexEnd', [0, 0, 0]), 6, 1);

        slope = vpit.topRadius / max(vpit.depth, eps);
        for facetIdx = 1:6
            angleDeg = vpit.orientationDeg + (facetIdx - 1) * 60;
            normalXY = [cosd(angleDeg), sind(angleDeg)];
            planeVector = [normalXY, -slope];
            planeNorm = norm(planeVector);
            planeOffset = (slope * surfaceZ_nm - vpit.topRadius - dot(normalXY, vpit.center)) / planeNorm;

            facets(facetIdx).index = facetIdx;
            facets(facetIdx).normalXY = normalXY;
            facets(facetIdx).planeNormalOutward = planeVector / planeNorm;
            facets(facetIdx).planeOffset = planeOffset;
            facets(facetIdx).planeNorm = planeNorm;
            facets(facetIdx).vertexStart = topVertices(facetIdx, :);
            facets(facetIdx).vertexEnd = topVertices(facetIdx + 1, :);
        end

        vpit.modelVersion = 2;
        vpit.facetType = defaults.facetType;
        vpit.avgInComposition = avgInComposition;
        vpit.facetAngleDeg = localFacetAngleDeg;
        vpit.surfaceZ_nm = surfaceZ_nm;
        vpit.apex = apex;
        vpit.topVertices = topVertices;
        vpit.facets = facets;
        vpit.vertexAnglesDeg = vertexAnglesDeg;
        vpit.semipolarInScale = defaults.semipolarInScale;
        vpit.semipolarThicknessScale = defaults.semipolarThicknessScale;
        vpit.carrierCaptureDistance_nm = defaults.carrierCaptureDistance_nm;
        vpit.carrierLateralDiffusionLength_nm = defaults.carrierLateralDiffusionLength_nm;
        vpit.shellThicknessByLayer_nm = shellThicknessByLayer_nm;
        builtVpits{i} = vpit;
    end
end

function state = resolvePointState(pos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ)
    defaults = getVpitsDefaults();
    if nargin < 5 || isempty(nominalSurfaceZ)
        nominalSurfaceZ = inferSurfaceZ(layerBoundaries, Vpits);
    end
    Vpits = ensureBuiltVpits(Vpits, layers, layerBoundaries);

    state = buildBaseState(pos_nm, layers, layerBoundaries);
    state.nominalSurfaceZ_nm = nominalSurfaceZ;
    state.vpitIndex = 0;
    state.facetIndex = 0;
    state.closestFacetDistance_nm = inf;
    state.closestFacetNormalXY = [0, 0];
    state.closestFacetNormal3D = [0, 0, 1];
    state.closestFacetPlaneNorm = 1;
    state.layerShellThickness_nm = 0;

    if isempty(Vpits) || isempty(layers) || state.layerIndex < 1 || state.layerIndex > size(layers, 1)
        if strcmp(state.materialKind, 'air')
            state.nFunc = @(lambda) 1.0;
            state.alphaFunc = @(lambda) defaults.defaultAirAbsorption;
        end
        return;
    end

    layerName = layers{state.layerIndex, 1};
    isActiveMQWLayer = contains(layerName, 'MQW-Well') || contains(layerName, 'MQW-Barrier');

    for vpitIdx = 1:numel(Vpits)
        vpit = Vpits{vpitIdx};
        [insidePit, facetIndex, inwardDistance_nm, facet] = classifyPointAgainstVpit(pos_nm, vpit);
        if ~insidePit
            positiveDistance_nm = nearestPositiveFacetDistance(pos_nm, vpit);
            if positiveDistance_nm < state.closestFacetDistance_nm
                state.closestFacetDistance_nm = positiveDistance_nm;
                state.closestFacetNormalXY = facet.normalXY;
                state.closestFacetNormal3D = facet.planeNormalOutward;
                state.closestFacetPlaneNorm = facet.planeNorm;
                state.vpitIndex = vpitIdx;
                state.facetIndex = facetIndex;
            end
            continue;
        end

        state.vpitIndex = vpitIdx;
        state.facetIndex = facetIndex;
        state.closestFacetDistance_nm = inwardDistance_nm;
        state.closestFacetNormalXY = facet.normalXY;
        state.closestFacetNormal3D = facet.planeNormalOutward;
        state.closestFacetPlaneNorm = facet.planeNorm;

        if isActiveMQWLayer && isKey(vpit.shellThicknessByLayer_nm, state.layerIndex)
            shellThickness_nm = vpit.shellThicknessByLayer_nm(state.layerIndex);
            state.layerShellThickness_nm = shellThickness_nm;
            if inwardDistance_nm <= shellThickness_nm
                state.materialKind = 'semipolar_shell';
                state.regionType = 'semipolar_shell';
                state.inComposition = state.inComposition * defaults.semipolarInScale;
                state.layerThickness_nm = layers{state.layerIndex, 3} * defaults.semipolarThicknessScale;
                state.bandgap_eV = calculateBandgapForComposition(state.inComposition);
                state.nFunc = @(lambda) localCalculateInGaNRefractiveIndex(lambda, state.inComposition);
                state.alphaFunc = @(lambda) localCalculateInGaNAbsorption(lambda, state.bandgap_eV);
            else
                state.materialKind = 'cavity';
                state.regionType = 'cavity';
                state.inComposition = 0;
                state.bandgap_eV = NaN;
                state.nFunc = @(lambda) 1.0;
                state.alphaFunc = @(lambda) defaults.defaultAirAbsorption;
            end
        else
            state.materialKind = 'cavity';
            state.regionType = 'cavity';
            state.inComposition = 0;
            state.bandgap_eV = NaN;
            state.nFunc = @(lambda) 1.0;
            state.alphaFunc = @(lambda) defaults.defaultAirAbsorption;
        end
        break;
    end
end

function crossing = firstBoundaryCrossing(oldPos_m, newPos_m, layers, layerBoundaries, Vpits, nominalSurfaceZ)
    defaults = getVpitsDefaults();
    oldPos_nm = oldPos_m * 1e9;
    newPos_nm = newPos_m * 1e9;
    direction_nm = newPos_nm - oldPos_nm;

    crossing = struct( ...
        'crossed', false, ...
        'point_m', oldPos_m, ...
        'normal', [0, 0, 0], ...
        'boundaryType', "", ...
        'oldState', resolvePointState(oldPos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ), ...
        'newState', resolvePointState(newPos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ), ...
        'newLayer', 0, ...
        't', 1);
    crossing.newLayer = crossing.newState.layerIndex;

    if statesEquivalent(crossing.oldState, crossing.newState)
        return;
    end

    directionNorm = norm(direction_nm);
    if directionNorm <= eps
        crossing.crossed = true;
        crossing.normal = [0, 0, 1];
        return;
    end

    directionUnit = direction_nm / directionNorm;
    candidateTs = [];
    candidateTypes = strings(0, 1);
    candidateNormals = zeros(0, 3);

    if oldPos_nm(3) ~= newPos_nm(3) && ~isempty(layerBoundaries)
        planeBoundaries_nm = layerBoundaries(2:end-1) * 1e9;
        for planeIdx = 1:numel(planeBoundaries_nm)
            t = (planeBoundaries_nm(planeIdx) - oldPos_nm(3)) / (newPos_nm(3) - oldPos_nm(3));
            if t <= 0 || t > 1
                continue;
            end
            point_nm = oldPos_nm + t * direction_nm;
            oldState = resolvePointState(point_nm - directionUnit * defaults.stateEpsilon_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
            newState = resolvePointState(point_nm + directionUnit * defaults.stateEpsilon_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
            if statesEquivalent(oldState, newState)
                continue;
            end
            candidateTs(end + 1, 1) = t;
            candidateTypes(end + 1, 1) = "layer_plane";
            candidateNormals(end + 1, :) = [0, 0, sign(newPos_nm(3) - oldPos_nm(3))];
        end
    end

    Vpits = ensureBuiltVpits(Vpits, layers, layerBoundaries);
    for vpitIdx = 1:numel(Vpits)
        vpit = Vpits{vpitIdx};
        for facetIdx = 1:numel(vpit.facets)
            facet = vpit.facets(facetIdx);
            planeNormal = facet.planeNormalOutward;
            planeOffset = facet.planeOffset;
            denom = dot(planeNormal, direction_nm);
            if abs(denom) > eps
                tFacet = -(dot(planeNormal, oldPos_nm) + planeOffset) / denom;
                if tFacet > 0 && tFacet <= 1
                    point_nm = oldPos_nm + tFacet * direction_nm;
                    if isPointWithinFacetDepth(point_nm, vpit)
                        oldState = resolvePointState(point_nm - directionUnit * defaults.stateEpsilon_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
                        newState = resolvePointState(point_nm + directionUnit * defaults.stateEpsilon_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
                        if ~statesEquivalent(oldState, newState)
                            candidateTs(end + 1, 1) = tFacet;
                            candidateTypes(end + 1, 1) = "facet_surface";
                            candidateNormals(end + 1, :) = orientFacetNormal(facet.planeNormalOutward, oldState, newState);
                        end
                    end
                end
            end

            shellThickness_nm = 0;
            if crossing.oldState.layerIndex > 0 && isKey(vpit.shellThicknessByLayer_nm, crossing.oldState.layerIndex)
                shellThickness_nm = vpit.shellThicknessByLayer_nm(crossing.oldState.layerIndex);
            elseif crossing.newState.layerIndex > 0 && isKey(vpit.shellThicknessByLayer_nm, crossing.newState.layerIndex)
                shellThickness_nm = vpit.shellThicknessByLayer_nm(crossing.newState.layerIndex);
            end

            if shellThickness_nm <= 0
                continue;
            end

            innerOffset = planeOffset + shellThickness_nm;
            denomInner = dot(planeNormal, direction_nm);
            if abs(denomInner) <= eps
                continue;
            end

            tInner = -(dot(planeNormal, oldPos_nm) + innerOffset) / denomInner;
            if tInner <= 0 || tInner > 1
                continue;
            end

            point_nm = oldPos_nm + tInner * direction_nm;
            if ~isPointWithinFacetDepth(point_nm, vpit)
                continue;
            end

            oldState = resolvePointState(point_nm - directionUnit * defaults.stateEpsilon_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
            newState = resolvePointState(point_nm + directionUnit * defaults.stateEpsilon_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
            if statesEquivalent(oldState, newState)
                continue;
            end
            candidateTs(end + 1, 1) = tInner;
            candidateTypes(end + 1, 1) = "shell_inner";
            candidateNormals(end + 1, :) = orientFacetNormal(facet.planeNormalOutward, oldState, newState);
        end
    end

    if isempty(candidateTs)
        crossing.crossed = true;
        crossing.normal = [0, 0, sign(newPos_nm(3) - oldPos_nm(3))];
        return;
    end

    [crossing.t, minIdx] = min(candidateTs);
    crossing.crossed = true;
    crossing.point_m = (oldPos_nm + crossing.t * direction_nm) * 1e-9;
    crossing.boundaryType = candidateTypes(minIdx);
    crossing.normal = candidateNormals(minIdx, :);
    crossing.oldState = resolvePointState((crossing.point_m * 1e9) - directionUnit * defaults.stateEpsilon_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
    crossing.newState = resolvePointState((crossing.point_m * 1e9) + directionUnit * defaults.stateEpsilon_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
    crossing.newLayer = crossing.newState.layerIndex;
end

function projectedPos_nm = projectCarrierSource(pos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ, sourceState)
    defaults = getVpitsDefaults();
    if nargin < 6 || isempty(sourceState)
        sourceState = resolvePointState(pos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
    end

    projectedPos_nm = pos_nm;
    Vpits = ensureBuiltVpits(Vpits, layers, layerBoundaries);
    if isempty(Vpits)
        return;
    end

    layerIdx = sourceState.layerIndex;
    if layerIdx < 1 || layerIdx > size(layers, 1)
        return;
    end

    if sourceState.vpitIndex < 1 || sourceState.facetIndex < 1
        captureInfo = findNearestCaptureFacet(pos_nm, layerIdx, Vpits, nominalSurfaceZ);
        if ~captureInfo.isValid
            return;
        end
        vpit = Vpits{captureInfo.vpitIndex};
        facet = vpit.facets(captureInfo.facetIndex);
    else
        vpit = Vpits{sourceState.vpitIndex};
        facet = vpit.facets(sourceState.facetIndex);
    end

    layerCenter_nm = mean(layerBoundaries(layerIdx:layerIdx + 1)) * 1e9;
    apothem_nm = currentApothemAtZ(vpit, layerCenter_nm);
    normalXY = facet.normalXY;
    tangentXY = [-normalXY(2), normalXY(1)];
    anchorXY = vpit.center + normalXY * (apothem_nm + 5);

    outwardShift_nm = exprnd(defaults.carrierLateralDiffusionLength_nm / 3);
    tangentialShift_nm = randn() * defaults.carrierLateralDiffusionLength_nm / 5;
    projectedXY = anchorXY + normalXY * outwardShift_nm + tangentXY * tangentialShift_nm;
    projectedPos_nm = [projectedXY, layerCenter_nm];

    for iter = 1:20
        projectedState = resolvePointState(projectedPos_nm, layers, layerBoundaries, Vpits, nominalSurfaceZ);
        if strcmp(projectedState.materialKind, 'bulk')
            return;
        end
        projectedPos_nm(1:2) = projectedPos_nm(1:2) + normalXY * 15;
    end
end

function state = buildBaseState(pos_nm, layers, layerBoundaries)
    state = struct( ...
        'layerIndex', 0, ...
        'layerName', 'Air', ...
        'materialKind', 'air', ...
        'regionType', 'air', ...
        'inComposition', 0, ...
        'layerThickness_nm', 0, ...
        'bandgap_eV', NaN, ...
        'isMQW', false, ...
        'nFunc', [], ...
        'alphaFunc', []);

    if isempty(layerBoundaries) || isempty(layers)
        return;
    end

    surfaceZ_nm = layerBoundaries(end - 1) * 1e9;
    bottomZ_nm = layerBoundaries(1) * 1e9;
    if pos_nm(3) > surfaceZ_nm || pos_nm(3) < bottomZ_nm
        return;
    end

    layerIdx = find(pos_nm(3) <= layerBoundaries(2:end - 1) * 1e9, 1);
    if isempty(layerIdx)
        layerIdx = size(layers, 1);
    end

    state.layerIndex = layerIdx;
    state.layerName = layers{layerIdx, 1};
    state.materialKind = 'bulk';
    state.regionType = 'bulk';
    state.inComposition = layers{layerIdx, 2};
    state.layerThickness_nm = layers{layerIdx, 3};
    state.bandgap_eV = calculateBandgapForComposition(state.inComposition);
    state.isMQW = contains(state.layerName, 'MQW');

    if size(layers, 2) >= 6 && ~isempty(layers{layerIdx, 5})
        state.nFunc = layers{layerIdx, 5};
        state.alphaFunc = layers{layerIdx, 6};
    end
end

function avgInComposition = estimateActiveInComposition(layers, ~)
    if nargin < 1 || isempty(layers)
        avgInComposition = 0.15;
        return;
    end
    mqwMask = contains(string(layers(:, 1)), "MQW-Well");
    if any(mqwMask)
        avgInComposition = mean(cell2mat(layers(mqwMask, 2)));
    else
        avgInComposition = 0.15;
    end
end

function facetAngleDeg = calculateFacetAngleForComposition(avgInComposition, facetType, defaultFacetAngleDeg)
    aGaN = 0.3189;
    cGaN = 0.5185;
    aInN = 0.3545;
    cInN = 0.5703;

    aAlloy = avgInComposition * aInN + (1 - avgInComposition) * aGaN;
    cAlloy = avgInComposition * cInN + (1 - avgInComposition) * cGaN;

    switch facetType
        case '10-11'
            facetAngleDeg = atand(cAlloy / (sqrt(3) * aAlloy));
        case '11-22'
            facetAngleDeg = atand(cAlloy / (sqrt(8 / 3) * aAlloy));
        otherwise
            facetAngleDeg = defaultFacetAngleDeg;
    end
end

function surfaceZ_nm = inferSurfaceZ(layerBoundaries, Vpits)
    if nargin >= 1 && ~isempty(layerBoundaries)
        surfaceZ_nm = layerBoundaries(end - 1) * 1e9;
        return;
    end
    if nargin >= 2 && ~isempty(Vpits) && isfield(Vpits{1}, 'surfaceZ_nm')
        surfaceZ_nm = Vpits{1}.surfaceZ_nm;
        return;
    end
    surfaceZ_nm = 0;
end

function Vpits = ensureBuiltVpits(Vpits, layers, layerBoundaries)
    if nargin < 1 || isempty(Vpits)
        Vpits = cell(0, 1);
        return;
    end
    if ~iscell(Vpits)
        Vpits = {Vpits};
    end
    if isempty(Vpits)
        return;
    end
    if ~isfield(Vpits{1}, 'facets')
        Vpits = buildVpits(Vpits, layers, layerBoundaries);
    end
end

function [insidePit, facetIndex, inwardDistance_nm, facet] = classifyPointAgainstVpit(pos_nm, vpit)
    insidePit = false;
    facetIndex = 1;
    inwardDistance_nm = inf;
    facet = vpit.facets(1);

    depthBelowSurface_nm = vpit.surfaceZ_nm - pos_nm(3);
    if depthBelowSurface_nm < 0 || depthBelowSurface_nm > vpit.depth
        return;
    end

    eqValues = zeros(numel(vpit.facets), 1);
    inwardDistances = zeros(numel(vpit.facets), 1);
    for idx = 1:numel(vpit.facets)
        eqValues(idx) = evaluateFacetPlane(vpit.facets(idx), pos_nm);
            inwardDistances(idx) = -eqValues(idx);
    end

    insidePit = all(eqValues <= 1e-9);
    [inwardDistance_nm, facetIndex] = min(inwardDistances);
    facet = vpit.facets(facetIndex);
end

function positiveDistance_nm = nearestPositiveFacetDistance(pos_nm, vpit)
    positiveDistance_nm = inf;
    for idx = 1:numel(vpit.facets)
        eqValue = evaluateFacetPlane(vpit.facets(idx), pos_nm);
        if eqValue > 0
            positiveDistance_nm = min(positiveDistance_nm, eqValue);
        end
    end
end

function eqValue = evaluateFacetPlane(facet, pos_nm)
    eqValue = dot(facet.planeNormalOutward, pos_nm) + facet.planeOffset;
end

function isInside = isPointWithinFacetDepth(pos_nm, vpit)
    depthBelowSurface_nm = vpit.surfaceZ_nm - pos_nm(3);
    isInside = depthBelowSurface_nm >= -1e-6 && depthBelowSurface_nm <= vpit.depth + 1e-6;
end

function normal = orientFacetNormal(baseNormal, oldState, newState)
    if stateRank(newState) >= stateRank(oldState)
        normal = baseNormal;
    else
        normal = -baseNormal;
    end
end

function rank = stateRank(state)
    switch state.materialKind
        case 'cavity'
            rank = 0;
        case 'semipolar_shell'
            rank = 1;
        case 'bulk'
            rank = 2;
        otherwise
            rank = 3;
    end
end

function tf = statesEquivalent(stateA, stateB)
    tf = stateA.layerIndex == stateB.layerIndex && ...
         strcmp(stateA.materialKind, stateB.materialKind) && ...
         stateA.vpitIndex == stateB.vpitIndex && ...
         stateA.facetIndex == stateB.facetIndex;
end

function captureInfo = findNearestCaptureFacet(pos_nm, ~, Vpits, nominalSurfaceZ)
    captureInfo = struct('isValid', false, 'distance_nm', inf, 'vpitIndex', 0, 'facetIndex', 0);
    if isempty(Vpits)
        return;
    end

    for vpitIdx = 1:numel(Vpits)
        vpit = Vpits{vpitIdx};
        depthBelowSurface_nm = nominalSurfaceZ - pos_nm(3);
        if depthBelowSurface_nm < 0 || depthBelowSurface_nm > vpit.depth
            continue;
        end

        for facetIdx = 1:numel(vpit.facets)
            eqValue = evaluateFacetPlane(vpit.facets(facetIdx), pos_nm);
            if eqValue <= 0
                continue;
            end
            distance_nm = eqValue / vpit.facets(facetIdx).planeNorm;
            if distance_nm < captureInfo.distance_nm
                captureInfo.isValid = true;
                captureInfo.distance_nm = distance_nm;
                captureInfo.vpitIndex = vpitIdx;
                captureInfo.facetIndex = facetIdx;
            end
        end
    end
end

function apothem_nm = currentApothemAtZ(vpit, z_nm)
    depthBelowSurface_nm = vpit.surfaceZ_nm - z_nm;
    apothem_nm = vpit.topRadius * max(0, 1 - depthBelowSurface_nm / max(vpit.depth, eps));
end

function Eg = calculateBandgapForComposition(inComposition)
    EgGaN = 3.3032;
    EgInN = 0.6086;
    Eg = EgInN * inComposition + EgGaN * (1 - inComposition) - 1.43 * inComposition * (1 - inComposition);
end

function [layers, layerBoundaries] = localBuildLayers(layer_params, params)
    mqw_pairs = max(1, round(params.mqw_pairs));
    mqw_barriers = localResolveMqwBarrierCount(params, mqw_pairs);
    mqwTemplateCount = nnz(strcmp(layer_params(:, 1), 'MQW-Barrier'));
    maxLayerRows = size(layer_params, 1) + mqwTemplateCount * max(0, mqw_pairs + mqw_barriers);
    layers = cell(maxLayerRows, size(layer_params, 2));
    layerWriteIdx = 0;
    skip_mqw_well = false;
    for i = 1:size(layer_params, 1)
        if skip_mqw_well
            skip_mqw_well = false;
            continue;
        end

        if strcmp(layer_params{i, 1}, 'MQW-Barrier')
            barrierRow = layer_params(i, :);
            wellRow = layer_params(i + 1, :);
            if mqw_barriers >= mqw_pairs
                layerWriteIdx = layerWriteIdx + 1;
                layers(layerWriteIdx, :) = barrierRow;
                remainingBarriers = mqw_barriers - 1;
                for j = 1:mqw_pairs
                    layerWriteIdx = layerWriteIdx + 1;
                    layers(layerWriteIdx, :) = wellRow;
                    if remainingBarriers > 0
                        layerWriteIdx = layerWriteIdx + 1;
                        layers(layerWriteIdx, :) = barrierRow;
                        remainingBarriers = remainingBarriers - 1;
                    end
                end
            else
                layerWriteIdx = layerWriteIdx + 1;
                layers(layerWriteIdx, :) = wellRow;
                remainingWells = mqw_pairs - 1;
                remainingBarriers = mqw_barriers;
                while remainingWells > 0 || remainingBarriers > 0
                    if remainingBarriers > 0
                        layerWriteIdx = layerWriteIdx + 1;
                        layers(layerWriteIdx, :) = barrierRow;
                        remainingBarriers = remainingBarriers - 1;
                    end
                    if remainingWells > 0
                        layerWriteIdx = layerWriteIdx + 1;
                        layers(layerWriteIdx, :) = wellRow;
                        remainingWells = remainingWells - 1;
                    end
                end
            end
            skip_mqw_well = true;
        else
            layerWriteIdx = layerWriteIdx + 1;
            layers(layerWriteIdx, :) = layer_params(i, :);
        end
    end
    layers = layers(1:layerWriteIdx, :);

    layerBoundaries = zeros(1, size(layers, 1) + 1);
    layerBoundaries(1) = -layers{1, 3} * 1e-9;
    for i = 1:size(layers, 1)
        thickness = layers{i, 3} * 1e-9;
        layerBoundaries(i + 1) = layerBoundaries(i) + thickness;
    end
    layerBoundaries(end + 1) = Inf;
end

function mqw_barriers = localResolveMqwBarrierCount(params, mqw_pairs)
    if isstruct(params) && isfield(params, 'mqw_barriers') && ~isempty(params.mqw_barriers)
        mqw_barriers = params.mqw_barriers;
    elseif isstruct(params) && isfield(params, 'mqwBarrierCount') && ~isempty(params.mqwBarrierCount)
        mqw_barriers = params.mqwBarrierCount;
    else
        mqw_barriers = mqw_pairs - 1;
    end
    mqw_barriers = max(0, round(double(mqw_barriers)));
end

function n_InGaN = localCalculateInGaNRefractiveIndex(lambda, inComposition)
    EgGaN = 3.3032;
    EgInN = 0.6086;
    lambda_nm = lambda * 1e9;
    E_g_InGaN = EgGaN * (1 - inComposition) + EgInN * inComposition - 1.4 * inComposition * (1 - inComposition);
    E_photon = 1240 ./ lambda_nm;
    delta_E = E_g_InGaN - EgGaN;
    E_shifted = E_photon - delta_E;

    shifted_lambda_nm = zeros(size(lambda_nm));
    valid_E = E_shifted ~= 0;
    shifted_lambda_nm(valid_E) = 1240 ./ E_shifted(valid_E);
    shifted_lambda_nm(~valid_E) = Inf;
    n_InGaN = localCalculateGaNRefractiveIndex(shifted_lambda_nm * 1e-9);
end

function n = localCalculateGaNRefractiveIndex(lambda)
    EgGaN = 3.3032;
    EgInN = 0.6086;
    y = 0;
    lambda_nm = lambda * 1e9;
    h_nu_eV = 1240 ./ lambda_nm;
    Eg_eV = EgInN * y + EgGaN * (1 - y) - 1.4 * y * (1 - y);
    ratio = h_nu_eV ./ Eg_eV;

    a = 9.82661 - 8.21608 * y - 31.5902 * y^2;
    b = 2.73591 + 0.84249 * y - 6.29321 * y^2;
    term = a * (ratio) .^ (-2) .* (2 - sqrt(1 + ratio) - sqrt(1 - ratio));
    term(ratio >= 1) = a * (ratio(ratio >= 1)) .^ (-2) .* (2 - sqrt(1 + ratio(ratio >= 1)));
    n = sqrt(term + b);
end

function alpha = localCalculateInGaNAbsorption(lambda, bandgap)
    lambda_nm = lambda * 1e9;
    E = 1240 ./ lambda_nm;
    E_B = 0.993 + 0.719 * bandgap;
    delta_E = 0.06;
    alpha_i = 5e3;
    alpha_0 = 5e6;
    exponent = (E_B - E) / delta_E;
    alpha = alpha_i + (alpha_0 - alpha_i) ./ (1 + exp(exponent));
end
