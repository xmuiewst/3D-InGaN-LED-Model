function varargout = tiffops(action, varargin)
switch lower(string(action))
    case "analyzetiffstructure"
        [varargout{1:nargout}] = analyzeTiffStructure(varargin{:});
    case "parseomexmlforwavelength"
        varargout{1} = parseOmeXmlForWavelength(varargin{:});
    case "tryalternativewavelengthsources"
        varargout{1} = tryAlternativeWavelengthSources(varargin{:});
    otherwise
        error('precies:tiffops:InvalidAction', 'Unsupported action: %s', action);
end
end
function [imageTypes, layerDimensions] = analyzeTiffStructure(tiffInfo)
    numLayers = length(tiffInfo);
    imageTypes = cell(numLayers, 1);
    layerDimensions = cell(numLayers, 1);
    
    for i = 1:numLayers
        if isfield(tiffInfo(i), 'PageName') && ~isempty(tiffInfo(i).PageName)
            imageTypes{i} = tiffInfo(i).PageName;
        elseif isfield(tiffInfo(i), 'ImageDescription') && ~isempty(tiffInfo(i).ImageDescription)
            desc = tiffInfo(i).ImageDescription;
            nameMatch = regexp(desc, '<Image[^>]*Name="([^"]+)"', 'tokens');
            if ~isempty(nameMatch)
                imageTypes{i} = nameMatch{1}{1};
            else
                imageTypes{i} = sprintf('Layer_%d', i);
            end
        else
            imageTypes{i} = sprintf('Layer_%d', i);
        end

        layerDimensions{i} = [tiffInfo(i).Width, tiffInfo(i).Height];
    end
end

function wavelengthAxis = parseOmeXmlForWavelength(omeXML)
    wavelengthAxis = [];
    
    if isempty(omeXML)
        warning('OME-XML is empty');
        return;
    end

    rangeTokens = regexp(omeXML, 'SpectrumRange="([0-9.]+)\s*-\s*([0-9.]+)\s*nm"', 'tokens', 'once');
    if ~isempty(rangeTokens)
        startWavelength = str2double(rangeTokens{1});
        endWavelength = str2double(rangeTokens{2});
        if isfinite(startWavelength) && isfinite(endWavelength) && startWavelength < endWavelength
            wavelengthAxis = linspace(startWavelength, endWavelength, 1024)';
            return;
        end
    end
    
    try
        omeStruct = xml2struct(omeXML);

        if isfield(omeStruct, 'OME') && isfield(omeStruct.OME, 'Instrument')
            instruments = omeStruct.OME.Instrument;
            if iscell(instruments)
                instruments = instruments{1};
            end
            
            if isfield(instruments, 'Spectrograph')
                spectrograph = instruments.Spectrograph;
                if isfield(spectrograph, 'Attributes') && isfield(spectrograph.Attributes, 'wavelength')
                    centerWavelength = str2double(spectrograph.Attributes.wavelength) * 1e9;
                    wavelengthAxis = linspace(centerWavelength-200, centerWavelength+200, 1024)';
                    return;
                end
            end
            
            if isfield(omeStruct.OME, 'Image')
                images = omeStruct.OME.Image;
                if iscell(images)
                    images = images{1};
                end
                
                if isfield(images, 'Pixels') && isfield(images.Pixels, 'Channel')
                    channels = images.Pixels.Channel;
                    if iscell(channels)
                        channels = channels{1};
                    end
                    
                    if isfield(channels, 'Attributes') && isfield(channels.Attributes, 'EmissionWavelength')
                        centerWavelength = str2double(channels.Attributes.EmissionWavelength);
                        wavelengthAxis = linspace(centerWavelength-200, centerWavelength+200, 1024)';
                        return;
                    end
                end
            end
        end
    catch
    end
    
    patterns = {
        '"Spectrograph":\s*{.*?"wavelength":\s*([0-9.e+-]+)';
        '<Channel.*?EmissionWavelength="([^"]+)"';
        '([Ww]avelength)[^0-9]*([0-9.e+-]+)';
        '([0-9]+)\s*[nm]';
        'center.*?([0-9]+)\s*nm';
        'lambda.*?([0-9]+)'
    };
    
    for i = 1:length(patterns)
        tokens = regexp(omeXML, patterns{i}, 'tokens');
        if ~isempty(tokens)
            try
                if i == 1 || i == 2
                    wavelengthValue = str2double(tokens{1}{1});
                    if i == 1
                        wavelengthValue = wavelengthValue * 1e9; 
                    end
                else
                    wavelengthValue = str2double(tokens{1}{end});
                end
                
                if ~isnan(wavelengthValue) && wavelengthValue > 0
                    wavelengthAxis = linspace(wavelengthValue-200, wavelengthValue+200, 1024)';
                    return;
                end
            catch
                continue;
            end
        end
    end

    extraSettingsPattern = '"Spectrograph":\s*{.*?"wavelength":\s*([0-9.e+-]+)';
    tokens = regexp(omeXML, extraSettingsPattern, 'tokens');
    
    if ~isempty(tokens)
        try
            centerWavelength = str2double(tokens{1}{1}) * 1e9;
            wavelengthAxis = linspace(centerWavelength-200, centerWavelength+200, 1024)';
            return;
        catch
        end
    end

    warning('Could not determine wavelength from OME-XML, using default 400-800nm range');
    wavelengthAxis = linspace(400, 800, 1024)';
end

function wavelengthAxis = tryAlternativeWavelengthSources(omeXML, tiffInfo, numWavelengthLayers)
    wavelengthAxis = [];
    for i = 1:length(tiffInfo)
        if isfield(tiffInfo(i), 'PageName') && contains(tiffInfo(i).PageName, 'Spectrum', 'IgnoreCase', true)
            if isfield(tiffInfo(i), 'UnknownTags')
                tags = tiffInfo(i).UnknownTags;
                for j = 1:length(tags)
                    if contains(tags(j).Name, 'wave', 'IgnoreCase', true) && ~isempty(tags(j).Value)
                        try
                            if isnumeric(tags(j).Value)
                                wavelengthValue = tags(j).Value;
                            else
                                wavelengthValue = str2double(tags(j).Value);
                            end
                            if ~isnan(wavelengthValue) && wavelengthValue > 0
                                wavelengthAxis = linspace(wavelengthValue-200, wavelengthValue+200, numWavelengthLayers)';
                                return;
                            end
                        catch
                        end
                    end
                end
            end
        end
    end

    if ~isempty(omeXML)
        rangePattern = '([0-9]+)\s*[-]\s*([0-9]+)\s*nm';
        rangeMatch = regexp(omeXML, rangePattern, 'tokens');
        if ~isempty(rangeMatch)
            try
                startWL = str2double(rangeMatch{1}{1});
                endWL = str2double(rangeMatch{1}{2});
                wavelengthAxis = linspace(startWL, endWL, numWavelengthLayers)';
                return;
            catch
            end
        end
    end
end

function s = xml2struct(xmlFileOrString)
    if exist(xmlFileOrString, 'file')
        xmlDoc = xmlread(xmlFileOrString);
    else
        xmlDoc = xmlread(string(xmlFileOrString));
    end
    s = parseChildNodes(xmlDoc);
end

function children = parseChildNodes(node)
    children = [];
    if node.hasChildNodes
        childNodes = node.getChildNodes;
        numChildNodes = childNodes.getLength;
        children = struct('Name', {}, 'Attributes', {}, 'Data', {}, 'Children', {});
        
        for i = 1:numChildNodes
            theChild = childNodes.item(i-1);
            if theChild.getNodeType == theChild.ELEMENT_NODE
                children(i) = makeStructFromNode(theChild);
            end
        end
    end
end

function nodeStruct = makeStructFromNode(theNode)
    nodeStruct = struct(...
        'Name', char(theNode.getNodeName),...
        'Attributes', parseAttributes(theNode),...
        'Data', '',...
        'Children', parseChildNodes(theNode));
    
    if any(strcmp(methods(theNode), 'getData'))
        nodeStruct.Data = char(theNode.getData());
    else
        nodeStruct.Data = '';
    end
end

function attributes = parseAttributes(theNode)
    attributes = [];
    if theNode.hasAttributes
        theAttributes = theNode.getAttributes;
        numAttributes = theAttributes.getLength;
        attributes = struct('Name', {}, 'Value', {});
        for i = 1:numAttributes
            attrib = theAttributes.item(i-1);
            attributes(i).Name = char(attrib.getName);
            attributes(i).Value = char(attrib.getValue);
        end
    end
end
