function activityPackage = run_3d_modeling(matFile, gridPoints)
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
oldDir = pwd;
cd(rootDir);
cleanupObj = onCleanup(@() cd(oldDir));
if nargin < 2 || isempty(gridPoints)
    gridPoints = 100;
end
if nargin < 1 || isempty(matFile)
    [file, path] = uigetfile({'*.mat', 'MAT files (*.mat)'}, 'Select Activity-Structure MAT file', ...
        fullfile(rootDir, 'results', 'three_dim_model_package.mat'));
    if isequal(file, 0) || isequal(path, 0)
        activityPackage = [];
        return;
    end
    matFile = fullfile(path, file);
end
activityPackage = precies.workflow('loadActivityPackage', matFile);
legacyVolume = precies.workflow('buildLegacyVolume', activityPackage, gridPoints);
activityVolume = precies.workflow('buildActivityStructureVolume', activityPackage, gridPoints);
precies.workflow('showLegacySliceViewer', legacyVolume);
precies.workflow('showActivityStructureVolume', activityVolume);
fprintf('3D modeling loaded from %s\n', matFile);
end
