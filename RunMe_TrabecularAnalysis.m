%% Use to run trabecular analysis code. Choose the folder where the images are stored.
function RunMe_TrabecularAnalysis
datafolder = uigetdir();
filePattern = fullfile(datafolder, '*.tif');

if isempty(filePattern)
    error("The file you have selected has no .tif files in it")
end
GenerateBWStack(filePattern);
pause(5);
% threeDrender = make3D(BWimages);






