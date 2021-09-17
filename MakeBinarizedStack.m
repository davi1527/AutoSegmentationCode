datafolder = uigetdir();
filePattern = fullfile(datafolder, '*.bmp');
files = dir(filePattern);
filename = files(round(length(files)));
for i = 1:length(files)
    image = imread([files(i).folder, '\', files(i).name]);
    binaryimage = imbinarize(image, 0.3);
    if ~isfile(['\\nas01.itap.purdue.edu\puhome\desktop\BNLOutputs\SampleTrab', '_BW.tif'])
            imwrite(binaryimage, ['\\nas01.itap.purdue.edu\puhome\desktop\BNLOutputs\SampleTrab', '_BW.tif'], 'tif', 'Compression', 'none')
            binaryimagestack = binaryimage;
            fprintf('Done.\n')
        else
            imwrite(binaryimage, ['\\nas01.itap.purdue.edu\puhome\desktop\BNLOutputs\SampleTrab', '_BW.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append')
            fprintf('Done.\n')
            binaryimagestack = cat(3, binaryimagestack, binaryimage);
    end
end
