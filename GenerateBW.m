function GenerateBWStack(filePattern)

mkdir BWImages 
mkdir TifImages
files = dir(filePattern);
filename = files(round(length(files)/2)).name;
rectangle_image = imread([files(1).folder, '\', filename]);
imshow(rectangle_image)
rectangular_ROI = drawrectangle;
xmin = round(rectangular_ROI.Vertices(1,1));
xmax = round(rectangular_ROI.Vertices(4,1));
ymin = round(rectangular_ROI.Vertices(1,2));
ymax = round(rectangular_ROI.Vertices(2,2));
BW3D = zeros(ymax-ymin+1, xmax-xmin+1, length(files));
close
% for j = 1 : length(files)-1
%     convertimage = imread([files(j).folder,'\',files(j).name]);
%     filename_ = files(j).name;
%     trimmedname_=filename_(1:strfind(filename,'.bmp')-1);
%     imwrite(convertimage, [files(j).folder,'\TifImages\', trimmedname_ ,'_tif.tif'], 'tif', 'Compression', 'none')
%     
% end
% datafolder = uigetdir();
% filedirectory = fullfile(datafolder, '*.tif');
% files_ = dir(filedirectory);
for i = 1 : length(files)
    filename = files(i).name;
    rawimage_ = imread([files(i).folder,'\',filename],'PixelRegion',{[ymin,ymax],[xmin,xmax]});
    smoothedimage = imgaussfilt(rawimage_, 2);
    threshold_weight = graythresh(smoothedimage);
    if threshold_weight < 0.15
        binaryimage = imbinarize(smoothedimage, 0.2);
    else
        binaryimage = imbinarize(smoothedimage, threshold_weight-0.05);
    end
    windowSize=2;  
    kernel=ones(windowSize)/windowSize^2;
    result=conv2(single(binaryimage),kernel,'same');
    result=result>0.5;
    binaryimage(~result)=0; 
    flippedbinaryimage = ~binaryimage;
    %figure, imshow(binaryimage);
    %imshow(smoothedbinary)
    trimmedname=filename(1:strfind(filename,'tif')-1);
%     imwrite(binaryimage, [files(i).folder,'\BWImages\', trimmedname, '_BW.bmp'])
    if i == 1
        imwrite(binaryimage, [files(i).folder, trimmedname, 'TifStack', '_BW.tif'], 'tif', 'Compression', 'none')
        imwrite(flippedbinaryimage, [files(i).folder, trimmedname, 'TifStackRealFlipped', '_BW.tif'], 'tif', 'Compression', 'none')
        
    else
        imwrite(binaryimage, [files(i).folder, 'TifStackReal', '_BW.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append')
        imwrite(flippedbinaryimage, [files(i).folder, trimmedname, 'TifStackFlipped', '_BW.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append')
    end
%     imwrite(rawimage_, [files(i).folder, 'CroppedRawImage', trimmedname, '.tif'], 'tif', 'Compression', 'none')
%     BW3D(:,:,i) = binaryimage(:,:);
    disp(i)
end 