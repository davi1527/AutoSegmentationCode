%% 3D Trabecular Image Stacking
% Written by: Zachary Davis
% Chan Laboratory at Purdue University
% Version 1.0
% Last edit: 5/11/2021

function GenerateBWStack(filePattern)

%% Reading in Files and drawing ROI using imrect function 
mkdir TifImages
files = dir(filePattern);

%% Converting from .bmp files to .tif files. *ONLY USE IF RAW RECONSTRUCTIONS ARE .bmp FILES*
% for j = 1 : length(files)
%     fprintf('Writing %d of %d...', j, length(files))
%     convertimage = imread([files(j).folder,'\',files(j).name]);
%     filename_ = files(j).name;
%     trimmedname_=filename_(1:strfind(filename,'.bmp')-1);
%     imwrite(convertimage, [files(j).folder,'\TifImages\', trimmedname_ ,'_tif.tif'], 'tif', 'Compression', 'none')
%     fprintf('Done.\n')
% end
% datafolder = uigetdir();
% filedirectory = fullfile(datafolder, '*.tif');
% files = dir(filedirectory);

%% Finding what slice shows the ROI to where you can manually identify it
% Goes through the pictures and opens based on thresholdvalue. Let the program know whether you can see 
% the ROI when prompted in command window typing Y or N.
% sizing_image = imread([files(1).folder, '\', files(1).name]);
% image_stack = zeros(size(sizing_image,1), size(sizing_image,2), length(files));
% for i = 1:length(files)
%     filename_ = files(i).name;
%     image_ = imread([files(i).folder,'\',filename_]);
%     thresholded_value = graythresh(image_);
%     big_binary_image = imbinarize(image_, thresholded_value);
%     image_stack(:,:,i) = big_binary_image;
%     disp(i);
% end
uiwait(msgbox('Draw rectangle over region of interest. Look at Command Window'));
k = 1;
while 1
    rectangle_image = imread([files(1).folder, '\', files(k).name]);
    threshold = graythresh(rectangle_image);
    if threshold > 0.25
        imshow(rectangle_image)
        x = input('Can you see the region of interest (Y/N)?', 's');
        lower_x = lower(x);
        if lower_x == 'y'
            break
        elseif lower_x == 'n'
            k = k+25;
        else
            disp('Please enter Y or N')
        end
    else
        k = k+25;
        disp(k)
        disp(threshold)
    end
    close 
end

%% Manually drawing rectangle completely covering ROI.
%Give a little extra space on the sides to guarantee full ROI is being
%included
imshow(rectangle_image)
rectangular_ROI = drawrectangle;
xmin = round(rectangular_ROI.Vertices(1,1));
xmax = round(rectangular_ROI.Vertices(4,1));
ymin = round(rectangular_ROI.Vertices(1,2));
ymax = round(rectangular_ROI.Vertices(2,2));
close
%% Reading in each file in the folder. Only the ROI drawn before using imrect will be read to save space. 
%The order is: read file --> smooth using gaussian filter --> exclude blank images based on threshold value --> smooth images --> binarize --> stack into tif stack 
k = 1;
uiwait(msgbox('Draw two points to align femoral neck. 2nd point closest to head. Look at Command Window'));
while 1
    alignment_image = imread([files(1).folder, '\', files(k).name]);
    imshow(alignment_image)
    x = input('Can you see the region of interest (Y/N)?', 's');
    lower_x = lower(x);
        if lower_x == 'y'
            break
        elseif lower_x == 'n'
                k = k+25;
        else
            disp('Please enter Y or N')
        end
        close 
end
[x, y] = getpts;
alpha_ = atan2d(y(1)-y(2),x(1)-x(2));
close
for i = 1 : length(files)
    fprintf('Processing %d of %d...', i, length(files))
    filename = files(i).name;
    rawimage_ = imread([files(i).folder,'\',filename],'PixelRegion',{[ymin-25,ymax+25],[xmin-25,xmax+25]});
%     smoothedimage = imgaussfilt(rawimage_, 2);
    initialrotation = imrotate(rawimage_, alpha_, 'bilinear', 'crop');
    threshold_weight = graythresh(rawimage_);
    if threshold_weight < 0.12
        disp('excluded')
    else
        binaryimage = imbinarize(rawimage_, threshold_weight-0.05);
        rotatedimage = imrotate(binaryimage, alpha_, 'bilinear', 'crop');
        windowSize=2;  
        kernel=ones(windowSize)/windowSize^2;
        result=conv2(single(rotatedimage),kernel,'same');
        result=result>0.5;
        rotatedimage(~result)=0; 
%         flippedbinaryimage = ~binaryimage;
        %figure, imshow(binaryimage);
        %imshow(smoothedbinary)
        trimmedname=filename(1:3);
%       imwrite(binaryimage, [files(i).folder,'\BWImages\', trimmedname, '_BW.bmp'])
        if ~isfile([files(i).folder, trimmedname, 'TifStack', '_BW.tif'])
%             BW3D = rotatedimage;
            rawimagestack = initialrotation;
            imwrite(rotatedimage, [files(i).folder, trimmedname 'TifStack', '_BW.tif'], 'tif', 'Compression', 'none')
%             imwrite(flippedbinaryimage, [files(i).folder, trimmedname, 'TifStackFlipped', '_BW.tif'], 'tif', 'Compression', 'none')  
            fprintf('Done.\n')
        else
%             imwrite(rotatedimage, [files(i).folder, trimmedname, 'TifStack', '_BW.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append')
%             imwrite(flippedbinaryimage, [files(i).folder, trimmedname, 'TifStackFlipped', '_BW.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append')
%             BW3D = cat(3, BW3D, rotatedimage);
            rawimagestack = cat(3, rawimagestack, initialrotation);
            fprintf('Done.\n')
        end
        
%     imwrite(rawimage_, [files(i).folder, 'CroppedRawImage', trimmedname, '.tif'], 'tif', 'Compression', 'none')
    end
   
end
uiwait(msgbox('Draw two points to align femoral neck. 2nd point closest to head. Look at Command Window'));
k = 1;
while 1
    YZ_AlignmentImage = permute(rawimagestack(:,k,:), [1,3,2]);
    imshow(YZ_AlignmentImage)
    x = input('Can you see the region of interest (Y/N)?', 's');
    lower_x = lower(x);
        if lower_x == 'y'
            break
        elseif lower_x == 'n'
                k = k+25;
        else
            disp('Please enter Y or N')
        end
        close 
end
[y_yz, z_yz] = getpts;
gamma_ = atan2d(z_yz(1)-z_yz(2),y_yz(1)-y_yz(2));
while 1
    XZ_AlignmentImage = permute(rawimagestack(k,:,:), [2,3,1]);
    imshow(XZ_AlignmentImage)
    x = input('Can you see the region of interest (Y/N)?', 's');
    lower_x = lower(x);
        if lower_x == 'y'
            break
        elseif lower_x == 'n'
                k = k+25;
        else
            disp('Please enter Y or N')
        end
        close 
end
[x_xz, z_xz] = getpts;
beta_ = atan2d(z_xz(1)-z_xz(2),x_xz(1)-x_xz(2));
gg  = imrotate3(rawimagestack, beta_, [0 1 0], 'cubic', 'loose', 'FillValues', 0);
AlignedImageStack = imrotate3(gg, gamma_, [1 0 0], 'cubic', 'loose' , 'FillValues', 0);


%% Reducing Field of view to only include the Region of interest
k = 1;
uiwait(msgbox('Draw rectangle over region of interest. Look at Command Window'));
while 1
    alignment_image = AlignedImageStack(:,:,k);
    imshow(alignment_image)
    x = input('Can you see the region of interest (Y/N)?', 's');
    lower_x = lower(x);
        if lower_x == 'y'
            break
        elseif lower_x == 'n'
                k = k+25;
        else
            disp('Please enter Y or N')
        end
        close 
end
rectangle_trim_points = drawrectangle;
xmin = round(rectangle_trim_points.Vertices(1,1));
xmax = round(rectangle_trim_points.Vertices(4,1));
AlignedImageStack(:,1:xmin,:) = [];
AlignedImageStack(:,xmax-xmin:size(AlignedImageStack,2),:) = [];

% TrimmedAlignedStack = zeros(size(AlignedImageStack,1), xmax-xmin+1, size(AlignedImageStack,3));
% for i = 1:size(AlignedImageStack,3)
%     TrimmedAlignedStack(:,:,i) = AlignedImageStack(:, xmin:xmax, i);
% end

%% Radial acquisition of trabecular points that are maximum distance away from centroid of cortical shell.
% Gaussian Filter --> threshold --> Binarize --> acquire radial points
ThreeDmask = zeros(size(AlignedImageStack));
MaskedImageStack = zeros(size(AlignedImageStack));
% GaussFiltStack = smooth3(AlignedImageStack, 'gaussian', 5);
stl = strel('sphere', 2);
stl2 = strel('sphere', 5);
ThreeDThreshold = graythresh(AlignedImageStack);
Base_Aligned_BW3D = imbinarize(AlignedImageStack, ThreeDThreshold);
ErrodedStack = imerode(Base_Aligned_BW3D, stl);
DilatedStack = imdilate(ErrodedStack, stl2);
% errodedThreshold = graythresh(ErrodedStack);
% Binarized_Erroded_Stack = imbinarize(ErrodedStack, errodedThreshold);
% Erroded_Threshold = graythresh(ErrodedStack);
% Processing_Stack = imbinarize(ErrodedStack, Erroded_Threshold);
for i = 1:size(AlignedImageStack,3)
    disp(i)
    BinaryImage = DilatedStack(:,:,i);
    threshold_weight = graythresh(AlignedImageStack(:,:,i));
%     BinaryImage = Processing_Stack(:,:,i);
    if threshold_weight < 0.12
        ThreeDmask(:,:,i) = 0;
        MaskedImageStack(:,:,i) = 0;
%         if exist('Aligned_BW3D', 'var') == 0
%             Aligned_BW3D = BinaryImage;
%         else 
%             Aligned_BW3D = cat(3, Aligned_BW3D, BinaryImage);
%         end
    else 
%         windowSize=2;             Binary Smoothing: Use only if needed
%         kernel=ones(windowSize)/windowSize^2;
%         result=conv2(single(BinaryImage),kernel,'same');
%         result = result > 0.5;
%         BinaryImage(~result) = 0;
        [masked_image,splinemask, spline_points] = IntersectionPointFinder(BinaryImage, i);
        ThreeDmask(:,:,i) = splinemask;
        MaskedImageStack(:,:,i) = masked_image;
        if exist('ThreeDSpline', 'var') == 0
            ThreeDSpline = spline_points;
        else
            ThreeDSpline = cat(1, ThreeDSpline, spline_points);
        end
%         if exist('Aligned_BW3D', 'var') == 0
%             Aligned_BW3D = BinaryImage;
%         else 
%             Aligned_BW3D = cat(3, Aligned_BW3D, BinaryImage);
%         end
    end
    
end
% thresh_ = graythresh(GaussFiltMask);
% binaryGaussFiltMask = imbinarize(GaussFiltMask, thresh_); 
% FilteredMaskedImageStack = binaryGaussFiltMask .* Aligned_BW3D;
%% Writing the output 3D stack into a 3D tiff stack to input into ImageJ
N = 7;
kernel = ones(N, N, N) / N^3;
blurryImage = convn(ThreeDmask, kernel, 'same');
newBinaryImage = blurryImage > 0.5;
Erroded_Mask = imdilate(newBinaryImage, strel('sphere', 5));
hhh1 = Erroded_Mask .* Base_Aligned_BW3D;
for i = 1: size(ThreeDmask,3)
    fprintf('Writing %d of %d...', i, size(ThreeDmask,3))
%     thresholded_mask = imbinarize(newBinaryImage(:,:,i), 0.5);
    Stack1 = Base_Aligned_BW3D(:,:,i) .* Erroded_Mask(:,:,i);
    fused_image = imfuse(Base_Aligned_BW3D(:,:,i), Stack1);
    binaryImage1 = imbinarize(hhh1(:,:,i), 0.6);
     if ~isfile([files(i).folder, trimmedname, 'Box7BW3D', '_BW.tif'])
            imwrite(binaryImage1, [files(i).folder, trimmedname, 'Box7BW3D', '_BW.tif'], 'tif', 'Compression', 'none')
%             imwrite(Stack1, [files(i).folder, trimmedname 'FilteredStackBox9BufferedPost', '_BW.tif'], 'tif', 'Compression', 'none')
            imwrite(fused_image, [files(i).folder, trimmedname 'Box7BW3DFused', '_BW.tif'], 'tif', 'Compression', 'none')
            fprintf('Done.\n')
        else
            imwrite(binaryImage1, [files(i).folder, trimmedname, 'Box7BW3D', '_BW.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append')
%             imwrite(Stack1, [files(i).folder, trimmedname 'FilteredStackBox9BufferedPost', '_BW.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append')
            imwrite(fused_image, [files(i).folder, trimmedname 'Box7BW3DFused', '_BW.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append')
            fprintf('Done.\n')
     end
end
%% Elipsoid fit code (not necessary)
% [center, radii, evecs, v, chi2] = ellipsoid_fit(elipsepoints,'0');
% xcenter = center(1);
% ycenter = center(2);
% zcenter = center(3);
% xrad = radii(1);
% yrad = radii(2);
% zrad = radii(3);
% mask = zeros(size(BW3D));
%     for x=1:size(BW3D,2)
%         for y=1:size(BW3D,1)
%             for z=1:size(BW3D,3)
%                 if ( ((x-xcenter)/xrad)*((x-xcenter)/xrad) + ((y-ycenter)/yrad)*((y-ycenter)/yrad) + ((z-zcenter)/zrad)*((z-zcenter)/zrad) < 1 )
%                     mask(x,y,z) = 1; %//set elements within ellipsoid to 1
%                 end
%             end
%         end
%     end
% maskedimage = BW3D.*mask;

end 