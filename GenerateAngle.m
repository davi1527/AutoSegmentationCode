function [alignmentangle] = GenerateAngle(ThreeDImage)
k = 1;
while 1
    alignment_image = ThreeDImage(:,:,k);
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
imshow(ThreeDImage(:,:,k))
f = msgbox('Select 2 points: Second point being closest to Acetabular head');
[x, y] = getpts;
alignmentangle = atan2d(y(1)-y(2),x(1)-x(2));
end