function [elipsepoints] = generatepoints(ThreeDmatrix)
k = 1;
while 1
    if k > size(ThreeDmatrix, 3)
        break
    else
        imshow(ThreeDmatrix(:,:,k))
        [xi, yi] = getpts;
        zi = zeros(length(xi), 1);
        zi(:,1) = k;
        addedpoints = cat(2, xi, yi, zi);
        if isempty(addedpoints)
        else
            if exist('elipsepoints', 'var') == 0
            elipsepoints = addedpoints;
            else 
            elipsepoints = cat(1,elipsepoints, addedpoints);
            end
        end
        
    end
    k = k + 15;
end
% imshow(ThreeDmatrix(:,:,1))
% [xi, yi] = getpts;
% zi = zeros(length(xi), 1);
% zi(:,1) = 1;
% elipsepoints = cat(2, xi, yi, zi);
% for i = 2:round(size(ThreeDmatrix,3)/10)-1
%     imshow(ThreeDmatrix(:,:,i*10))
%     [xi, yi] = getpts;
%     zi = zeros(length(xi), 1);
%     zi(:,i) = 1;
%     newpoints_ = cat(2, xi, yi, zi);
%     elipsepoints = cat(1, elipsepoints, newpoints_);
% end
% 
% end
