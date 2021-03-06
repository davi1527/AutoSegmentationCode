%% Intersection Point finder. Finds the trabecular points that are furthest distance from the centroid using line intersection.

function [masked_Image, dilationmask, spline_points] = IntersectionPointFinder(Image, iteration)
% small_outlines = bwpropfilt(Image, 'Perimeter', 1);
% small_point_removing_mask = imcomplement(small_outlines);
% Image = Image.*small_point_removing_mask;
ImageOutline = bwboundaries(Image, 'noholes');
if size(ImageOutline,1) > 1
    Image = bwpropfilt(Image, 'Perimeter', 1);
    ImageOutline = bwboundaries(Image, 'noholes');
end
if size(ImageOutline, 1) == 0
    fprintf('Image %d does not have trabeculae', iteration)
    masked_Image = zeros(size(Image));
    dilationmask = zeros(size(Image));
    spline_points = [];
else
    Shape = ImageOutline{1};
    shapepoints = polyshape(Shape(:,2), Shape(:,1));
    [center_x, center_y] = centroid(shapepoints);
    ImageHoles = bwboundaries(Image, 'holes');
    counter = 1;
    while 1
        if size(ImageHoles,1) < counter
            break
        else
            ll = polyshape(ImageHoles{counter});
            sizing_points = ImageHoles{counter};
            if size(sizing_points,1) == size(ImageOutline{1},1) 
               ImageHoles(counter,:) = [];
            elseif size(ImageHoles{counter},1) <= 1
               ImageHoles(counter,:) = [];
            else 
            end
        end
        counter = counter + 1;
    end
    hold on
    for j = 1:length(ImageOutline)
    boundary_ = ImageOutline{j};
    plot(boundary_(:,2), boundary_(:,1), 'r', 'LineWidth', 2);
    end
    hold on
    for n = 1:length(ImageHoles)
    trabeculae = ImageHoles{n};
    plot(trabeculae(:,2), trabeculae(:,1), 'b', 'LineWidth', 2);
    end
    Length = 225;
    % for i = 0:2:360
    % x_point = center_x+(Length*cosd(i));
    % y_point = center_y+(Length*sind(i));
    % Line1 = line([center_point(1)/2, x_point], [center_point(2)/2, y_point], 'Color', 'g', 'LineWidth',2);
    % % end
    % plot(h)
    % xdataholes = f(:,2);
    % ydataholes = f(:,1);
    % h = polyshape(xdataholes,ydataholes);
    for j = 1:length(ImageHoles)
    f = ImageHoles{j};
    xdataholes = f(:,2);
    ydataholes = f(:,1);
    if exist('total_x', 'var') == 0
        total_x = xdataholes;
        total_y = ydataholes;
    else
        total_x = cat(1, total_x, xdataholes);
        total_y = cat(1, total_y, ydataholes);
    end
        for i = 0:2:360
            x_point = center_x+(Length*cosd(i));
            y_point = center_y+(Length*sind(i));
            xdataline = [center_x, x_point];
            ydataline = [center_y, y_point];
    %         Line1 = line([center_x, x_point], [center_y, y_point], 'Color', 'g', 'LineWidth',2);
    %         xdataline = get(Line1, 'XData');
    %         ydataline = get(Line1, 'YData');
            [imagepoints_x, imagepoints_y] = polyxpoly(xdataline,ydataline, xdataholes, ydataholes, 'unique');
            if isempty(imagepoints_x)
            else
                distance = zeros(length(imagepoints_x),1);
                for g = 1:length(imagepoints_x)
                    distance(g) = sqrt((center_x-imagepoints_x(g)).^2+(center_y-imagepoints_y(g)).^2);
                    if length(imagepoints_x) > 1
                        distance_between = zeros(length(imagepoints_x)-1, 1);
                        for k = 1:length(imagepoints_x)-1
                            distance_between(k) = sqrt((imagepoints_x(k) - imagepoints_x(k+1)).^2 + (imagepoints_y(k)-imagepoints_y(k+1)).^2); 
                        end
                    else
                    end
                end
                
                if exist('distance_between', 'var') == 1
                   distance_between_thresh = find(distance_between > 160);
                   if isempty(distance_between_thresh)
                   else
                   imagepoints_x(distance_between_thresh) = [];
                   imagepoints_y(distance_between_thresh) = [];
                   distance(distance_between_thresh) = [];
                   end
                else  
                end
                MaxValue = max(distance);
                index = find(distance == MaxValue);
                if exist('splinepoints', 'var') == 0
                    splinepoints = [imagepoints_x(index(1)), imagepoints_y(index(1)), iteration, distance(index(1))];
                else
                    splinepoints = cat(1,splinepoints, [imagepoints_x(index(1)), imagepoints_y(index(1)), iteration, distance(index(1))]);

                end
            end
        end
    end
    if exist('splinepoints', 'var') == 0
        fprintf('Slice %d does not have trabeculae\nDefin', iteration);
        masked_Image = zeros(size(Image));
        dilationmask = zeros(size(Image));
        spline_points = [];
    else
        distance_threshold = max(splinepoints(:,4))*1/2;
        refinedpoints = find(splinepoints(:,4)>distance_threshold);
        newsplinepoints = splinepoints(refinedpoints,:);
        points1 = transpose(newsplinepoints);
        splineinput = [points1(1,:); points1(2,:)];
        tangent_y = splineinput(2,:) - center_y;
        tangent_x = splineinput(1,:) - center_x;
        angle_ = atan2d(tangent_y, tangent_x);
        angle_h = cat(1, splineinput, angle_);
        BB = transpose(angle_h); 
        BBB = sortrows(BB, 3);
        sorted_matrix = transpose(BBB);
        splineinput_refined = [sorted_matrix(1,:); sorted_matrix(2,:)];
        spline_points = [transpose(sorted_matrix(1,:)) transpose(sorted_matrix(2,:))];
        spline_points(:,3) = iteration;
        cscvn(splineinput_refined)
        RawPolygon = polyshape(sorted_matrix(1,:), sorted_matrix(2,:));
%         RefinedPolygon = polybuffer(RawPolygon, 2, 'JointType', 'round');
        polygon_x = RawPolygon.Vertices(:,1);
        polygon_y = RawPolygon.Vertices(:,2);
%         bp = boundary(total_x, total_y);
      
        try 
            smoothed_y = smooth(polygon_y, 5);
            smoothed_x = smooth(polygon_x, 5);
        catch 
            input('g');
        end
        if ~isempty(smoothed_x) && ~isempty(smoothed_y)
            [~, ~, k]= curvature([smoothed_x smoothed_y]);
            tangent_y2 = smoothed_y - center_y;
            tangent_x2 = smoothed_x - center_x;
            angle_matrix = atan2(tangent_y2, tangent_x2);
            for i = 1:size(angle_matrix,1)
                while abs(k(i,1)) > 0.25
                    smoothed_x(i) = smoothed_x(i) + (2*cos(angle_matrix(i)));
                    smoothed_y(i) = smoothed_y(i) + (2*sin(angle_matrix(i)));
                    [~, ~, k]= curvature([smoothed_x smoothed_y]);
                end
                while abs(k(i,2)) > 0.25
                    smoothed_x(i) = smoothed_x(i) + (2*cos(angle_matrix(i)));
                    smoothed_y(i) = smoothed_y(i) + (2*sin(angle_matrix(i)));
                    [~, ~, k]= curvature([smoothed_x smoothed_y]);
                end
            end
        else
        end
        gg = isnan(smoothed_x);
        if sum(gg) > 0
            index = find(gg == 1);
            smoothed_x(index) = [];
            smoothed_y(index) = [];
        end
        spline_mask = poly2mask(smoothed_x, smoothed_y, size(Image,1), size(Image,2));
%         refinedsplinemask = poly2mask(newsplinepoints(:,1), newsplinepoints(:,2), size(Image,1), size(Image,2));
%         boundarymask  = poly2mask(total_x(bp), total_y(bp), size(Image,1), size(Image,2));
%         spline_mask(boundarymask ~= 0) = true;
%         spline_mask(refinedsplinemask ~= 0) = true;
        spline_mask = imfill(spline_mask, 'holes');
        errosion = imerode(spline_mask, strel('disk', 5));
        dilationmask = imdilate(errosion, strel('disk', 10));
        masked_Image = dilationmask .* Image;
    end
end
end