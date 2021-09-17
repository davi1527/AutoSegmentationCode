
 matrix_alpha = zeros(size(binaryimagestack));
 matrix_beta = zeros(size(binaryimagestack));
 matrix_gamma = zeros(size(binaryimagestack));
 angles_matrix = zeros(length(coordinate_trab), 3);
 
for i = 1:length(coordinate_trab)
    disp(i)
    x1 = coordinate_trab(i,1);
    x2 = coordinate_trab(i,4);
    y1 = coordinate_trab(i,2);
    y2 = coordinate_trab(i,5);
    z1 = coordinate_trab(i,3);
    z2 = coordinate_trab(i,6);
    xpoints = linspace(round(x1), round(x2), 50);
    ypoints = linspace(round(y1), round(y2), 50);
    zpoints = linspace(round(z1), round(z2), 50);
    alpha = atan2d((y2-y1), (x2-x1));
    beta = atan2d((x2-x1), (y2-y1));
    gamma = atan2d((z2-z1), (y2-y1));
    true_alpha = alpha+180;
    true_beta = beta+180;
    true_gamma = gamma+180;
    angles_matrix(i,:) = [alpha beta gamma];
    for j = 1:50
        disp(j)
        if zpoints(j) > 3.5 && xpoints(j) > 3.5 && ypoints(j) > 3.5
            matrix_alpha(round(ypoints(j))-3:round(ypoints(j))+3, round(xpoints(j))-3:round(xpoints(j))+3, round(zpoints(j))-3:round(zpoints(j))+3) = true_alpha/360;
            matrix_beta(round(ypoints(j))-3:round(ypoints(j))+3, round(xpoints(j))-3:round(xpoints(j))+3, round(zpoints(j))-3:round(zpoints(j))+3) = true_beta/360;
            matrix_gamma(round(ypoints(j))-3:round(ypoints(j))+3, round(xpoints(j))-3:round(xpoints(j))+3, round(zpoints(j))-3:round(zpoints(j))+3) = true_gamma/360;
        end
    end
    color_matrix = [0.9, true_gamma/360, 0.9];
    hold on
    plot3([x1 x2], [y1 y2], [z1 z2], 'Color', color_matrix, 'LineWidth', 4);
end

for i = 1:size(binaryimagestack,3)
    for j = 1:size(binaryimagestack, 2)
        for k = 1:size(binaryimagestack,1)
            if binaryimagestack(k,j,i) == 1
                hold on;
                plot3(j, k, i, 'Color', [1, 1, 1])
            end
        end
    end
end
            