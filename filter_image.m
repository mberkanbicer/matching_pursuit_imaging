org_image = img;
tmp_img = org_image;

[row, col] = size(tmp_img);

degree = 2;
c = degree + 1;

hm = 3;
for h = 1 : hm
    row_img = reshape(tmp_img, [row * col 1]);

    for m = c : (row - c + 1)
        for n = c : (col - c + 1)
            active_pix = [m n];
            active_pix_row = (m - 1) * row + n;

            neighbors = find_neighbor(active_pix, degree);

            [krow, kcol] = size(neighbors);

            for p = 1 : krow
                for r = 1 : kcol
                    if p ~= c && r ~=c
                        neighbor = neighbors{p, r};
                        neighbor_row = (neighbor(1, 1) - 1) * col + neighbor(1, 2);

                        row_img(active_pix_row) = row_img(active_pix_row) + row_img(neighbor_row) / (((c + 1) * (c + 1)) + h - 1);
                    end
                end
            end        
        end
    end
    
    tmp_img = reshape(row_img, [row col]);
end

figure
[matimg, p] = matplot2(1 : cols, 1 : rows, flipud(tmp_img), 40);
xlabel('x-axis (cm)')
ylabel('y-axis (cm)')
title('Filtered Image')
axis equal
axis tight

colormap jet
colorbar
set(gca, 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman')
set(gca,'XTick', linspace(1, imag_dims, axisvc))
set(gca,'XTickLabel', axes_scale)
set(gca,'YTick', linspace(1, imag_dims, axisvc))
set(gca,'YTickLabel', -axes_scale)

drawnow