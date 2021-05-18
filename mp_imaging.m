%% Matching Pursuit Algorithm for Microwave Breast Imaging

close all
clear all
clc

%% measure time
tic

%% initialize geometry defaults
disp('Process: Geometry Definations..')
as_dist         = 2;
radii_skin      = 7;
radii_tumor     = 1;
tumxy_1         = [2.5593 4.8920]; % randomly chosen (x, y)
tumxy_2         = [0.1609 -5.5168]; % randomly chosen (x, y)
tumxy_3         = [-0.6747 0.1048]; % randomly chosen (x, y)
antenna_count   = 360;
skin_point_num  = 64; % number of skin scatterers
tumor_point_num = 8; % number of tumor scatterers

radii_anlayer   = radii_skin + as_dist; % radius of the antenna layer
center_skin     = [radii_anlayer radii_anlayer];
tumorc_relxy_1  = [center_skin(1,1) + tumxy_1(1,1) center_skin(1,2) + tumxy_1(1,2)]; % relative center of tumor
tumorc_relxy_2  = [center_skin(1,1) + tumxy_2(1,1) center_skin(1,2) + tumxy_2(1,2)]; % relative center of tumor
tumorc_relxy_3  = [center_skin(1,1) + tumxy_3(1,1) center_skin(1,2) + tumxy_3(1,2)]; % relative center of tumor

pos_anlayer     = pcircle(center_skin, radii_anlayer, antenna_count); % all positions of the antenna
pos_skin        = pcircle(center_skin, radii_skin, skin_point_num); % positions of the skin scatterers

disp('Geometry Information:')
disp('Tumor #1')
disp(['R: ' num2str(radii_anlayer + abs(tumxy_1(1) + 1j*tumxy_1(2))) ', Angle: ' num2str(angle(tumxy_1(1) + 1j*tumxy_1(2)) * 180 / pi)])
disp('')
disp('Tumor #2')
disp(['R: ' num2str(radii_anlayer + abs(tumxy_2(1) + 1j*tumxy_2(2))) ', Angle: ' num2str(angle(tumxy_2(1) + 1j*tumxy_2(2)) * 180 / pi)])
disp('')
disp('Tumor #3')
disp(['R: ' num2str(radii_anlayer + abs(tumxy_3(1) + 1j*tumxy_3(2))) ', Angle: ' num2str(angle(tumxy_3(1) + 1j*tumxy_3(2)) * 180 / pi)])

%% initialize radar defaults
disp('Process: Radar Definations..')
c    = 299792458; % velocity of the light
fmin = 3.1 * 10^9;
fmax = 10.6 * 10^9;
f    = linspace(fmin, fmax, (fmax - fmin) / (5 * 10^7) + 1); % frequency array

epsR = 1; % relative permittivity of the air (for the sake of simplicity, medium is considered to be homogeneous, filled with air)
v    = c / sqrt(epsR); % velocity of the wave in medium

k    = (2 * pi * f) / v; % wavenumber

df   = f(2) - f(1); % frequency resolution
N    = length(f);   % number of frequencies
dr   = v / (2*N*df);% range resolution
R    = 0 : dr : dr * (N - 1); % axis scaling
%% generate scattering electric field data
[~, any_count]  = size(pos_anlayer);
[~, sy_count]   = size(pos_skin);
[~, f_count]    = size(f);
es              = zeros(f_count, any_count);
tumor_count     = 3;

disp('Process: Calculating Es of Tumor ..')
for m = 1 : any_count
    anxy = [pos_anlayer(1, m) pos_anlayer(2, m)];
    
    if tumor_count == 1
        txy_1  = [tumorc_relxy_1(1, 1) tumorc_relxy_1(1, 2)];
        
        es(:, m) = 2 * (exp(-1i * 0.02 * k * norm(anxy - txy_1)).');
    elseif tumor_count == 2
        txy_1  = [tumorc_relxy_1(1, 1) tumorc_relxy_1(1, 2)];
        txy_2 = [tumorc_relxy_2(1, 1) tumorc_relxy_2(1, 2)];
        
        es(:, m) = 2 * (exp(-1i * 0.02 * k * norm(anxy - txy_1)).')...
                 + 3 * (exp(-1i * 0.02 * k * norm(anxy - txy_2)).');
    elseif tumor_count == 3
        txy_1  = [tumorc_relxy_1(1, 1) tumorc_relxy_1(1, 2)];
        txy_2 = [tumorc_relxy_2(1, 1) tumorc_relxy_2(1, 2)];
        txy_3 = [tumorc_relxy_3(1, 1) tumorc_relxy_3(1, 2)];
        
        es(:, m) = 2 * (exp(-1i * 0.02 * k * norm(anxy - txy_1)).')...
                 + 3 * (exp(-1i * 0.02 * k * norm(anxy - txy_2)).')...
                 + 1 * (exp(-1i * 0.02 * k * norm(anxy - txy_3)).');
    end
end
es = awgn(es, 10); % add some noise
%%
disp('Process: Meshing and Matching..')
%% meshing process
disp('Process: Mesh Generation..')
mesh_distance = 0.05;
fd = @(p) sqrt(sum(p .^ 2, 2)) - 1;
[p, t] = distmesh2d(fd, @huniform, mesh_distance, [-1, -1;1, 1], []);
title('Generated Mesh')

disp(['Number of points: ' num2str(size(p))])
disp('Process: Mesh Generated..')

points = p * (radii_anlayer - 0.1);
points(:, 1) = points(:, 1) + center_skin(1, 1);
points(:, 2) = points(:, 2) + center_skin(1, 2);

[anrows, ancols] = size(pos_anlayer);
[prows, pcols] = size(points);

if anrows == 2 && ancols > 2
    pos_anlayer = pos_anlayer';
    
    [anrows, ancols] = size(pos_anlayer);
end

%% Matching Pursuit based Imaging
disp('Process: Matching..')

[esrow, escol]   = size(es);

fi_mp = [];

disp('Process: Calculating Matching-Es..')
for m = 1:prows
    ptxy = [points(m, 1) points(m, 2)];

    for n = 1 : anrows
        anxy = [pos_anlayer(n, 1) pos_anlayer(n, 2)];
        fi_mp(:, n, m) = exp(-1i * 0.02 * k * norm(ptxy - anxy));
    end

    L = sum(sum(es .* conj(fi_mp(:, :, m))));

    Limg(m, :) = [ptxy L];
end
   

%% plot images
disp('Process: Imaging..')
axisvc = 5;
axes_scale  = linspace(-(radii_anlayer), (radii_anlayer), axisvc);

figure
imag_dims = 161;
mult = (imag_dims - 1) / 16;
img = zeros(imag_dims);

img_pixels = round(Limg(:, 1:2) * mult);
img_vals   = Limg(:, 3);

for m = 1:prows
    img(img_pixels(m, 1) + 1, img_pixels(m, 2) + 1) = abs(img_vals(m, 1));
end

img = rot90(flipud(img), 3);

[rows, cols] = size(img);
[matimg, p] = matplot2(1 : cols, 1 : rows, flipud(img), 40);
xlabel('x-axis (cm)')
ylabel('y-axis (cm)')
title('Unfiltered Image')
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

%% filter image
filter_image

%% end of the time measurement
toc