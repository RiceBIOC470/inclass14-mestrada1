%Inclass 14

%Work with the image stemcells_dapi.tif in this folder

image = 'stemcells_dapi1.tif';
img = imread(image);
imshow(img, []);

% (1) Make a binary mask by thresholding as best you can

img_bw = img > 300;
imshow(img_bw, []);

% (2) Try to separate touching objects using watershed. Use two different
% ways to define the basins. 
% (A) With erosion of the mask 

% Create binary image of "fused" nuclei only
conn_comp = bwconncomp(img_bw);
data = regionprops(conn_comp, 'Area');
area = [data.Area];
fused_cand = area > mean(area) + std(area);
sublist = conn_comp.PixelIdxList(fused_cand);
sublist = cat(1, sublist{:});
mask = false(size(img_bw));
mask(sublist) = 1;
imshow(mask);

% Perform erosion 
s = round(1.2*sqrt(mean(area))/pi);
nucmin = imerode(mask, strel('disk', s));
imshow(nucmin);

outside = ~imdilate(mask, strel('disk', 1));
imshow(outside);

basin = imcomplement(bwdist(outside));
basin = imimposemin(basin, nucmin | outside);
pcolor(basin); shading flat;

ws = watershed(basin);
imshow(ws); colormap('jet'); caxis([0 20]);

new_mask = ws > 1 | (img_bw - mask);
imshow(new_mask);

% (B) with a distance transform. Which works better in this case?

D = bwdist(~img_bw);
imshow(D, []);
D = -D;
D(~img_bw) = -Inf;
imshow(D, []);
ws = watershed(D);
rgb = label2rgb(ws, 'jet', [.5 .5 .5]);
imshow(rgb);

% The distance transform method tends to lead to oversegmentation. Defining
% the basins with erosion of the mask seems to work better in this case.