function [Zt,dctZ,I,minZ,maxZ] = ImPreProc(I0)
% function [Zt,dctZ,I] = ImPreProc(I0)
% Processes an image by converting
% to a matrix and scaling the data between 0 and 1 for quantization.  The
% output of this function is a 3d matrix of dimensions 8x8xM, where M is
% the number of 8x8 pixel blocks that the image is made from.

% input: an image read in by inmread
% outputs Zt, a 3d matrix of size 8x8xM

% testing value was 32
N = 10;
points = 1:2^N;
I = I0(points+600,points+500,:);              % Crop Image
I = imresize(I, 0.5);                   % Resize image

% figure(1)
% imshow(I0)                              % Display original image
% title('Original Image')

Z = rgb2gray(I);                        % Convert to grayscale
Z = im2double(Z);                       % Convert image to type double

figure(2)
imshow(Z)                               %display the grayscaled cropped/scaled image
title('Cropped Image')

m = size(Z,1);                          % Get number of rows
n = size(Z,2);                          % number of columns

blockdim = 8; % size of block to perform DCT on
fun = @dct2;
dctZ= blkproc(Z,[blockdim,blockdim],fun);             %Perform DCT on image

figure(3)
imshow(dctZ)                            %Display the DCT of image
title('DCT of Cropped Image')

% scale the DCT image matrix from 0-1
minZ = min(dctZ(:)); % min value in entire dctZ matrix
maxZ = max(dctZ(:)); % max value in entire dctZ matrix
scaledZ = (dctZ-minZ)./(maxZ-minZ);

% reshapes scaledZ into a 3D matrix. the first two dimensions are the size
% of the dct block, the 3rd dimension is the number of dct blocks in the
% image. 
% for example, if image is size 12x8 and we use a block size of 4x4
% then the image fits (12/4)*(8/4) = 3*2 = 6 blocks
% thus the reshaped matrix will be size 4x4x6
Zt = reshape(scaledZ, [blockdim,blockdim,(m*n)/blockdim^2]);
