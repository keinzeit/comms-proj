function [Zt dctZ I] = ImPreProc(I0)
%function imports the image to be used, processes the image by converting
%to a matrix and scaling the data between 0 and 1 for quantization.  The
%output of this function is a 3d array of 8x8 binary matrices

% no inputs
% outputs Zt, a 3d array of 8x8 matrices

I = I0(801:816,401:416,:);              % Crop Image
I = imresize(I, 0.5);                  % Resize image

%figure(1)
%imshow(I0)                              % Display original image

Z = rgb2gray(I);                        % Convert to grayscale
Z = im2double(Z);                       % Convert image to type double

m = size(Z,1);                          % Get number of rows
n = size(Z,2);                          % number of columns

%figure(2)
%imshow(Z)                               %display the cropped/scaled image

fun = @dct2;
dctZ= blkproc(Z,[8,8],fun);             %Perform DCT on image

%figure(3)
%imshow(dctZ)                            %Display the DCT of image

            % scale the DCT image matrix from 0-1
minZ = min(dctZ(:)); maxZ = max(dctZ(:));
scaledZ = (dctZ-min(dctZ(:)))./(max(dctZ(:)-min(dctZ(:))));


Zt = reshape(scaledZ, [8,8,(m*n)/64]);  % Outputs a 3D array of scaled data
