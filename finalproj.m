%% EE107 Final Project
% Kenny Yau, Jon Manninen
%
% This project is a simulation of the model of communication systems for
% images, such as that used for digital HDTV. It consists of a number of
% parts:
%     - Image pre-processing
%     - Conversion to bit-stream
%     - Modulation
%     - Channel Effects
%     - Noise
%     - Matched Filter
%     - Receiver Equalization
%     - Sampling and Detection
%     - Converion to Image
%     - Image post-processing

% Stuff to do:

%% Image pre-process:
% Find own image
%{
filename = 'heartstringtug.jpeg';

I = imread(filename);
I = imrotate(I,270,'bilinear');
figure(1); imshow(I); title('Original Image')
Z = rgb2gray(I); % convert to grayscale
figure(2); imshow(Z); title('Grayscaled Image')
%}

% convert image to double data type using im2double
% make sure image is m-by-n, where both m and n are divisible by 8
% perform DCT image in 8x8 blocks using blkproc
% scale all values to fit between 0 and 1
% quantize the scaled, transformed image using code from hw 7
% use reshape and permute to convert the scaled image? into what?
%     - reshape reorganizes a matrix with certain dimensions and size into a matrix of another dimension and size
%     - we are organizing the image data into 8x8 blocks by stacking them on top of each other
%     - use permute to make sure first two dimensions specify a single 8x8 block
% 
%% Conversion to bit-stream:
% we transmit N DCT blocks at a time
% Convert N DCT blocks into a single column vector
% convert each element in this column vector into an 8-bit binary number using de2bi, stored as a row
%     - ex [5] -> [1 0 1 0 0 0 0 0]

% ex start with N=2 blocks of 8x8 DCT blocks

% rows = 8; cols = 8;
% dct1 = randi(20,rows,cols)
% dct2 = randi(20,rows,cols)
% dct(:,:,1) = dct1;
% dct(:,:,2) = dct2;

% Assumptions: dct is a bunch of 8x8 dct blocks stacked on top of each other
% and that the values in dct are quantized to integer values

% get dimensions of dct
% [M,N,B] = size(dct);
% 
% len = M*N;
% bs = zeros(M*N*B,1); %preallocate memory for vector
% for i = 1:B
%     bs_temp = reshape(dct(:,:,i),len,1); % reshape B MxN dct blocks into a long column vector
%     start = len*(i-1)+1
%     stop  = len*i
%     bs(start:stop) = bs_temp; % put values into correct location in bs
% end
% 
% bs = de2bi(bs,8); % convert quantized values into 8-bit numbers. right MSB; use correct flag for left MSB 



%% final project
I = imread('file.jpeg');                % Import Image
I = imrotate(I, 270);                   % Rotate Image
I = imresize(I, 0.25);                  % Resize image

%figure(1)
%imshow(I)                               % Display original image

Z = rgb2gray(I);                        % Convert to grayscale
Z = im2double(Z);                       % Convert image to type double

m = size(Z,1);                          % Get number of rows
n = size(Z,2);                          % number of columns

%fun = @(block_struct)(block_struct.data);
fun = @dct2;
dctZ= blkproc(Z,[8,8],fun);             %Perform DCT on image

figure(2)
imshow(dctZ)                            %Display the DCT of image

% scale the DCT image matrix from 0-1
scaledZ = (dctZ-min(dctZ(:)))./(max(dctZ(:)-min(dctZ(:))));
min(scaledZ(:));                        % the min is 0
max(scaledZ(:));                        % the max 1

Z3d = reshape(scaledZ, [8,8,(m*n)/64]); % Create 3D array

%% Quantizer
qbits = 16;                     % Number of Quantizer bits

if qbits == 8
   Zt=im2uint8(Z3d);                    % Quantize to 2^8 levels
elseif qbits == 16
   Zt=im2uint16(Z3d);                   % Quantize to 2^16 levels
end       

%% Conversion to bit-stream:
% we transmit N DCT blocks at a time
% Convert N DCT blocks into a single column vector
% convert each element in this column vector into an 8-bit binary number using de2bi, stored as a row
%     - ex [5] -> [1 0 1 0 0 0 0 0]
 
% ex start with N=2 blocks of 8x8 DCT blocks
 
% Assumptions: dct is a bunch of 8x8 dct blocks stacked on top of each other
% and that the values in dct are quantized to integer values
 
% get dimensions of dct
[M,N,B] = size(Zt);
 
len = M*N;
bs = zeros(M*N*B,1); %preallocate memory for vector
for i = 1:B
    bs_temp = reshape(Zt(:,:,i),len,1); % reshape B MxN dct blocks into a long column vector
    start = len*(i-1)+1;
    stop  = len*i;
    bs(start:stop) = bs_temp;  % put values into correct location in bs
end
 
bs = de2bi(bs,qbits); % convert quantized values into 8-bit numbers. right MSB; use correct flag for left MSB

%% Modulation
T = 32; % number of samples per pulse
K = 4; % number of symbols overlap per pulse

% half-sine pulse
t1 = linspace(0,T); 
sPulse = sin(pi*t1/T); 

% SRRC Pulse
A = norm(sPulse,2);
a = 0.5; % roll-off factor
t2 = linspace(-K*T,K*T,2*K*(T)+1);
srrcPulse = A*rcosdesign(a,K,2*T,'sqrt');

PAM_level = 2;
R = qbits/(qbits/log2(PAM_level)); % number of quantization bits

% calculate amplitudes for all pulses
Am = zeros(M*N*B,R);
for i = 1:M*N*B
    for j = 1:R % perform 16-PAM
        head = 1+4*(j-1);
        tail = 4*j;
        Am(i,j) = bi2de(bs(i,head:tail));
    end
end
sp1 = Am(1,1)*sPulse; srrcp1 = Am(1,1)*srrcPulse;
plot(t1,sp1,t2,srrcp1)





% TEST = zeros(1,64*5+1);
% TEST(1:257) = TEST(1:257) + Am(1,1)*srrcPulse;
% TEST(65:321) = TEST(65:321) + Am(1,2)*srrcPulse
% plot(TEST)
% hold on
% plot(1:257,Am(1,1)*srrcPulse)
% hold on
% plot(65:321,Am(1,2)*srrcPulse)



