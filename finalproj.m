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

%% Image Pre-processing
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

%figure(2)
%imshow(dctZ)                            %Display the DCT of image

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
 
% bs = de2bi(bs,qbits); % convert quantized values into 8-bit numbers. right MSB; use correct flag for left MSB
bs = [1 0 1 1 1 0 0 0 0 1 1 0 0 0 1 1
      0 1 1 1 0 0 0 0 1 1 0 0 0 1 0 0
      0 1 1 0 0 0 0 0 1 0 0 1 1 0 1 0
      0 1 1 1 0 1 1 0 0 1 1 1 0 1 0 1];

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
for i = 1:4% M*N*B
    for j = 1:qbits/R % perform 16-PAM
        head = 1+log2(PAM_level)*(j-1);
        tail = log2(PAM_level)*j;
        Am(i,j) = bi2de(bs(i,head:tail));
    end
end
if PAM_level ==2
    Am(Am==0) = -1;
end

sp1 = Am(1,1)*sPulse; srrcp1 = Am(2,1)*srrcPulse;
plot(t1,sp1,t2,srrcp1)

signalSP = zeros(1,size(Am,1)*size(Am,2)*T);