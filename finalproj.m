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
I0 = imread('file.jpeg');               % Import Image
I0 = imrotate(I0, 270);                 % Rotate Image
I = I0(801:864,401:464,:);            % Crop Image
I = imresize(I, 0.25);                  % Resize image

% figure(1)
% imshow(I)                               % Display original image

Z = rgb2gray(I);                        % Convert to grayscale
Z = im2double(Z);                       % Convert image to type double

m = size(Z,1);                          % Get number of rows
n = size(Z,2);                          % number of columns

imshow(Z)
%fun = @(block_struct)(block_struct.data);
fun = @dct2;
dctZ= blkproc(Z,[8,8],fun);             %Perform DCT on image

%figure(2)
%imshow(dctZ)                            %Display the DCT of image

% scale the DCT image matrix from 0-1
minZ = min(dctZ(:)); maxZ = max(dctZ(:));
scaledZ = (dctZ-min(dctZ(:)))./(max(dctZ(:)-min(dctZ(:))));
% min(scaledZ(:));                        % the min is 0
% max(scaledZ(:));                        % the max 1

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
bsd = zeros(M*N*B,1); %preallocate memory for vector
for i = 1:B
    bs_temp = reshape(Zt(:,:,i),len,1); % reshape B MxN dct blocks into a long column vector
    start = len*(i-1)+1;
    stop  = len*i;
    bsd(start:stop) = bs_temp;  % put values into correct location in bs
end
 
bs = de2bi(bsd,qbits); % convert quantized values into 8-bit numbers. right MSB; use correct flag for left MSB
% bs = [1 0 1 1 1 0 0 0 0 1 1 0 0 0 1 1
%       0 1 1 1 0 0 0 0 1 1 0 0 0 1 0 0
%       0 1 1 0 0 0 0 0 1 0 0 1 1 0 1 0
%       0 1 1 1 0 1 1 0 0 1 1 1 0 1 0 1];

%% Modulation
T = 32; % technically the bit duration, but for now it is number of samples per bit duration
K = 4; % half the number of bit durations per SRRC pulse

% half-sine pulse
t1 = linspace(0,T,T+1); % I have pulses include the first value of the next pulse
sPulse = sin(pi*t1/T); 

% SRRC Pulse
A = norm(sPulse,2);
a = 0.5; % roll-off factor
t2 = linspace(-K*T,K*T,2*K*T+1);
srrcPulse = A*rcosdesign(a,K,2*T,'sqrt'); % energy of rcosdesign is one, just multiply by energy of half-sine pulse

PAM_level = 2;
R = log2(PAM_level); % number of bits used to change amplitude of pulse shape
pulsesPerRow = qbits/R; % number of pulses made from a row of Am

% % calculate amplitudes for all pulses
% Am = zeros(M*N*B,R);
% for i = 1:M*N*B
%     for j = 1:pulsesPerRow % perform 16-PAM
%         head = 1+log2(PAM_level)*(j-1);
%         tail = log2(PAM_level)*j;
%         Am(i,j) = bi2de(bs(i,head:tail));
%     end
% end
% if PAM_level == 2
%     Am(Am==0) = -1;
% end

% % modulate bk's with half-sine pulses
% modBsSP = zeros(1,size(Am,1)*size(Am,2)*T+1);
% tModBsSP = (1:length(modBsSP))/T;
% for i = 1:M*N*B*qbits
%     j = fix((i-1)/pulsesPerRow)+1;
%     k = mod(i-1,pulsesPerRow)+1;
%     start = (i-1)*T+1;
%     stop  = i*T+1;
%     modBsSP(start:stop) = Am(j,k)*sPulse; %consider removing modBsSP for half-sine pulse  
% end

% calculate amplitudes for all pulses
Am = zeros(M*N*B,R);
for i = 1:M*N*B
    for j = 1:pulsesPerRow % perform 16-PAM
        head = 1+log2(PAM_level)*(j-1);
        tail = log2(PAM_level)*j;
        Am(i,j) = bi2de(bs(i,head:tail));
    end
end
if PAM_level == 2
    Am(Am==0) = -1;
end

modBsSP = zeros(1,size(Am,1)*size(Am,2)*T+1); 
tModBsSP = (1:length(modBsSP))/T;
for i = 1:M*N*B*qbits
    j = fix((i-1)/pulsesPerRow)+1;
    k = mod(i-1,pulsesPerRow)+1;
    start = (i-1)*T+1;
    stop  = i*T+1;
    modBsSP(start:stop) = Am(j,k)*sPulse; %consider removing modBsSP for half-sine pulse
    
end

% proof-of-concept: plot the first 16 symbols
figure(11)
plot(tModBsSP(1:T*16),modBsSP(1:T*16))
title('First 16 Symbols of Modulation w/ Half-Sine Pulse')
xlabel('Time (s)')
ylabel('Amplitude')

% modulate bk's with srrc pulses
modStreamSRRC = zeros(1,(size(Am,1)*size(Am,2)+2*K-1)*T+1);
% tModBsSP = (1:length(modStreamSRRC))/T;
for i = 1:M*N*B*qbits
    j = fix((i-1)/pulsesPerRow)+1;
    k = mod(i-1,pulsesPerRow)+1;
    start = (i-1)*T+1;
    stop  = (i-1)*T+2*K*T+1;
    modStreamSRRC(start:stop) = modStreamSRRC(start:stop) + Am(j,k)*srrcPulse; %consider removing modBsSP for half-sine pulse  
end

whos modBsSP

%% Channel
delay = zeros(1, 31);
h = [1 delay 1/2 delay 3/4 delay -2/7];

y = conv(modBsSP, h);

freqz(h)
stem(h)

noise = rand(1,length(y));
yn = y+ noise;

%% Matched filter
g = sPulse;







