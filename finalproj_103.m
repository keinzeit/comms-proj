%% EE107 Final Project
% Jon Manninen, Kenny Yau

%%%
% This project is a simulation of the model of communication systems for
% images, such as that used for digital HDTV. It consists of a number of
% parts:
%     - Image pre-processing (ImagePreProcess_gray.m)
%     - Conversion to bit-stream
%     - Modulation
%     - Channel Effects
%     - Noise
%     - Matched Filter
%     - Receiver Equalization
%     - Sampling and Detection
%     - Conversion to Image
%     - Image post-processing

%% Image Pre-processing
qBits = 16;
[Zq,blockrow,blockcol,imrow,imcol,minZ,maxZ] = ImagePreProcess_gray('file.jpeg',qBits);

% figure(11); imshow(I); title('Original Image')                       % Show Processed Image        
% figure(12); imshow(dctZ); title('DCT of Image')                    % Show DCT of Image

%% Conversion to Bit Stream

% create bit stream out of quantized data
[bitstream,R,C,N] = blocks2bitStream(qBits,Zq);
% store number of total bits sent
numBits = length(bitstream);


%% Modulation 
PAM_level = 2; % PAM levels
spb = 32;      % Samples per bit

% Modulate Bit Stream using Half Sine Pulses
t1 = linspace(0,spb,spb+1);                 % Pulses include the first value of the next pulse
HS_Pulse = sin(pi*t1/spb);
modStreamHS = pamModulate(bitstream,HS_Pulse,spb);

% Modulation w/ SRRC Pulse
alpha = 0.5;
K = 4;
SRRC_Pulse = srrcPulse(alpha,spb,K);
SRRC_Pulse = normalizePulse(SRRC_Pulse,HS_Pulse);
modStreamSRRC = pamModulate(bitstream,SRRC_Pulse,spb);

% pulse time and frequency plots
figure(3); plot(HS_Pulse); title('Half-Sine Pulse');
figure(4); freqz(HS_Pulse); title('Frequency Response of Half-Sine Pulse')
figure(5); plot(SRRC_Pulse); title('Square-Root Raised-Cosine Pulse')
figure(6); freqz(SRRC_Pulse); title('Frequency Response of SRRC Pulse')

% modulated signal plots
figure(7); freqz(modStreamHS); title('Frequency Response of Modulated Signal - HS Pulse')
figure(8); freqz(modStreamSRRC); title('Frequency Response of Modulated Signal - SRRC Pulse')

% random section of 10 bit-durations of each modulated signal
figure(9); plot(modStreamHS((1:10*T)+100*T)); title('Modulated Half-Sine Pulse')
figure(10); plot(modStreamSRRC((1:10*T)+100*T)); title('Modulated SRRC Pulse')

tEye = (2^10+1):2^14; % portion of signal we use to draw eye diagram
% eye diagrams for the output of the modulator
eyediagram(modStreamHS(tEye),T ,1,16); title('Eye Diagram of Modulated Signal - Half-Sine Pulse')
eyediagram(modStreamSRRC(tEye),T,1); title('Eye Diagram of Modulated Signal - SRRC Pulse')

%% Channel Effects
noise = 0.05; % noise to be added after channel (power?)

[SRRC_Channel_out,HS_Channel_out,Noisy_HS,Noisy_SRRC] = channel(modStreamSRRC, modStreamHS, noise);

% eye diagrams of Modulated signals at the output of the Channel
eyediagram(SRRC_Channel_out(tEye), 32);title('Eye Diagram of SRRC Signal after Channel')
eyediagram(HS_Channel_out(tEye), 32, 1, 16); title('Eye Diagram of HS Signal after Channel')

% eye diagrams of Modulated signals with noise out of the Channel
eyediagram(Noisy_HS(tEye), 32, 1, 16); title('Eye Diagram of Noisey HS Signal')
eyediagram(Noisy_SRRC(tEye), 32); title('Eye Diagram of Noisey SRRC Signal')

%% Matched Filter
% Pass Noisy Signals through the Matched Filter
[HS_MF_Out,SRRC_MF_Out] = Matched_Filter(Noisy_HS, Noisy_SRRC, SRRC_Pulse, T);

eyediagram(HS_MF_Out(tEye), 32); title('Eye Diagram of HS Signal after Matched Filter')
eyediagram(SRRC_MF_Out(tEye), 32); title('Eye Diagram of SRRC Signal after Matched Filter')

%% ZF Equalizer
% Pass Signal Output from the Matched Filter through the Zero-Forcing
% Equalizer

% note: ZF_equalizer is now written to filter output from any kind of
% pulse-shaping function
[HS_ZF_Equalizer_Out] = ZF_equalizer(HS_MF_Out);
[SRRC_ZF_Equalizer_Out] = ZF_equalizer(SRRC_MF_Out);

eyediagram(HS_ZF_Equalizer_Out(tEye), 32); title('Eye Diagram Output of ZF Equalizer HS')
eyediagram(SRRC_ZF_Equalizer_Out(tEye), 32); title('Eye Diagram Output of ZF Equalizer SRRC')

%% MMSE Equalizer

% take FFT of channel impulse response
channel_h = [1 1/2 3/4 -2/7];
channel_up = upsample(channel_h,T);
N = 2.^(nextpow2(length(HS_MF_Out)));
channel_FFT = fft(channel_h,N);

% create MMSE Equalizer w/o noise
Qmmse = qmmse(channel_up,N);
Qmmse0 = qmmse(channel_h,N); % used to plot 1 period of the equalizer frequency response
figure; freqz(ifft(Qmmse0))
title('Frequency Response of the MMSE Equalizer (no noise)')

% use the equalizer we just created and pass the matched filter output
% through it
HS_MMSE_Out = MMSE_Equalizer(Qmmse, HS_MF_Out);
SRRC_MMSE_Out = MMSE_Equalizer(Qmmse, SRRC_MF_Out);


%Eye diagrams of SRRC and HS from the MMSE Equalizer
eyediagram(HS_MMSE_Out(tEye), 32); title('HS Output of MMSE Equalizer')
eyediagram(SRRC_MMSE_Out(tEye), 32); title('SRRC Output of MMSE Equalizer')

%% MMSE Equalizer (with noise)

% create MMSE Equalizer including noise
Qmmse = qmmse(channel_up,N,noise);
Qmmse0 = qmmse(channel_h,N,noise); % used to plot 1 period of the equalizer frequency response
figure; freqz(ifft(Qmmse0)) 
title('Frequency Response of the MMSE Equalizer (with noise=0.05)')

% use the equalizer we just created and pass the matched filter output
% through it
HS_MMSE_Out = MMSE_Equalizer(Qmmse, HS_MF_Out);
SRRC_MMSE_Out = MMSE_Equalizer(Qmmse, SRRC_MF_Out);


%Eye diagrams of SRRC and HS from the MMSE Equalizer
eyediagram(HS_MMSE_Out(tEye), 32); title('HS Output of MMSE Equalizer (with noise)')
eyediagram(SRRC_MMSE_Out(tEye), 32); title('SRRC Output of MMSE Equalizer (with noise)')


%% Sampling and Detection - HS ZF Equalizer 
 
numBits; % this was found from way earlier
HS_ZF_sampledSignal = zeros(1,numBits);
currBit = 1;

% HS
% sample
while (currBit<=numBits)
   HS_ZF_sampledSignal(currBit) = HS_ZF_Equalizer_Out(currBit*T);
   currBit = currBit + 1;
end

% decision
HS_ZF_decidedSignal(HS_ZF_sampledSignal>0) = 1;
HS_ZF_decidedSignal(HS_ZF_sampledSignal<=0) = 0;

% Conversion to Image
newZq_HS_ZF = bitstream2blocks(HS_ZF_decidedSignal,qBits,R,C);
% Image Post-processing
newZ_HS_ZF = ImagePostProcess_gray(newZq_HS_ZF,blockrow,blockcol,imrow,imcol,minZ,maxZ);
figure
imshow(newZ_HS_ZF)
text = sprintf('Recovered Image for Half-Sine ZF Equalizer (noise = %1.2f)',noise);
title(text)

%% Sampling and Detection - SRRC ZF Equalizer
 
numBits; % this was found from way earlier
SRRC_ZF_sampledSignal = zeros(1,numBits);
SRRC_ZF_decidedSignal = zeros(1,numBits);
currBit = 1;

% SRRC
% sample
while (currBit<=numBits)
   SRRC_ZF_sampledSignal(currBit) = SRRC_ZF_Equalizer_Out(currBit*T+224);
   currBit = currBit + 1;
end

% decision
SRRC_ZF_decidedSignal(SRRC_ZF_sampledSignal>0) = 1;
SRRC_ZF_decidedSignal(SRRC_ZF_sampledSignal<=0) = 0;

% Conversion to Image
newZq_SRRC_ZF = bitstream2blocks(SRRC_ZF_decidedSignal,qBits,R,C);
% Image Post-processing
newZ_SRRC_ZF = ImagePostProcess_gray(newZq_SRRC_ZF,blockrow,blockcol,imrow,imcol,minZ,maxZ);
figure
imshow(newZ_SRRC_ZF)
text = sprintf('Recovered Image for SRRC ZF Equalizer (noise = %1.2f)',noise);
title(text)

%% Sampling and Detection - HS MMSE Equalizer 
 
numBits; % this was found from way earlier
HS_MMSE_sampledSignal = zeros(1,numBits);
HS_MMSE_decidedSignal = zeros(1,numBits);
currBit = 1;

% HS
% sample
while (currBit<=numBits)
   HS_MMSE_sampledSignal(currBit) = HS_MMSE_Out(currBit*T);
   currBit = currBit + 1;
end

% decision
HS_MMSE_decidedSignal(HS_MMSE_sampledSignal>0) = 1;
HS_MMSE_decidedSignal(HS_MMSE_sampledSignal<=0) = 0;

% Conversion to Image
newZq_HS_MMSE = bitstream2blocks(HS_MMSE_decidedSignal,qBits,R,C);
% Image Post-processing
newZ_HS_MMSE = ImagePostProcess_gray(newZq_HS_MMSE,blockrow,blockcol,imrow,imcol,minZ,maxZ);
figure
imshow(newZ_HS_MMSE)
text = sprintf('Recovered Image for Half-Sine MMSE Equalizer (noise = %1.2f)',noise);
title(text)

%% Sampling and Detection - SRRC MMSE Equalizer 
 
numBits; % this was found from way earlier
SRRC_MMSE_sampledSignal = zeros(1,numBits);
SRRC_MMSE_decidedSignal = zeros(1,numBits);
currBit = 1;

% SRRC
% sample
while (currBit<=numBits)
   SRRC_MMSE_sampledSignal(currBit) = SRRC_MMSE_Out(currBit*T+224);
   currBit = currBit + 1;
end

% decision
SRRC_MMSE_decidedSignal(SRRC_MMSE_sampledSignal>0) = 1;
SRRC_MMSE_decidedSignal(SRRC_MMSE_sampledSignal<=0) = 0;

% Conversion to Image
newZq_SRRC_MMSE = bitstream2blocks(SRRC_MMSE_decidedSignal,qBits,R,C);
% Image Post-processing
newZ_SRRC_MMSE = ImagePostProcess_gray(newZq_SRRC_MMSE,blockrow,blockcol,imrow,imcol,minZ,maxZ);
figure
imshow(newZ_SRRC_MMSE)
text = sprintf('Recovered Image for SRRC MMSE Equalizer (noise = %1.2f)',noise);
title(text)
