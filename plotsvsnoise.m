%% Recovered Images After Different levels of noise

% run this script to observe output images without drawing eyediagrams
% please first run finalproj_103.m to get the modulated signal
% relevant variables
modStreamHS;
modStreamSRRC;
tEye;
noise = 0.5; % this one can be changed

% channel
[SRRC_Channel_out,HS_Channel_out,Noisy_HS,Noisy_SRRC] = channel(modStreamSRRC, modStreamHS, noise);

% matched filter out
[HS_MF_Out,SRRC_MF_Out] = Matched_Filter(Noisy_HS, Noisy_SRRC, SRRC_Pulse, spb);

% ZF
[HS_ZF_Equalizer_Out] = ZF_equalizer(HS_MF_Out);
[SRRC_ZF_Equalizer_Out] = ZF_equalizer(SRRC_MF_Out);

% MMSE
channel_h = [1 1/2 3/4 -2/7];
channel_up = upsample(channel_h,spb);
N = 2.^(nextpow2(length(HS_MF_Out)));
% channel_FFT = fft(channel_h,N);

Qmmse = qmmse(channel_up,N,noise);
Qmmse0 = qmmse(channel_h,N,noise); % used to plot 1 period of the equalizer frequency response
% figure; freqz(ifft(Qmmse0)) 
% title('Frequency Response of the MMSE Equalizer (with noise=0.05)')
% figure; impz(Qmmse0)
% text = sprintf('Impulse Response of the MMSE Equalizer (noise = %1.2f)',noise);
% title(text)

% use the equalizer we just created and pass the matched filter output
% through it
HS_MMSE_Out = MMSE_Equalizer(Qmmse, HS_MF_Out);
SRRC_MMSE_Out = MMSE_Equalizer(Qmmse, SRRC_MF_Out);

%% Sampling and Detection - HS ZF Equalizer 
 
numBits; % this was found from way earlier
HS_ZF_sampledSignal = zeros(1,numBits);
currBit = 1;

% HS
% sample
while (currBit<=numBits)
   HS_ZF_sampledSignal(currBit) = HS_ZF_Equalizer_Out(currBit*spb);
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
text = sprintf('Recovered Image for Half-Sine ZF Equalizer (noise = %1.3f)',noise);
title(text)

%% Sampling and Detection - SRRC ZF Equalizer
 
numBits; % this was found from way earlier
SRRC_ZF_sampledSignal = zeros(1,numBits);
SRRC_ZF_decidedSignal = zeros(1,numBits);
currBit = 1;

% SRRC
% sample
while (currBit<=numBits)
   SRRC_ZF_sampledSignal(currBit) = SRRC_ZF_Equalizer_Out(currBit*spb+224);
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
text = sprintf('Recovered Image for SRRC ZF Equalizer (noise = %1.3f)',noise);
title(text)

%% Sampling and Detection - HS MMSE Equalizer 
 
numBits; % this was found from way earlier
HS_MMSE_sampledSignal = zeros(1,numBits);
HS_MMSE_decidedSignal = zeros(1,numBits);
currBit = 1;

% HS
% sample
while (currBit<=numBits)
   HS_MMSE_sampledSignal(currBit) = HS_MMSE_Out(currBit*spb);
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
text = sprintf('Recovered Image for Half-Sine MMSE Equalizer (noise = %1.3f)',noise);
title(text)

%% Sampling and Detection - SRRC MMSE Equalizer 
 
numBits; % this was found from way earlier
SRRC_MMSE_sampledSignal = zeros(1,numBits);
SRRC_MMSE_decidedSignal = zeros(1,numBits);
currBit = 1;

% SRRC
% sample
while (currBit<=numBits)
   SRRC_MMSE_sampledSignal(currBit) = SRRC_MMSE_Out(currBit*spb+224);
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
text = sprintf('Recovered Image for SRRC MMSE Equalizer (noise = %1.3f)',noise);
title(text)

%% Calculate SNRs

SNR_HS_ZF = snr(HS_ZF_Equalizer_Out,noise*randn(1,length(HS_ZF_Equalizer_Out)))
SNR_SRRC_ZF = snr(SRRC_ZF_Equalizer_Out,noise*randn(1,length(SRRC_ZF_Equalizer_Out)))

SNR_HS_MMSE = snr(HS_MMSE_Out,noise*randn(1,length(HS_MMSE_Out)))
SNR_SRRC_MMSE = snr(SRRC_MMSE_Out,noise*randn(1,length(SRRC_MMSE_Out)))

