%% final test
% Jon Manninen, Kenny Yau

%%%
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
%     - Conversion to Image
%     - Image post-processing

%% Image Pre-processing
I0 = imread('file.jpeg');                   % Import Image
I0 = imrotate(I0, 270);                     % Rotate Image
[Zt,dctZ,I] = ImPreProc(I0);                % Send image for processing

figure(1); imshow(I); title('Original Image')                       % Show Processed Image        
figure(2); imshow(dctZ); title('DCT of Image')                    % Show DCT of Image

%% Conversion to Bit Stream
qBits = 16;                                 % Quantization bits
Zq = quantizer(qBits,Zt);                   % Send DCT image to be quantized

% Create Bit Stream out of Quantized Data
[bitStream,M,N,B] = Convert_to_Bitstream(qBits, Zq);

T = 32;                                     % Number of samples per bit duration
PAM_level = 2;                              % PAM levels

%% Modulation 
% Modulate Bit Stream using Half Sine Pulses
[modStreamHS,HS_Pulse,Am] = ModulatedSinePulse(bitStream, T, qBits, PAM_level, M, N, B);

% Modulate Bit Stream using SRRC Pulses
alpha = 0.5;                                % Roll-Off Factor
K = 4;                                      % # of bit durations each side of the SRRC Pulse should last
[SRRC_Pulse,modStreamSRRC] = ModulatedSRRCPulse(bitStream,alpha,T,K,Am,M,N,B,qBits,PAM_level);

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

%eye diagrams for the output of the modulator
eyediagram(modStreamHS,T ,1,16); title('Eye Diagram of Modulated Signal - Half-Sine Pulse')
eyediagram(modStreamSRRC,T,1); title('Eye Diagram of Modulated Signal - SRRC Pulse')

%% Channel Effects
noise = 0.05;                                % noise to be added after channel (power?)

[SRRC_Channel_out,HS_Channel_out,Noisy_HS,Noisy_SRRC] = channel(modStreamSRRC, modStreamHS, noise);

% eye diagrams of Modulated signals at the output of the Channel
eyediagram(SRRC_Channel_out, 32);title('Eye Diagram of SRRC Signal after Channel')
eyediagram(HS_Channel_out, 32, 1, 16); title('Eye Diagram of HS Signal after Channel')

% eye diagrams of Modulated signals with noise out of the Channel
eyediagram(Noisy_HS, 32, 1, 16); title('Eye Diagram of Noisey HS Signal')
eyediagram(Noisy_SRRC, 32); title('Eye Diagram of Noisey SRRC Signal')

%% Matched Filter
% Pass Noisy Signals through the Matched Filter
[HS_MF_Out,SRRC_MF_Out] = Matched_Filter(Noisy_HS, Noisy_SRRC, SRRC_Pulse, T);

eyediagram(HS_MF_Out, 32); title('Eye Diagram of HS Signal after Matched Filter')
eyediagram(SRRC_MF_Out, 32); title('Eye Diagram of SRRC Signal after Matched Filter')

%% ZF Equalizer
% Pass Signal Output from the Matched Filter through the Zero-Forcing
% Equalizer

% note: ZF_equalizer is now written to filter output from any kind of
% pulse-shaping function
[HS_ZF_Equalizer_Out] = ZF_equalizer(HS_MF_Out);
[SRRC_ZF_Equalizer_Out] = ZF_equalizer(SRRC_MF_Out);

eyediagram(HS_ZF_Equalizer_Out, 32); title('Eye Diagram Output of ZF Equalizer HS')
eyediagram(SRRC_ZF_Equalizer_Out, 32); title('Eye Diagram Output of ZF Equalizer SRRC')

%% MMSE Equalizer
% half-sine pulse
% take FFT of channel impulse response
channel_h = [1 1/2 3/4 -2/7];
channel_up = upsample(channel_h,T);
N = 2.^(nextpow2(length(HS_ZF_Equalizer_Out)));
channel_FFT = fft(channel_h,N);

% create MMSE Equalizer w/o noise
Qmmse = qmmse(channel_up,N);
Qmmse0 = qmmse(channel_h,N); % used to plot 1 period of the equalizer frequency response
figure; freqz(ifft(Qmmse0))
title('Frequency Response of the MMSE Equalizer (no noise)')

% % create MMSE Equalizer including noise
% Qmmse = qmmse(channel_up,N,noise);
% Qmmse0 = qmmse(channel,N); % used to plot 1 period of the equalizer frequency response
% figure; freqz(ifft(Qmmse0)) 
% title('Frequency Response of the MMSE Equalizer (with noise=0.05)')

% use the equalizer we just created and pass the matched filter output
% through it

MMSE_HS_Out = MMSE_Equalizer(Qmmse, HS_MF_Out);
MMSE_SRRC_Out = MMSE_Equalizer(Qmmse, SRRC_MF_Out);


%Eye diagrams of SRRC and HS from the MMSE Equalizer

eyediagram(MMSE_HS_Out, 32); title('HS Output of MMSE Equalizer')
eyediagram(MMSE_SRRC_Out, 32); title('SRRC Output of MMSE Equalizer')







