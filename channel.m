function [SRRC_Channel_out HS_Channel_out Noisy_HS Noisy_SRRC] = Channel(modStreamSRRC, modBsSP, noise)

%Function imports the modulated bit stream and passes them through the
%channel for transmission.  The cannel is created, plotted, and modulated
%signals are passed through.  Noise is also added in this function

%inputs modStreamSRRC modBsSP, noise
%outputs HS_Channel_out, SRRC_Channel_out, Noisy_HS, Noisy_SRRC

delay = zeros(1, 31);                           %vector of 31 zeros to space between values in h
h = [1 delay 1/2 delay 3/4 delay -2/7];         % channel filter FIR

%Plot Impulse and Frequency Response of Channel
figure(101); freqz(h); suptitle('Frequency Response of the Channel')                                       
figure(102); stem(h); suptitle('Impulse Response of Channel')   

%Create Outputs by Convolving Signals with the Channel
HS_Channel_out = conv(modBsSP, h);                           
SRRC_Channel_out = conv(modStreamSRRC, h);   

%% Noise 
sigma = sqrt(noise);

sineNoise = sigma*randn(1,length(HS_Channel_out));       %noise
srrcNoise = sigma*randn(1, length(SRRC_Channel_out));      

Noisy_HS = HS_Channel_out+ sineNoise;                    %signals with noise
Noisy_SRRC = SRRC_Channel_out+ srrcNoise;
