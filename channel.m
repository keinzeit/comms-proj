function [ySRRC ySine] = Channel(modBsSRRC, modBsSP, noise)

%Function imports the modulated bit stream dna passes them through the
%channel for transmission.  The channel is also created in this function
%and noise is added to the signals at the end.  This function also plots
%eye diagrams of both signals before and after the noise is added

%inputs modBsSRRC and modBsSP
%outputs ySRRC and ySine with added noise

delay = zeros(1, 31);                           %vector of 31 zeros to space between values in h
h = [1 delay 1/2 delay 3/4 delay -2/7];         % channel filter FIR

figure
freqz(h)                                        %Plots frequency response of h

figure
stem(h)                                         %plots step response of h

ysine = conv(modBsSP, h);                       %channel output 

ysrrc = conv(modStreamSRRC, h);                 %channel output

figure
plot(ysine)

figure
plot(ysrrc)

%% Noise 
sigma = sqrt(noise);

sineNoise = sigma*randn(1,length(ysine));       %noise
srrcNoise = sigma*randn(1, length(ysrrc));      

ysn = ysine+ sineNoise;                         % signals with noise
ysrrcn = ysrrc+ srrcNoise;
