function [ZF_Equalizer_Out_HS ZF_Equalizer_Out_SRRC] = ZF_equalizer(HS_MF_Out, SRRC_MF_Out)

%Function is a zero-forcing equalizer which undoes the effects of the
%Channel.  The Equalizer compensates for channel distortion by inverting
%the channel response.  This function also plots the impulse and frequency
%response of the equalizer

delay = zeros(1,31);                            %Matrix of zeros to delay between pulses
h = [1 delay 1/2 delay 3/4 delay -2/7];         %Coefficients of the Channel

L1 = length(HS_MF_Out); P1 = 2^nextpow2(L1);
L2 = length(SRRC_MF_Out); P2 = 2^nextpow2(L2);

SRRC_FFT = fft(SRRC_MF_Out, P2);                %FFT of the SRRC signal
HS_FFT = fft(HS_MF_Out, P1);                    %FFT of HS signal

Hf_HS = fft(h, length(HS_MF_Out));              %FFT of Channel for HS signal
Hf_SRRC = fft(h, length(SRRC_MF_Out));          %FFt of Channel for SRRC signal

ZFE_HS = 1./Hf_HS;                              %Equalizer for HS signal
ZFE_SRRC = 1./Hf_SRRC;                          %Equalizer for SRRC signal

freqz(ZFE_HS);suptitle('Frequency Response of Zero Forcing Equalizer')
freqz(ZFE_SRRC);suptitle('Frequency Response of Zero Forcing Equalizer')

ZF_Equalizer_Out_HS = ifft(Hf_HS.*ZFE_HS);
ZF_Equalizer_Out_SRRC = ifft(Hf_SRRC.*ZFE_SRRC);

figure(114),stem(ZF_Equalizer_Out_HS);
suptitle('Impulse Response of the Zero-Forcing Equalizer')

figure(115),stem(ZF_Equalizer_Out_SRRC);
suptitle('Impulse Response of the Zero-Forcing Equalizer')







