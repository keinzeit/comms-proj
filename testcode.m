%% Test Code File

%% Applying MMSE Equalizer
% FFTs of MF Outputs
HS_MF_Out_FFT = fftshift(fft(HS_MF_Out,N));
figure(399); freqz(HS_MF_Out);
figure(401); freqz(ifft(HS_MF_Out_FFT));

HS_MMSE_Out_FFT = HS_MF_Out_FFT.*Qmmse;
HS_MMSE_Out = fftshift(ifft(HS_MMSE_Out_FFT));
figure(402); freqz(HS_MMSE_Out)


%% Calculating Qzf via FFT
% channel
h = [1 1/2 3/4 -2/7];
h_up = upsample(h,32); % upsampling h inserts 32-1 zeros between each element in h

% find frequency response of the channel using fft
% Output is centered around 0Hz
H = fftshift(fft(h,512)); 
% H_up = fftshift(fft(h_up,512));
figure(301)
plot(abs(H))
figure(302)
plot(abs(fftshift(H)))

Qzf = 1./H;
zf = ifft(Qzf);
figure(202), plot(zf)
title('Impulse Response of ZF Equalizer - HS')
figure(203), freqz(zf)
title('Frequency Response of ZF Equalizer - HS')

H_Qzf = H.*Qzf;
h_zf = ifft(H_Qzf);
figure(204); freqz(h_zf)


%% Plotting the signal at various stages in the communication chain all on the same figure

figure(101)
plot(modStreamHS); hold on
plot(HS_Channel_out);
plot(Noisy_HS);
plot(HS_MF_Out);
plot(HS_ZF_Equalizer_Out)
title('Signal At Various Points in the Comms Chain - HS')
legend('modStreamHS','HS\_Channel\_out','Noise\_HS','HS\_MF\_Out','ZF\_Equalizer\_Out\_HS')

figure(102)
plot(modStreamSRRC); hold on
plot(SRRC_Channel_out);
plot(Noisy_SRRC);
plot(SRRC_MF_Out);
plot(SRRC_ZF_Equalizer_Out)
title('Signal At Various Points in the Comms Chain - SRRC')
legend('modStreamSRRC','SRRC\_Channel\_out','Noise\_SRRC','SRRC\_MF\_Out','ZF\_Equalizer\_Out\_SRRC')