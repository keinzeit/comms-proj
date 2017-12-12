function [MMSE_Out] = MMSE_Equalizer(Qmmse,MF_Out)

%Function to pass the ourputs of the Matched filter through the MMSE Equalizer created in it own function
%This function will undo the effects of the channel and will output a data
%stream ready to be sampled

%inputs Qmmse, HS_MF_Out and SRRC_MF_Out
%Outputs MMSE_HS_Out and MMSE_SRRC_Out

% Transform Matched FIlter Outputs to Frequency Domain

LS = 2^(nextpow2(length(MF_Out))); 
S_FFT = fft(MF_Out, LS);

QMMSE_Out = (S_FFT.*Qmmse);

MMSE_Out = ifft(QMMSE_Out);


return