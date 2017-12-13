%% Plotting the signal at various stages in the communication chain all on the same figure

% run this script to see what the signal looks like at various points in
% the communication chain
% please first run finalproj_103.m to get the modulated signal

figure(101)
plot(modStreamHS(tEye)); hold on
plot(HS_Channel_out(tEye));
plot(Noisy_HS(tEye));
plot(HS_MF_Out(tEye));
plot(HS_ZF_Equalizer_Out(tEye))
plot(HS_MMSE_Out(tEye))
title('Signal At Various Points in the Comms Chain - HS')
legend('modStreamHS','HS\_Channel\_out','Noise\_HS','HS\_MF\_Out','ZF\_Equalizer\_Out\_HS','HS\_MMSE\_Out')

figure(102)
plot(modStreamSRRC(tEye)); hold on
plot(SRRC_Channel_out(tEye));
plot(Noisy_SRRC(tEye));
plot(SRRC_MF_Out(tEye));
plot(SRRC_ZF_Equalizer_Out(tEye))
plot(SRRC_MMSE_Out(tEye))
title('Signal At Various Points in the Comms Chain - SRRC')
legend('modStreamSRRC','SRRC\_Channel\_out','Noise\_SRRC','SRRC\_MF\_Out','ZF\_Equalizer\_Out\_SRRC','SRRC\_MMSE\_Out')