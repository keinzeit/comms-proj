%% final test
I0 = imread('file.jpeg');                   % Import Image
I0 = imrotate(I0, 270);                     % Rotate Image
[Zt dctZ I] = ImPreProc(I0);                % Send image for processing

%figure(1); imshow(I)                       % Show Processed Image        
%figure(2); imshow(dctZ)                    % Show DCT of Image

qBits = 16;                                 % Quantization bits
Zq = quantizer(qBits, Zt);                  % Send DCT image to be quantized


% Create Bit Stream out of Quantized Data
[bitStream M N B] = Convert_to_Bitstream(qBits, Zq);

T = 32;                                     % Sampling Period
PAM_level = 2;                              % PAM levels

% Modulate Bit Stream using Half Sine Pulses
[modBsSP sPulse Am] = ModulatedSinePulse(bitStream, T, qBits, PAM_level, M, N, B);



alpha = 0.5;                                % Roll- Off Factor
K = 4;

% Modulate Bit Stream using SRRC Pulses
[SRRCPulse modStreamSRRC] = ModulatedSRRCPulse(bitStream, alpha, T, K, Am, M, N, B, qBits, PAM_level);

% pulse time and frequency plots
figure(3); plot(sPulse)                     %Plots single Half-sine pulse
figure(4); freqz(sPulse)
figure(5); plot(SRRCPulse)                  %Plots single SRRC pulse                      
figure(6); freqz(SRRCPulse)
figure(7); plot(modBsSP)                    %Modulated signal with HS pulse
figure(8); freqz(modBsSP)
figure(9); plot(modStreamSRRC)              %Modulated Signal with SRRC pulse
figure(10); freqz(modStreamSRRC)

%eye diagrams for the output of the modulator
eyediagram(modBsSP, 32 , 1, 16)
eyediagram(modStreamSRRC, 32, 1)

noise = 0.1;                                %noise to be added after channel

[SRRC_Channel_out HS_Channel_out Noisy_HS Noisy_SRRC] = channel(modStreamSRRC, modBsSP, noise);

%eye diagrams of Modulated signals at the output of the Channel
eyediagram(SRRC_Channel_out, 32)
eyediagram(HS_Channel_out, 32, 1, 16)

%eye diagrams of Modulated signals with noise out of the Channel
eyediagram(Noisy_HS, 32, 1, 16)
eyediagram(Noisy_SRRC, 32)

%Pass Noisy Signals through the Matched Filter
[HS_MF_Out SRRC_MF_Out] = Matched_Filter(Noisy_HS, Noisy_SRRC, SRRCPulse, T);

eyediagram(HS_MF_Out, 32, 1, 16)
eyedigram(SRRC_MF_Out, 32)


